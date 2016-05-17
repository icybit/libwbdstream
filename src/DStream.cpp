#include "DStream.h"
#include "wbdstream/wbdstreamapi.h"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <unordered_map>
#include <vector>

#include "CharacteristicVector.h"
#include "Cluster.h"
#include "Common.h"
#include "Key.h"
#include "Parameters.h"

typedef std::unordered_map<Key, CharacteristicVector * > Gridlist;
typedef std::unordered_map<uint32_t, Cluster *> Clusters;


DSTREAM_PUBLIC void dstream_clusterize(uint8_t * buffer, uint32_t buffer_size)
{
	uint64_t time_now;
	Gridlist grid_list;
	Clusters clusters;
	float d_m, d_l;
	Deserialize(buffer, buffer_size, time_now, grid_list);
	ReassembleClusters(grid_list, clusters);
	CalculateDensityParams(d_m, d_l);
	AdjustClustering(grid_list, clusters, time_now, d_m, d_l);
	MergeChangesToBuffer(buffer, buffer_size, grid_list);

	for (auto it = grid_list.begin(); it != grid_list.end(); ++it)
	{
		delete it->second;
	}
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		delete it->second;
	}
	/*
	Clusters need to be pushed to key-value store;
	*/
}

DSTREAM_PUBLIC int dstream_calculate_gap_time() {
	/* Simplyfied:
	Original method counts the time for sparse grid to become dense. Basically, how many hits in a row it needs to become dense.
	In that case gap time is too small to recluster every time.
	Now we count only decay factor - time for dense grid to decay to sparse.*/
	int gap = 0;
	double dense_to_sparse = log(Parameters::c_l / Parameters::c_m) / log(Parameters::decay_factor);/* (11) - research paper*/
	gap = (int)floor(dense_to_sparse);
	return gap;
}

DSTREAM_PUBLIC double * dstream_calculate_xy_coords(double dx, double dy)
{
	double * coords = new double[2];
	CalculateXYCoords(dx, dy, coords[0], coords[1]);
	return coords;
}

DSTREAM_PUBLIC double * dstream_calculate_xy_distances(double x, double y)
{
	double * distances = new double[2];
	CalculateXYDistances(x, y, distances[0], distances[1]);
	return distances;
}

void AdjustClustering(Gridlist & grid_list, Clusters & clusters, uint64_t time_now, float d_m, float d_l)
{
	uint8_t status;
	GridTuple grid;
	CharacteristicVector * vect;
	Key key;
	Key key_neighbor;
	uint32_t label;

	RemoveSporadic(grid_list, time_now);
	UpdateDensities(grid_list, time_now, d_m, d_l);
	for (auto it = grid_list.begin(); it != grid_list.end(); ++it)
	{
		if (it->second->is_changed())
		{
			key = it->first;
			vect = it->second;
			status = vect->get_status();
			label = vect->get_label();
			/* SPARSE */
			if (status < TRANSITIONAL_FROM_SPARSE && status > SPORADIC)
			{
				if (label != NO_CLASS)
				{
					bool is_connected;
					if (clusters.count(label))
					{
						clusters[label]->RemoveElement(key.x_, key.y_);
					}
					vect->set_label(NO_CLASS);

					is_connected = CheckIfConnected(clusters, grid_list, label);
					/* If cluster was named after this grid, we have to rename it.*/
					if (vect->get_status() == SPARSE_FROM_DENSE && label == LabelHash(key) && is_connected)
					{
						RenameCluster(clusters, grid_list, label);
					}
				}
			}
			/* TRANSITIONAL*/
			else if (status < DENSE_FROM_SPARSE)
			{
				if (label == NO_CLASS)
				{
					CallClusteringOnGrid(key, grid_list, clusters);
				}
				else
				{
					bool renamed = true;
					if (status == TRANSITIONAL_FROM_DENSE)
					{
						/* if cluster was named after this grid we have to rename cluster - grid might be removed*/
						if (label == LabelHash(key))
						{
							renamed = RenameCluster(clusters, grid_list, label);
							label = vect->get_label();
						}
						/* If rename was not successful, that means there were no dense grids - cluster disassembled.*/
						if (renamed)
						{
							CheckNeighborsConnectivity(key, grid_list, clusters);
							if (!CheckForDenseNeighbors(key, grid_list))
							{
								clusters[label]->RemoveElement(key.x_, key.y_);
								vect->set_label(NO_CLASS);
							}
						}
					}
				}
			}
			/* DENSE*/
			else
			{
				/* If it doesn't belong to cluster, create its distinct cluster*/
				if (label == NO_CLASS)
				{
					label = LabelHash(key);
					if (clusters.count(label))
					{
						RenameCluster(clusters, grid_list, label);
					}
					clusters[label] = new Cluster(label);
					clusters[label]->AddElement(GridTuple(key.x_, key.y_, vect->get_status()));
					vect->set_label(label);
					CallClusteringOnGrid(key, grid_list, clusters);
				}
				/* If it belongs to cluster, update neighbors neighboring scores*/
				else
				{
					/* This whole branch might be unnecessary. If dense grid has label - it means, it was added to cluster. 
					 If it was added to cluster, CallClusteringOnGrid was already called on it.
					 DENSE_FROM_TRANSITIONAL && DENSE*/
					if (vect->get_status() > DENSE_FROM_SPARSE)
					{
						CallClusteringOnGrid(key, grid_list, clusters);
					}
					else
					{
						/* DENSE FROM SPARSE - logically, it should have NO_CLASS label;
						 If it has non-NO_CLASS label, it means, it was clustered by its neighbour.
						 So, no worries.*/
					}
				}
			}
			it->second->set_unchanged();
		}
	}
	RemoveEmptyClusters(clusters);
}

/* Calculate Dm and Dl;*/
void CalculateDensityParams(float & d_m, float & d_l)
{
	float denumerator = Parameters::total_grids * (1 - Parameters::decay_factor);
	d_m = Parameters::c_m / denumerator;
	d_l = Parameters::c_l / denumerator;
}

void CallClusteringOnGrid(Key grid, Gridlist & grid_list, Clusters & clusters)
{
	CharacteristicVector * vect;
	CharacteristicVector * vect_neighbor;
	Key key;
	Key key_neighbor;
	std::vector<Key> grids;
	Key neighbors[4];

	uint32_t label;
	uint32_t label_neighbor;

	grids.push_back(grid);
	for (uint32_t i = 0; i < grids.size(); i++)
	{
		key = Key(grids[i].x_, grids[i].y_);
		vect = grid_list[key];
		label = vect->get_label();
		GetNeighbors(key.x_, key.y_, neighbors);
		for (int j = 0; j < 4; j++)
		{
			key_neighbor = neighbors[j];
			if (grid_list.count(key_neighbor))
			{
				vect_neighbor = grid_list[key_neighbor];
				label_neighbor = vect_neighbor->get_label();
				if (label_neighbor != NO_CLASS && label != NO_CLASS && label_neighbor != label)
				{
					if (clusters[label_neighbor]->get_size() >= clusters[label]->get_size())
					{
						MergeClusters(*clusters[label_neighbor], *clusters[label], grid_list, clusters);
						label = vect->get_label();
					}
					else {
						MergeClusters(*clusters[label], *clusters[label_neighbor], grid_list, clusters);
					}
				}
				else if (vect->get_status() >= DENSE_FROM_SPARSE && label_neighbor != label && vect_neighbor->get_status() >= TRANSITIONAL_FROM_SPARSE)
				{
					clusters[label]->AddElement(GridTuple(key_neighbor.x_, key_neighbor.y_, vect_neighbor->get_status()));
					grids.push_back(key_neighbor);
					vect_neighbor->set_label(label);
				}
			}
		}

	}
}

void CalculateXYCoords(double dx, double dy, double & x, double & y)
{
	y = 360 * dy / (2 * M_PI * R_EARTH);
	double ry = y * M_PI / 180;
	double r = sin(M_PI_2 - ry) * R_EARTH;
	x = 360 * dx / (2 * M_PI * r);
}

void CalculateXYDistances(double x, double y, double & dx, double & dy)
{
	dy = (2 * M_PI * R_EARTH * y) / 360;
	double ry = y * M_PI / 180;
	double r = sin(M_PI_2 - ry) * R_EARTH;
	dx = (2 * M_PI * r * x) / 360;
}

/* Checks if cluster is connected. Splits it if not.*/
bool CheckIfConnected(Clusters & clusters, Gridlist & grid_list, uint32_t label)
{
	Key key;
	Key key_neighbor;
	Cluster * cluster = clusters[label];
	bool is_connected = false;

	std::vector<int> checklist(cluster->get_size(), 0);
	uint32_t counter = 0;

	for (uint32_t i = 0; i < checklist.size(); i++)
	{
		if (checklist[i] == 0)
		{
			counter++;
			key = Key(cluster->GetElement(i).x, cluster->GetElement(i).y);
			grid_list[key]->set_label(counter);
			checklist[i] = counter;

			RelabelNeighbors(key, grid_list, label, counter);
			/* can be split max into 4 clusters, no worries about cycle inside cycle*/
			for (uint32_t j = i + 1; j < cluster->get_size(); j++)
			{
				key = Key(cluster->GetElement(j).x, cluster->GetElement(j).y);
				if (!grid_list.count(key))
				{
					std::cout << "STOP" << std::endl;
				}
				if (grid_list[key]->get_label() == counter)
				{
					checklist[j] = counter;
				}
			}
		}
	}
	if (counter > 1)
	{
		SplitCluster(clusters, grid_list, checklist, label, counter);
		is_connected = false;
	}
	else
	{
		for (uint32_t i = 0; i < cluster->get_size(); i++)
		{
			key = Key(cluster->GetElement(i).x, cluster->GetElement(i).y);
			grid_list[key]->set_label(label);
		}
		is_connected = true;
	}
	return is_connected;
}

bool CheckForDenseNeighbors(Key & key, Gridlist & grid_list)
{
	Key neighbors[4];
	Key neighbor;
	uint32_t label;
	bool has_neighbor = false;

	label = grid_list[key]->get_label();
	GetNeighbors(key.x_, key.y_, neighbors);
	for (int i = 0; i < 4; i++)
	{
		neighbor = neighbors[i];
		if (grid_list.count(neighbor))
		{
			if (grid_list[neighbor]->get_status() >= DENSE_FROM_SPARSE && grid_list[neighbor]->get_label() == label)
			{
				has_neighbor = true;
				break;
			}
		}
	}
	return has_neighbor;
}

/* Checks if grid's neighbors are still connected to the cluster.
 If not - removes them from cluster and resets their labels.*/
void CheckNeighborsConnectivity(Key grid, Gridlist & grid_list, Clusters & clusters)
{
	Key key_neighbor;
	Key neighbors[4];
	CharacteristicVector * vect_neighbor;
	uint32_t label;

	label = grid_list[grid]->get_label();
	GetNeighbors(grid.x_, grid.y_, neighbors);
	for (int i = 0; i < 4; i++)
	{
		key_neighbor = neighbors[i];
		if (grid_list.count(key_neighbor))
		{
			vect_neighbor = grid_list[key_neighbor];
			if (vect_neighbor->get_label() == label && vect_neighbor->get_status() >= TRANSITIONAL_FROM_SPARSE
				&& vect_neighbor->get_status() < DENSE_FROM_SPARSE)
			{
				if (!CheckForDenseNeighbors(key_neighbor, grid_list))
				{
					clusters[label]->RemoveElement(key_neighbor.x_, key_neighbor.y_);
					vect_neighbor->set_label(NO_CLASS);
				}
			}
		}
	}
}

std::string CreateClustersJsonStringGrids(Clusters & clusters, float lat, float lng) {

	std::stringstream ss;
	std::string str;

	ss << std::fixed << std::setprecision(6);

	ss << "{\"channel\":\"" << CLUSTERS << "\",";
	ss << "\"clusters\": [";

	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		ss << "{\"ID\":" << it->first << ",";
		ss << "\"grids\": [";
		for (auto it_vect = it->second->get_begin_iterator(); it_vect != it->second->get_end_iterator(); ++it_vect)
		{
			ss << "{\"lat\":" << lat + it_vect->y / 1000 * 0.8 << ",\"lng\":" << lng + it_vect->x / 1000
				/*<< ",\"density\":" << it_vect->density_status*/
				<< "},";
		}
		ss.seekp(-1, ss.cur);
		ss << "]},";
	}
	if (clusters.size() > 0)
	{
		ss.seekp(-1, ss.cur);
	}
	ss << "]}";
	str = ss.str();
	return str;
}

std::string CreateCoordsJsonString(float lat, float lng)
{
	std::stringstream ss;
	std::string str;

	ss << std::fixed << std::setprecision(6);

	ss << "{\"channel\":\"" << COORDS << "\",\"lat\":" << lat << ",\"lng\":" << lng << "}";

	str = ss.str();

	return str;
}

std::string CreateGridsJsonString(Gridlist & grid_list, float lat, float lng) {
	std::stringstream ss;
	std::string str;
	Key key;
	ss << std::fixed << std::setprecision(6);

	ss << "{\"channel\":\"" << GRIDS << "\",";
	ss << "\"grids\": [";

	for (auto it = grid_list.begin(); it != grid_list.end(); ++it)
	{
		key = Key(it->first.x_, it->first.y_);
		ss << "{\"lat\":" << lat + it->first.y_ / 1000 * 0.8 << ",\"lng\":" << lng + it->first.x_ / 1000
			<< ",\"density\":" << it->second->get_density() << ",\"status\":" << (int)it->second->get_status()
			<< ",\"label\":" << it->second->get_label() << ",\"hash\":" << LabelHash(key) << "},";

	}
	ss.seekp(-1, ss.cur);
	ss << "]}";
	str = ss.str();
	return str;
}

std::string CreatePathJsonString(std::vector<Key> & path, uint32_t label)
{
	std::stringstream ss;
	std::string str;

	ss << std::fixed << std::setprecision(6);

	ss << "{\"ID\":\"" << label << "\",";
	ss << "\"path\": [";

	for (auto it = path.begin(); it != path.end(); ++it)
	{
		ss << "[" << it->y_ << "," << it->x_ << "],";
	}
	ss.seekp(-1, ss.cur);
	ss << "]}";

	str = ss.str();
	return str;
}

/* Function to estimate density threshold for specific grid to determine whether it is sporadic or not */
float DensityThresholdFunction(uint64_t time_updated, uint64_t time_now)
{
	float threshold = 0.0f;
	/* (27) - research paper */
	float numerator = float(Parameters::c_l * (1 - pow(Parameters::decay_factor, time_now - time_updated + 1)));
	float denumerator = Parameters::total_grids * (1 - Parameters::decay_factor);
	threshold = numerator / denumerator;
	return threshold;
}

void Deserialize(uint8_t * buffer, uint32_t buffer_size, uint64_t & time_now, Gridlist & grid_list)
{
	uint32_t index = 0;

	float x, y;
	Key key;
	uint8_t char_vect[CHAR_VECT_SIZE];

	/* time_now */
	memcpy(&time_now, &buffer[index], sizeof(uint64_t));
	index += sizeof(uint64_t);

	while (index < buffer_size)
	{
		memcpy(&x, &buffer[index], sizeof(float));
		index += sizeof(float);
		memcpy(&y, &buffer[index], sizeof(float));
		index += sizeof(float);
		index += sizeof(nullptr);
		memcpy(char_vect, &buffer[index], CHAR_VECT_SIZE);
		index += CHAR_VECT_SIZE;
		key = Key(x, y);
		grid_list[key] = new CharacteristicVector(char_vect);
	}
}

float EstimatedDensitiesSum(uint64_t time_now)
{
	float estimated_sum = float((1 - pow(Parameters::decay_factor, time_now + 1)) / (1 - Parameters::decay_factor)); // up from (8) - research paper
	return estimated_sum;
}

void GenerateRandomPair(float & x, float & y, int counter, std::vector<Key> & points)
{
	int x_diff, y_diff;
	uint32_t index = counter % (points.size() * 2);
	if (index >= points.size())
	{
		x = (float)(rand() % DIMENSIONS);
		y = (float)(rand() % DIMENSIONS);
	}
	else
	{
		x_diff = rand() % 3 - 1;
		switch (abs(x_diff)){
		case 0:
			y_diff = rand() % 3 - 1;
			break;
		case 1:
			y_diff = 0;
			break;
		}
		x = points[index].x_ + x_diff;
		y = points[index].y_ + y_diff;
	}

}

void GetNeighbors(float x, float y, Key neighbors[4]) {
	neighbors[0] = Key(x - STEP, y);
	neighbors[1] = Key(x + STEP, y);
	neighbors[2] = Key(x, y - STEP);
	neighbors[3] = Key(x, y + STEP);
}

bool IsSporadic(float density, uint64_t time_updated, uint64_t time_now)
{
	bool is_sporadic = false;
	float density_threshold;
	density_threshold = DensityThresholdFunction(time_updated, time_now);
	if (density < density_threshold) {
		is_sporadic = true;
	}
	return is_sporadic;
}

uint32_t LabelHash(Key & key)
{
	return std::hash<Key>()(key);
}

void MergeChangesToBuffer(uint8_t * buffer, uint32_t buffer_size, Gridlist & grid_list)
{
	CharacteristicVector * vect;
	CharacteristicVector vect_empty;
	float x, y;
	Key key;
	uint8_t char_vect[CHAR_VECT_SIZE];
	uint32_t index = 0;

	vect_empty = CharacteristicVector();

	/*Skips time_now*/
	index += sizeof(uint64_t);

	while (index < buffer_size)
	{
		memcpy(&x, &buffer[index], sizeof(float));
		index += sizeof(float);
		memcpy(&y, &buffer[index], sizeof(float));
		index += sizeof(float);
		index += sizeof(nullptr);
		key = Key(x, y);
		if (grid_list.count(key))
		{
			vect = grid_list[key];
			vect->Serialize(char_vect);
		}
		else
		{
			vect_empty.Serialize(char_vect);
		}

		memcpy(&buffer[index], char_vect, CHAR_VECT_SIZE);
		index += CHAR_VECT_SIZE;
	}
}

/* Merges two clusters and deletes the second one. */
void MergeClusters(Cluster & main_cluster, Cluster & secondary_cluster, Gridlist & grid_list, Clusters & clusters)
{
	GridTuple grid;
	Key grid_key;
	uint32_t label = main_cluster.get_label();
	uint32_t secondary_label = secondary_cluster.get_label();
	while (secondary_cluster.get_size() > 0)
	{
		grid = secondary_cluster.PopBack();
		grid_key = Key(grid.x, grid.y);
		if (grid_list.count(grid_key))
		{
			grid_list[grid_key]->set_label(label);
			main_cluster.AddElement(grid);
		}
		else
		{
			/* If it goes here, then there is a problem elsewhere. Non-existing grid should not be in cluster.*/
			std::cout << "WARNING: MergeClusters method tries to access char_vector that doesn't exist." << std::endl;
		}
	}
	delete clusters[secondary_label];
	clusters.erase(secondary_label);
}

void PrintClusters(Clusters & clusters)
{
	GridTuple grid;
	std::cout << std::endl;
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		std::cout << "Cluster #" << it->first << ":" << std::endl;
		for (uint32_t i = 0; i < it->second->get_size(); i++)
		{
			std::cout << "Grid:     " << grid.x << " : " << grid.y << std::endl;
		}
	}
	std::cout << std::endl;
}

void PrintPath(std::vector<Key> & path)
{
	int grid[10][10];
	int counter = 1;
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++)
		{
			grid[i][j] = 0;
		}
	}
	for (auto it = path.begin(); it != path.end(); ++it)
	{
		grid[(int)it->x_][(int)it->y_] = counter++;
	}
	std::cout << "--------------------------------------------------------" << std::endl;
	for (int i = 9; i >= 0; i--) {

		for (int j = 0; j < 12; j++)
		{
			std::cout << "|    ";
		}
		std::cout << std::endl << "| --";
		for (int j = 0; j < 10; j++)
		{
			if (grid[j][i] < 10)
			{
				std::cout << " ";
			}
			std::cout << grid[j][i] << " --";
		}
		std::cout << " |" << std::endl;
	}
	std::cout << "--------------------------------------------------------" << std::endl;
}

void PrintTable(Gridlist & grid_list, uint64_t time_now)
{
	Key key;
	std::cout << "------------------------------------------------------------------------" << std::endl;
	std::cout << std::setprecision(3);
	float sum = 0;
	for (int i = 0; i < DIMENSIONS; i++)
	{
		for (int j = 0; j < DIMENSIONS; j++) {
			key = Key(i, j);
			if (grid_list.count(key))
			{
				std::cout << "|" << std::setw(6) << grid_list[key]->get_density();
				sum += grid_list[key]->get_density();
			}
			else
			{
				std::cout << "|" << std::setw(6) << 0;
			}
		}
		std::cout << "|" << std::endl;
		std::cout << "-----------------------------------------------------------------------" << std::endl;
	}
	std::cout << std::setprecision(7) << std::endl;
	std::cout << "Time now: " << time_now << " (+1)." << std::endl;
	std::cout << "Sum: " << sum << std::endl;
	std::cout << "Estimated sum: " << EstimatedDensitiesSum(time_now) << std::endl;
}

void PrintGridsByClusters(Clusters & clusters, Gridlist & grid_list)
{
	float x, y;
	Key key;
	uint32_t label;
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		std::cout << "Cluster with label " << it->first << std::endl;
		for (auto it_vect = it->second->get_begin_iterator(); it_vect < it->second->get_end_iterator(); ++it_vect)
		{
			x = it_vect->x;
			y = it_vect->y;
			key = Key(x, y);
			label = grid_list[key]->get_label();
			std::cout << "    Grid (" << x << ":" << y << ") label is " << label << std::endl;
		}
		std::cout << std::endl;
	}
}

void ReassembleClusters(Gridlist & grid_list, Clusters & clusters)
{
	uint32_t label = NO_CLASS;
	for (auto it = grid_list.begin(); it != grid_list.end(); ++it)
	{
		label = it->second->get_label();
		if (label != NO_CLASS)
		{
			if (!clusters.count(label))
			{
				clusters[label] = new Cluster(label);

			}
			clusters[label]->AddElement(GridTuple(it->first.x_, it->first.y_, it->second->get_status()));
		}
	}
}

void RelabelNeighbors(Key key, Gridlist & grid_list, uint32_t old_label, uint32_t new_label)
{
	Key key_neighbor;
	Key neighbors[4];
	std::vector<Key> grids;
	grids.push_back(key);
	for (uint32_t i = 0; i < grids.size(); i++)
	{
		GetNeighbors(grids[i].x_, grids[i].y_, neighbors);
		for (int j = 0; j < 4; j++)
		{
			key_neighbor = neighbors[j];
			if (grid_list.count(key_neighbor))
			{
				if (grid_list[key_neighbor]->get_label() == old_label)
				{
					grid_list[key_neighbor]->set_label(new_label);
					grids.push_back(key_neighbor);
				}
			}
		}
	}
}

void RemoveEmptyClusters(Clusters & clusters)
{
	auto it = clusters.begin();
	while (it != clusters.end())
	{
		if (it->second->get_size() > 0)
		{
			it++;
		}
		else
		{
			delete it->second;
			it = clusters.erase(it);
		}
	}
}

void RemoveSporadic(Gridlist & grid_list, uint64_t time_now)
{
	uint64_t time_updated;
	float density = 0;
	auto it = grid_list.begin();

	while (it != grid_list.end())
	{
		if (it->second->get_status() == SPORADIC)
		{
			if (!it->second->is_changed())
			{
				delete it->second;
				it = grid_list.erase(it);
			}
			else
			{
				time_updated = it->second->get_time_updated();
				it->second->UpdateDensity(time_now);
				density = it->second->get_density();
				if (!IsSporadic(density, time_updated, time_now))
				{
					it->second->SetStatus(SPARSE);
				}
				it++;
			}
		}
		else if (it->second->get_status() < TRANSITIONAL_FROM_SPARSE)
		{
			time_updated = it->second->get_time_updated();
			it->second->UpdateDensity(time_now);
			density = it->second->get_density();
			if (IsSporadic(density, time_updated, time_now))
			{
				it->second->SetStatus(SPORADIC);
			}
			it++;
		}
		else
		{
			it++;
		}
	}
}

/* Renames cluster after it's first dense grid. If grid is not found - removes cluster. */
bool RenameCluster(Clusters & clusters, Gridlist & grid_list, uint32_t label)
{
	if (!clusters.count(label))
	{
		std::cout << "RenameCluster: cluster doesn't exist" << std::endl;
		return false;
	}
	Cluster * cluster = clusters[label];
	bool found = false;
	Key key;
	uint32_t new_label = NO_CLASS;
	uint32_t temp_label = NO_CLASS;

	for (auto it = cluster->get_begin_iterator(); it != cluster->get_end_iterator(); ++it)
	{
		key = Key(it->x, it->y);
		if (grid_list[key]->get_status() >= DENSE_FROM_SPARSE)
		{
			temp_label = LabelHash(key);
			/* Since all clusters are renamed after one of its DENSE grids, check is not needed.
			 Just to be safe!*/
			if (!clusters.count(temp_label) || temp_label == label)
			{
				new_label = temp_label;
				found = true;
				break;
			}
			else
			{
				RenameCluster(clusters, grid_list, temp_label);
				std::cout << "RenameCluster: cluster with that label already exists. DIFFERENT NAME" << std::endl;
			}
		}
	}
	for (auto it = cluster->get_begin_iterator(); it != cluster->get_end_iterator(); ++it)
	{
		key = Key(it->x, it->y);
		grid_list[key]->set_label(new_label);
	}
	if (found)
	{
		clusters[new_label] = clusters[label];
		clusters[new_label]->set_label(new_label);
	}
	else
	{
		delete clusters[label];
	}
	if (label != new_label)
	{
		clusters.erase(label);
	}
	return found;
}

void SplitCluster(Clusters & clusters, Gridlist & grid_list, std::vector<int> checklist, uint32_t label, int count)
{
	Key key;
	std::vector<Cluster *> new_clusters;
	std::vector<uint32_t> new_labels;
	for (int i = 0; i < count; i++)
	{
		new_labels.push_back(0);
		new_clusters.push_back(nullptr);
	}
	Cluster * cluster = clusters[label];

	int index = 0;
	uint32_t temp_label = 0;
	clusters[temp_label] = clusters[label];
	clusters.erase(label);
	for (uint32_t i = 0; i < checklist.size(); i++)
	{

		index = checklist[i] - 1;
		key = Key(cluster->GetElement(i).x, cluster->GetElement(i).y);
		if (new_labels[index] == 0)
		{

			new_labels[index] = LabelHash(key);
			if (clusters.count(new_labels[index]))
			{
				std::cout << "Name collision SplitCluster  " << clusters.count(new_labels[index]) << " " << new_labels[index] << std::endl;
				RenameCluster(clusters, grid_list, new_labels[index]);
			}
			new_clusters[index] = new Cluster(new_labels[index]);

		}
		grid_list[key]->set_label(new_labels[index]);
		new_clusters[index]->AddElement(GridTuple(key.x_, key.y_, grid_list[key]->get_status()));
	}
	delete clusters[temp_label];
	clusters.erase(temp_label);

	for (int i = 0; i < count; i++)
	{
		clusters[new_labels[i]] = new_clusters[i];
	}
	for (int i = 0; i < count; i++)
	{
		RenameCluster(clusters, grid_list, new_labels[i]);
	}
}

void UpdateDensities(Gridlist & grid_list, uint64_t time_now, float d_m, float d_l)
{
	float density;
	for (auto it = grid_list.begin(); it != grid_list.end(); ++it)
	{
		it->second->UpdateDensity(time_now);
		density = it->second->get_density();
		if (density >= d_m)
		{
			it->second->SetStatus(DENSE);
		}
		else if (density >= d_l)
		{
			it->second->SetStatus(TRANSITIONAL);
		}
		else if (it->second->get_status() != SPORADIC)
		{
			it->second->SetStatus(SPARSE);
		}
	}
}











