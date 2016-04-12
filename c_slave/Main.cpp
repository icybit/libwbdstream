#include "d-stream\CharacteristicVector.h"
#include "d-stream\Cluster.h"
#include "d-stream\Key.h"

#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

#include <iostream>
#include <sstream>

#include <algorithm>

#include "zmq.hpp"
#include <tuple>
#include <unordered_map>

//
#include <Windows.h>
//


using namespace std;

#define DIMENSIONS 40
#define TOTAL_GRIDS DIMENSIONS*DIMENSIONS
#define NO_OF_CYCLES 20000000

#define STEP 1

// Parameters controlling the treshold to recognize dense and sparse grids;
#define C_M 2.0f // min threshold for dense grid
#define C_L 0.8f // max threshold for sparse grid

#define COORDS "coords"
#define CLUSTERS "clusters"
#define GRIDS "grids"


//#define CENTER_X 25.290391f
//#define CENTER_Y 54.686668f

#define CENTER_X -122.506181
#define CENTER_Y 37.708818

//DIRECTIONS
#define EAST 1
#define SOUTH 2
#define WEST 3
#define	NORTH 4

typedef unordered_map<Key, CharacteristicVector * > Gridlist;
typedef unordered_map<unsigned int, Cluster *> Clusters;

void UpdateDensities(Gridlist & grid_list, unsigned __int64 time_now, float d_m, float d_l);
void AdjustClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now, float d_m, float d_l, int counter);
void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now);
float EstimatedDensitiesSum(unsigned __int64 time_now);
int CalculateGapTime();
float DensityThresholdFunction(unsigned __int64 time_updated, unsigned __int64 time_now);
void CalculateDensityParams(float & d_m, float & d_l);
bool IsSporadic(float density, unsigned __int64 time_updated, unsigned __int64 time_now);

void CallClusteringOnGrid(Key grid, Gridlist & grid_list, Clusters & clusters);
void CheckNeighborsConnectivity(Key grid, Gridlist & grid_list, Clusters & clusters);
void RemoveEmptyClusters(Clusters & clusters);
unsigned int LabelHash(Key & key);
void GetNeighbors(float x, float y, Key neighbors[4]);
bool CheckForDenseNeighbors(Key & key, Gridlist & grid_list);
void MergeClusters(Cluster & main_cluster, Cluster & secondary_cluster, Gridlist & grid_list, Clusters & clusters);
void SplitCluster(Clusters & clusters, Gridlist & grid_list, vector<int> checklist, unsigned int label, int count);
void RelabelNeighbors(Key key, Gridlist & grid_list, unsigned int old_label, unsigned int new_label);
bool RenameCluster(Clusters & clusters, Gridlist & grid_list, unsigned int label);
bool CheckIfConnected(Clusters & clusters, Gridlist & grid_list, unsigned int label);

void GenerateRandomPair(float & x, float & y, int counter, vector<Key> & points);
void PrintTable(Gridlist & grid_list, unsigned __int64 time_now);
void PrintClusters(Clusters & clusters);
void PrintGridsByClusters(Clusters & clusters, Gridlist & grid_list);

string CreateCoordsJsonString(float lat, float lng);
string CreateClustersJsonStringGrids(Clusters & clusters, float lat, float lng);
string CreateClustersJsonStringPaths(Clusters & clusters, float lat, float lng);
string CreateClustersJsonStringSelector(Clusters & clusters, float lat, float lng, bool grids);
string CreatePathJsonString(vector<Key> & path, unsigned int label);
string CreateGridsJsonString(Gridlist & grid_list, float lat, float lng);
string ClusterToPath(float lat, float lng, Cluster & cluster);
void PrintPath(vector<Key> & path);

void PushCoords(float x, float y, zmq::socket_t & publisher, float lat, float lng);
void PushClusters(Clusters & clusters, zmq::socket_t & publisher, bool grids);
void PushGrids(Gridlist & grid_list, zmq::socket_t & publisher);

void Deserialize(unsigned char * buffer, unsigned int buffer_size, Gridlist & grid_list);
void ReassembleClusters(Gridlist & grid_list, Clusters & clusters);

int main() {
	float x, y;
	float d_m, d_l; // threshold for dense grid, threshold for sparse grid;
	unsigned __int64 time_now = 0;
	unsigned __int64 gap = 0;
	Gridlist grid_list;
	Clusters clusters;
	Key key;
	float lat = CENTER_Y - 0.03f;
	float lng = CENTER_X - 0.1f;
	float x_diff;
	float y_diff;
	size_t str_length = 52;
	string str;

	bool push_as_polygons = false;

	int counter = 0;


	///////////////
	unsigned char buffer[300];
	unsigned int length = 0;
	CharacteristicVector  * vect_array[4];

	vect_array[0] = new CharacteristicVector(1);
	vect_array[0]->set_label(1111);
	vect_array[1] = new CharacteristicVector(2);
	vect_array[1]->set_label(2222);
	vect_array[2] = new CharacteristicVector(3);
	vect_array[2]->set_label(3333);
	vect_array[3] = new CharacteristicVector(4);
	vect_array[3]->set_label(4444);

	unsigned __int64 time_updated_;
	float density_;
	unsigned int label_; //cluster label
	unsigned __int8 status_; // sporadic/sparse/transitional/dense
	bool changed_;

	float xy;

	for (int i = 0; i < 4; i++)
	{
		CharacteristicVector * vect = vect_array[i];
		time_updated_ = vect->get_time_updated();
		density_ = vect->get_density();
		label_ = vect->get_label();
		status_ = vect->get_status();
		changed_ = vect->is_changed();

		memcpy(&buffer[length], &vect, sizeof(nullptr));
		length += sizeof(nullptr);

		memcpy(&buffer[length], &time_updated_, sizeof(unsigned __int64));
		length += sizeof(unsigned __int64);

		memcpy(&buffer[length], &density_, sizeof(float));
		length += sizeof(float);

		memcpy(&buffer[length], &label_, sizeof(unsigned int));
		length += sizeof(unsigned int);

		memcpy(&buffer[length], &status_, sizeof(unsigned __int8));
		length += sizeof(unsigned __int8);

		memcpy(&buffer[length], &changed_, sizeof(bool));
		length += sizeof(bool);

		xy = (float)i;
		memcpy(&buffer[length], &xy, sizeof(float));
		length += sizeof(float);
		memcpy(&buffer[length], &xy, sizeof(float));
		length += sizeof(float);

	}
	///////////////
	Deserialize(buffer, length, grid_list);
	ReassembleClusters(grid_list, clusters);
	///////////////
	grid_list.clear();
	clusters.clear();
	zmq::context_t context(1);
	zmq::socket_t publisher(context, ZMQ_PUB);
	publisher.bind("tcp://*:5556");
	// Handshake
	str = "Hello, my dear subscriber!";
	str_length = str.length();
	zmq::message_t message(str_length);
	memcpy(message.data(), str.c_str(), str_length);
	publisher.send(message);
	/*	auto start5 = std::chrono::high_resolution_clock::now();
		str = ClusterToPath(CENTER_Y, CENTER_X);
		auto end5 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff5 = end5 - start5;
		cout << diff5.count() << " s. TIME TO GET PATH " << endl << endl;
		str_length = str.length();
		zmq::message_t message(str_length);
		memcpy(message.data(), str.c_str(), str_length);
		publisher.send(message);
		*/

	// Start time

	//	

	CalculateDensityParams(d_m, d_l);
	gap = CalculateGapTime();
	// fake d_m - just for testing;
	//d_m = 1;
	cout << "Dm: " << d_m << ", Dl: " << d_l << ", gap: " << gap << endl;
	// fake gap time - for testing;
	//gap = 1050;
	srand(0);
	int timer = 0;
	vector<Key> points;
	for (int i = 0; i < 40; i++)
	{
		points.push_back(Key(rand() % DIMENSIONS, rand() % DIMENSIONS));
	}
	auto start0 = std::chrono::high_resolution_clock::now();


	for (int i = 0; i < NO_OF_CYCLES; i++)
	{
		GenerateRandomPair(x, y, i, points);
		key = Key(x, y);

		//PushCoords(x, y, publisher, CENTER_Y, CENTER_X);
		//Sleep(2);

		if (!grid_list.count(key))
		{
			grid_list[key] = new CharacteristicVector(time_now);
		}
		else
		{
			grid_list[key]->AddRecord(time_now);
		}
		if (time_now % gap == 0 && time_now != 0)
		{
			timer++;
			auto start1 = std::chrono::high_resolution_clock::now();
			AdjustClustering(grid_list, clusters, time_now, d_m, d_l, counter++);
			auto end1 = std::chrono::high_resolution_clock::now();
			//cout << "Clusters count: " <<  clusters.size() << endl;
			//cout << "Grids count: " << grid_list.size() << endl;
			std::chrono::duration<double> diff1 = end1 - start1;
			//cout << diff1.count() << " s. time to Adjust Clusters " << endl << endl;

			for (int j = 0; j < points.size() - 1; j++)
			{
				int diff_x = rand() % 3 - 1;
				int diff_y = rand() % 3 - 1;
				if (points[j].x_ + diff_x >= 0 && points[j].x_ + diff_x < DIMENSIONS)
				{
					points[j].x_ += diff_x;
				}
				if (points[j].y_ + diff_y >= 0 && points[j].y_ + diff_y < DIMENSIONS)
				{
					points[j].y_ += diff_y;
				}
			}
			//PrintGridsByClusters(clusters, grid_list);
			//PrintTable(grid_list, time_now);

			auto start2 = std::chrono::high_resolution_clock::now();
			//PushClusters(clusters, publisher, push_as_polygons);
			//PushGrids(grid_list, publisher);
			auto end2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff2 = end2 - start2;
			//cout << diff2.count() << " s. time to get paths " << endl << endl;
			//Sleep(5);
		}
		time_now++;
	}
	time_now--;
	zmq_close(socket);
	context.close();
	std::cout << "Press enter to exit." << endl;
	cin.ignore();
	return 0;
}

// Function to estimate density threshold for specific grid to determine whether it is sporadic or not
float DensityThresholdFunction(unsigned __int64 time_updated, unsigned __int64 time_now)
{
	float threshold = 0.0f;
	// (27) - research paper
	float numerator = float(C_L * (1 - pow(DECAY_FACTOR, time_now - time_updated + 1)));
	float denumerator = TOTAL_GRIDS * (1 - DECAY_FACTOR);
	threshold = numerator / denumerator;
	return threshold;
}

bool IsSporadic(float density, unsigned __int64 time_updated, unsigned __int64 time_now)
{
	bool is_sporadic = false;
	float density_threshold;
	density_threshold = DensityThresholdFunction(time_updated, time_now);
	if (density < density_threshold) {
		is_sporadic = true;
	}
	return is_sporadic;
}

int CalculateGapTime() {
	// Simplyfied:
	// Original method counts the time for sparse grid to become dense. Basically, how many hits in a row it needs to become dense.
	// In that case gap time is too small to recluster every time.
	// Now we count only decay factor - time for dense grid to decay to sparse.
	int gap = 0;
	double dense_to_sparse = log(C_L / C_M) / log(DECAY_FACTOR);// (11) - research paper
	gap = (int)floor(dense_to_sparse);
	return gap;
}

void UpdateDensities(Gridlist & grid_list, unsigned __int64 time_now, float d_m, float d_l)
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

void RemoveEmptyClusters(Clusters & clusters)
{
	for (auto it = clusters.begin(); it != clusters.end(); it)
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

unsigned int LabelHash(Key & key)
{
	return std::hash<Key>()(key);
}

void GetNeighbors(float x, float y, Key neighbors[4]) {
	neighbors[0] = Key(x - STEP, y);
	neighbors[1] = Key(x + STEP, y);
	neighbors[2] = Key(x, y - STEP);
	neighbors[3] = Key(x, y + STEP);
}

// Merges two clusters and deletes the second one.
void MergeClusters(Cluster & main_cluster, Cluster & secondary_cluster, Gridlist & grid_list, Clusters & clusters)
{
	GridTuple grid;
	Key grid_key;
	unsigned int label = main_cluster.get_label();
	unsigned int secondary_label = secondary_cluster.get_label();
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
			// If it goes here, then there is a problem elsewhere. Non-existing grid should not be in cluster.
			cout << "WARNING: MergeClusters method tries to access char_vector that doesn't exist." << endl;
		}
	}
	delete clusters[secondary_label];
	clusters.erase(secondary_label);
}

void RelabelNeighbors(Key key, Gridlist & grid_list, unsigned int old_label, unsigned int new_label)
{
	Key key_neighbor;
	Key neighbors[4];
	vector<Key> grids;
	grids.push_back(key);
	for (int i = 0; i < grids.size(); i++)
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

// Renames cluster after it's first dense grid. If grid is not found - removes cluster.
bool RenameCluster(Clusters & clusters, Gridlist & grid_list, unsigned int label)
{
	if (!clusters.count(label))
	{
		cout << "RenameCluster: cluster doesn't exist" << endl;
		return false;
	}
	Cluster * cluster = clusters[label];
	bool found = false;
	Key key;
	unsigned int new_label = NO_CLASS;
	unsigned int temp_label = NO_CLASS;

	for (auto it = cluster->get_begin_iterator(); it != cluster->get_end_iterator(); ++it)
	{
		key = Key(it->x, it->y);
		if (grid_list[key]->get_status() >= DENSE_FROM_SPARSE)
		{
			temp_label = LabelHash(key);
			// Since all clusters are renamed after one of its DENSE grids, check is not needed.
			// Just to be safe!
			if (!clusters.count(temp_label) || temp_label == label)
			{
				new_label = temp_label;
				found = true;
				break;
			}
			else
			{
				RenameCluster(clusters, grid_list, temp_label);
				cout << "RenameCluster: cluster with that label already exists. DIFFERENT NAME" << endl;
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

void SplitCluster(Clusters & clusters, Gridlist & grid_list, vector<int> checklist, unsigned int label, int count)
{
	Key key;
	vector<Cluster *> new_clusters;
	vector<unsigned int> new_labels;
	for (int i = 0; i < count; i++)
	{
		new_labels.push_back(0);
		new_clusters.push_back(nullptr);
	}
	Cluster * cluster = clusters[label];
	
	int index = 0;
	unsigned int temp_label = 0;
	clusters[temp_label] = clusters[label];
	clusters.erase(label);
	for (int i = 0; i < checklist.size(); i++)
	{

		index = checklist[i] - 1;
		key = Key(cluster->GetElement(i).x, cluster->GetElement(i).y);
		if (new_labels[index] == 0)
		{

			new_labels[index] = LabelHash(key);
			if (clusters.count(new_labels[index]))
			{
				cout << "Name collision SplitCluster  " << clusters.count(new_labels[index]) << " " << new_labels[index] << endl;
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

// Checks if cluster is connected. Splits it if not.
bool CheckIfConnected(Clusters & clusters, Gridlist & grid_list, unsigned int label)
{
	Key key;
	Key key_neighbor;
	Key neighbors[4];
	Cluster * cluster = clusters[label];
	bool is_connected = false;

	vector<int> checklist(cluster->get_size(), 0);
	unsigned int counter = 0;
	int index = 0;

	for (int i = 0; i < checklist.size(); i++)
	{
		if (checklist[i] == 0)
		{
			counter++;
			key = Key(cluster->GetElement(i).x, cluster->GetElement(i).y);
			grid_list[key]->set_label(counter);
			checklist[i] = counter;

			RelabelNeighbors(key, grid_list, label, counter);
			// can be split max into 4 clusters, no worries about cycle inside cycle;
			for (int j = i + 1; j < cluster->get_size(); j++)
			{
				key = Key(cluster->GetElement(j).x, cluster->GetElement(j).y);
				if (!grid_list.count(key))
				{
					cout << "STOP" << endl;
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
		for (int i = 0; i < cluster->get_size(); i++)
		{
			key = Key(cluster->GetElement(i).x, cluster->GetElement(i).y);
			grid_list[key]->set_label(label);
		}
		is_connected = true;
	}
	return is_connected;
}

void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now)
{
	unsigned __int64 time_updated;
	float density = 0;

	for (auto it = grid_list.begin(); it != grid_list.end(); it)
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

void CallClusteringOnGrid(Key grid, Gridlist & grid_list, Clusters & clusters)
{
	CharacteristicVector * vect;
	CharacteristicVector * vect_neighbor;
	Key key;
	Key key_neighbor;
	vector<Key> grids;
	Key neighbors[4];

	unsigned int label;
	unsigned int label_neighbor;

	grids.push_back(grid);
	for (int i = 0; i < grids.size(); i++)
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

// Checks if grid's neighbors are still connected to the cluster.
// If not - removes them from cluster and resets their labels.
void CheckNeighborsConnectivity(Key grid, Gridlist & grid_list, Clusters & clusters)
{
	Key key_neighbor;
	Key neighbors[4];
	CharacteristicVector * vect_neighbor;
	unsigned int label;

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

bool CheckForDenseNeighbors(Key & key, Gridlist & grid_list)
{
	Key neighbors[4];
	Key neighbor;
	unsigned int label;
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

void AdjustClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now, float d_m, float d_l, int counter)
{
	unsigned __int8 status;
	Key neighbors[4];
	GridTuple grid;
	CharacteristicVector * vect;
	CharacteristicVector * vect_neighbor;
	Key key;
	Key key_neighbor;
	unsigned int label;
	unsigned int label_neighbor;

	RemoveSporadic(grid_list, time_now);
	UpdateDensities(grid_list, time_now, d_m, d_l);
	//cout << counter << endl;
	for (auto it = grid_list.begin(); it != grid_list.end(); ++it)
	{	
		if (it->second->is_changed())
		{
			key = it->first;
			vect = it->second;
			status = vect->get_status();
			label = vect->get_label();
			// SPARSE
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
					// If cluster was named after this grid, we have to rename it.
					if (vect->get_status() == SPARSE_FROM_DENSE && label == LabelHash(key) && is_connected)
					{
						RenameCluster(clusters, grid_list, label);
					}
				}
			}
			// TRANSITIONAL
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
						// if cluster was named after this grid we have to rename cluster - grid might be removed
						if (label == LabelHash(key))
						{
							renamed = RenameCluster(clusters, grid_list, label);
							label = vect->get_label();
						}
						// If rename was not successful, that means there were no dense grids - cluster disassembled.
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
			// DENSE
			else
			{
				// If it doesn't belong to cluster, create its distinct cluster
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
				// If it belongs to cluster, update neighbors neighboring scores
				else
				{
					// This whole branch might be unnecessary. If dense grid has label - it means, it was added to cluster. 
					// If it was added to cluster, CallClusteringOnGrid was already called on it.
					// DENSE_FROM_TRANSITIONAL && DENSE
					if (vect->get_status() > DENSE_FROM_SPARSE)
					{
						CallClusteringOnGrid(key, grid_list, clusters);
					}
					else
					{
						// DENSE FROM SPARSE - logically, it should have NO_CLASS label;
						// If it has non-NO_CLASS label, it means, it was clustered by its neighbour.
						// So, no worries.
					}
				}
			}
			it->second->set_unchanged();
		}
		//cout << cnt << " " << counter << " After looking through everyone" << endl;
	}
	//cout << counter_adjusted << endl;
	RemoveEmptyClusters(clusters);
}

float EstimatedDensitiesSum(unsigned __int64 time_now)
{
	float estimated_sum = float((1 - pow(DECAY_FACTOR, time_now + 1)) / (1 - DECAY_FACTOR)); // up from (8) - research paper
	return estimated_sum;
}

// Calculate Dm and Dl;
void CalculateDensityParams(float & d_m, float & d_l)
{
	float denumerator = TOTAL_GRIDS * (1 - DECAY_FACTOR);
	d_m = C_M / denumerator;
	d_l = C_L / denumerator;
}

void PrintTable(Gridlist & grid_list, unsigned __int64 time_now)
{
	Key key;
	std::cout << "------------------------------------------------------------------------" << endl;
	std::cout << setprecision(3);
	float sum = 0;
	for (int i = 0; i < DIMENSIONS; i++)
	{
		for (int j = 0; j < DIMENSIONS; j++) {
			key = Key(i, j);
			if (grid_list.count(key))
			{
				std::cout << "|" << setw(6) << grid_list[key]->get_density();
				sum += grid_list[key]->get_density();
			}
			else
			{
				std::cout << "|" << setw(6) << 0;
			}
		}
		std::cout << "|" << endl;
		std::cout << "-----------------------------------------------------------------------" << endl;
	}
	std::cout << setprecision(7) << endl;
	std::cout << "Time now: " << time_now << " (+1)." << endl;
	std::cout << "Cycles: " << NO_OF_CYCLES << endl;
	std::cout << "Sum: " << sum << endl;
	std::cout << "Estimated sum: " << EstimatedDensitiesSum(time_now) << endl;
}

void PrintClusters(Clusters & clusters)
{
	GridTuple grid;
	std::cout << endl;
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		std::cout << "Cluster #" << it->first << ":" << endl;
		for (unsigned int i = 0; i < it->second->get_size(); i++)
		{
			std::cout << "Grid:     " << grid.x << " : " << grid.y << endl;
		}
	}
	std::cout << endl;
}

void PrintGridsByClusters(Clusters & clusters, Gridlist & grid_list)
{
	float x, y;
	Key key;
	unsigned int label;
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		int counter = 0;
		std::cout << "Cluster with label " << it->first << endl;
		for (auto it_vect = it->second->get_begin_iterator(); it_vect < it->second->get_end_iterator(); ++it_vect)
		{
			x = it_vect->x;
			y = it_vect->y;
			key = Key(x, y);
			label = grid_list[key]->get_label();
			std::cout << "    Grid (" << x << ":" << y << ") label is " << label << endl;
		}
		std::cout << endl;
	}
}

string CreateCoordsJsonString(float lat, float lng)
{
	std::stringstream ss;
	string str;

	ss << std::fixed << setprecision(6);

	ss << "{\"channel\":\"" << COORDS << "\",\"lat\":" << lat << ",\"lng\":" << lng << "}";

	str = ss.str();

	return str;
}

string CreateClustersJsonStringSelector(Clusters & clusters, float lat, float lng, bool grids)
{
	if (!grids)
	{
		return CreateClustersJsonStringGrids(clusters, lat, lng);
	}
	else
	{
		return CreateClustersJsonStringPaths(clusters, lat, lng);
	}
}

string CreateClustersJsonStringGrids(Clusters & clusters, float lat, float lng) {

	std::stringstream ss;
	string str;
	const char* separator = "";

	ss << std::fixed << setprecision(6);

	ss << "{\"channel\":\"" << CLUSTERS << "\",";
	ss << "\"clusters\": [";

	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		ss << "{\"ID\":" << it->first << ",";
		ss << "\"grids\": [";
		for (auto it_vect = it->second->get_begin_iterator(); it_vect != it->second->get_end_iterator(); ++it_vect)
		{
			ss << "{\"lat\":" << lat + it_vect->y / 1000 * 0.8 << ",\"lng\":" << lng + it_vect->x / 1000
				//<< ",\"density\":" << it_vect->density_status
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

string CreateClustersJsonStringPaths(Clusters & clusters, float lat, float lng)
{
	std::stringstream ss;
	string str;
	const char* separator = "";

	ss << std::fixed << setprecision(6);

	ss << "{\"channel\":\"" << CLUSTERS << "\",";
	ss << "\"clusters\": [";

	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		ss << ClusterToPath(lat, lng, *(it->second)) << ",";
	}
	if (clusters.size() > 0)
	{
		ss.seekp(-1, ss.cur);
	}
	ss << "]}";
	str = ss.str();
	return str;
}

string CreateGridsJsonString(Gridlist & grid_list, float lat, float lng) {
	std::stringstream ss;
	string str;
	const char* separator = "";

	ss << std::fixed << setprecision(6);

	ss << "{\"channel\":\"" << GRIDS << "\",";
	ss << "\"grids\": [";

	for (auto it = grid_list.begin(); it != grid_list.end(); ++it)
	{

			ss << "{\"lat\":" << lat + it->first.y_ / 1000 * 0.8 << ",\"lng\":" << lng + it->first.x_ / 1000
				<< ",\"density\":" << it->second->get_density() << ",\"status\":" << (int)it->second->get_status()
				<< ",\"label\":" << it->second->get_label() << ",\"hash\":" << LabelHash(Key(it->first)) << "},";
		
	}
	ss.seekp(-1, ss.cur);
	ss << "]}";
	str = ss.str();
	return str;
}

string CreatePathJsonString(vector<Key> & path, unsigned int label)
{
	std::stringstream ss;
	string str;

	ss << std::fixed << setprecision(6);

	ss << "{\"ID\":\"" << label << "\",";
	ss << "\"path\": [";

	for (auto it = path.begin(); it != path.end(); ++it)
	{
		ss << "[" << it->y_ << "," << it->x_<< "],";
	}
	ss.seekp(-1, ss.cur);
	ss << "]}";

	str = ss.str();
	return str;
}

string ClusterToPath(float lat, float lng, Cluster & cluster)
{
	float step = 1.0f;
	float half_step = step / 2;
	float x, y;
	vector<Key> points;
	vector<Key> points_vect;
	vector<int> indexes;
	unordered_map<Key, Key> points_map;
	vector<Key> path;
	int counter = 1;
	int last_direction = 0;
	int current_row = 0;
	int max_row;
	bool able_to_continue = true;

	Key temp_key;
	float temp_x;
	float temp_y;

	for (auto it = cluster.get_begin_iterator(); it != cluster.get_end_iterator(); ++it)
	{
		points.push_back(Key(it->x - half_step, it->y - half_step));
		points.push_back(Key(it->x - half_step, it->y + half_step));
		points.push_back(Key(it->x + half_step, it->y - half_step));
		points.push_back(Key(it->x + half_step, it->y + half_step));
	}


	std::sort(points.begin(), points.end());

	x = points[0].x_;
	y = points[0].y_;
	points_vect.push_back(Key(x, y));
	for (auto it = points.begin(); it != points.end(); ++it)
	{
		if ((abs(x - it->x_) < EPSILON) && (abs(y - it->y_) < EPSILON))
		{
			counter++;
		}
		else
		{
			if (counter % 2 == 1)
			{
				points_vect.push_back(Key(x, y));
			}
			counter = 1;
			x = it->x_;
			y = it->y_;
		}
	}
	if (counter % 2 == 1)
	{
		points_vect.push_back(Key(x, y));
	}

	for (auto it = points_vect.begin(); it != points_vect.end(); ++it)
	{
		Key key = Key(it->x_, it->y_);
		points_map[key] = key;
	}

	x = points_vect[0].x_;
	counter = 0;
	indexes.push_back(counter);
	for (auto it = points_vect.begin() + 1; it != points_vect.end(); ++it)
	{
		counter++;
		if (!(abs(it->x_ - x) < EPSILON))
		{
			int skips = (it->x_ - x) / (1 - EPSILON);
			for (int i = 0; i < skips; i++)
			{
				indexes.push_back(counter);
			}
			x = it->x_;
		}
	}
	indexes.push_back(points_vect.size());

	max_row = indexes.size() - 2;

	x = points_vect[indexes[1] - 1].x_;
	y = points_vect[indexes[1] - 1].y_;
	path.push_back(Key(x, y));
	points_map.erase(Key(x, y));
	while (able_to_continue)
	{
		if (last_direction % 2 == 0)
		{
			temp_x = x;
			able_to_continue = false;
			for (int i = current_row + 1; i <= max_row; i++)
			{
				temp_x += step;
				if (points_map.count(Key(temp_x, y)))
				{
					temp_key = Key(temp_x, y);
					points_map.erase(temp_key);
					path.push_back(temp_key);
					current_row = i;
					x = temp_x;
					able_to_continue = true;
					last_direction = EAST;
					break;
				}
			}
			if (!able_to_continue)
			{
				temp_x = x;
				for (int i = current_row - 1; i >= 0; i--)
				{
					temp_x -= step;
					if (points_map.count(Key(temp_x, y)))
					{
						Key temp_key = Key(temp_x, y);
						points_map.erase(temp_key);
						path.push_back(temp_key);
						current_row = i;
						x = temp_x;
						able_to_continue = true;
						last_direction = WEST;
						break;
					}
				}
			}
		}
		else
		{
			temp_y = y;
			able_to_continue = false;
			for (int i = indexes[current_row + 1] - 1; i >= indexes[current_row]; i--)
			{
				temp_key = points_vect[i];
				if (points_map.count(temp_key))
				{
					points_map.erase(temp_key);
					path.push_back(temp_key);
					last_direction = (y < temp_key.y_) ? NORTH : SOUTH;
					y = temp_key.y_;
					able_to_continue = true;
					break;
				}
			}
		}

	}
	for (auto it = path.begin(); it != path.end(); ++it)
	{
		it->x_ = lng + it->x_ / 1000;
		it->y_ = lat + it->y_ / 1000 * 0.6;
	}
	return CreatePathJsonString(path, cluster.get_label());
}

void PrintPath(vector<Key> & path)
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
	std::cout << "--------------------------------------------------------" << endl;
	for (int i = 9; i >= 0; i--) {

		for (int j = 0; j < 12; j++)
		{
			std::cout << "|    ";
		}
		std::cout << endl << "| --";
		for (int j = 0; j < 10; j++)
		{
			if (grid[j][i] < 10)
			{
				std::cout << " ";
			}
			std::cout << grid[j][i] << " --";
		}
		std::cout << " |" << endl;
	}
	std::cout << "--------------------------------------------------------" << endl;
}

void PushCoords(float x, float y, zmq::socket_t & publisher, float lat, float lng)
{
	size_t str_length;
	string str;
	float x_diff, y_diff;
	x_diff = x / 1000; // 0 to 1 -- /5 - 0 to 0.2
	y_diff = y / 1000 * 0.8; // 0 to 1 -- /16 - 0 to 0.06
	str = CreateCoordsJsonString(lat + y_diff, lng + x_diff);
	str_length = str.length();
	zmq::message_t message(str_length);
	memcpy(message.data(), str.c_str(), str_length);
	publisher.send(message);
}

void PushClusters(Clusters & clusters, zmq::socket_t & publisher, bool grids)
{
	string str;
	size_t str_length;

	str = CreateClustersJsonStringSelector(clusters, CENTER_Y, CENTER_X, grids);
	str_length = str.length();
	zmq::message_t message(str_length);
	memcpy(message.data(), str.c_str(), str_length);
	publisher.send(message);
}

void PushGrids(Gridlist & grid_list, zmq::socket_t & publisher)
{
	string str;
	size_t str_length;

	str = CreateGridsJsonString(grid_list, CENTER_Y, CENTER_X);
	str_length = str.length();
	zmq::message_t message(str_length);
	memcpy(message.data(), str.c_str(), str_length);
	publisher.send(message);
}

void GenerateRandomPair(float & x, float & y, int counter, vector<Key> & points)
{
	int x_diff, y_diff;
	int index = counter % (points.size() * 2);
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

void Deserialize(unsigned char * buffer, unsigned int buffer_size, Gridlist & grid_list)
{
	unsigned int index = 0;

	float x, y;
	Key key;

	size_t pointer_size = sizeof(nullptr);
	size_t float_size = sizeof(float);
	size_t char_vect_size = 18;
	unsigned char char_vect[18];
	cout << buffer[0];
	cout << buffer[500];
	for (int index = 0; index < buffer_size; index)
	{
		index += sizeof(nullptr);
		memcpy(&char_vect[0], &buffer[index], 18);
		index += 18;
		memcpy(&x, &buffer[index], sizeof(float));
		index += sizeof(float);
		memcpy(&y, &buffer[index], sizeof(float));
		index += sizeof(float);

		key = Key(x, y);
		grid_list[key] = new CharacteristicVector(char_vect);
	}
}

void ReassembleClusters(Gridlist & grid_list, Clusters & clusters)
{
	unsigned int label = NO_CLASS;
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