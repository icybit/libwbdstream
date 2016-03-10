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

#include "zmq.hpp"
#include <tuple>
#include <unordered_map>

//
#include <Windows.h>
//


using namespace std;

#define DIMENSIONS 10
#define TOTAL_GRIDS DIMENSIONS*DIMENSIONS
#define NO_OF_CYCLES 500

#define STEP 1

// Parameters controlling the treshold to recognize dense and sparse grids;
#define C_M 3.0f // min threshold for dense grid
#define C_L 0.8f // max threshold for sparse grid

#define OUTSIDE_GRID_BOUNDARY 4 // =<4 grid is OUTSIDE grid
#define VALUE_CAN_ADD OUTSIDE_GRID_BOUNDARY-VALUE_TRANSITIONAL
#define VALUE_TRANSITIONAL 3
#define VALUE_DENSE 1
#define VALUE_DIFF VALUE_TRANSITIONAL-VALUE_DENSE

#define COORDS "coords"
#define CLUSTERS "clusters"

#define CENTER_X 25.290391
#define CENTER_Y 54.686668

typedef unordered_map<Key, CharacteristicVector * > Gridlist;
typedef unordered_map<unsigned int, Cluster *> Clusters;

void UpdateDensities(Gridlist & grid_list, unsigned __int64 time_now, float d_m, float d_l);
void InitialClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now, float d_m, float d_l);
void AdjustClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now, float d_m, float d_l);
void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now);
float EstimatedDensitiesSum(unsigned __int64 time_now);
int CalculateGapTime();
float DensityThresholdFunction(unsigned __int64 time_updated, unsigned __int64 time_now);
void CalculateDensityParams(float & d_m, float & d_l);
bool IsSporadic(float density, unsigned __int64 time_updated, unsigned __int64 time_now);


void CreateDistinctClusters(Gridlist & grid_list, Clusters & clusters, float d_m);
void CallClusteringOnGrid(Key grid, Gridlist & grid_list, Clusters & clusters);
void RemoveEmptyClusters(Clusters & clusters);
unsigned int LabelHash(Key & key);
void GetNeighbors(float x, float y, tuple<float, float> neighbors[4]);
unsigned __int8 GetNeighborsScore(Gridlist & grid_list, tuple<float, float>(&neighbors_join)[4], unsigned int label, bool transitional);
void UpdateNeighborsScores(Gridlist & grid_list, tuple<float, float> neighbors[4], unsigned __int8 score, unsigned int label);
void MergeClusters(Cluster & main_cluster, Cluster & secondary_cluster, Gridlist & grid_list);
void ResetGridsStatusToUnchanged(Gridlist & grid_list);

void GenerateRandomPair(float & x, float & y, int counter);
void PrintTable(Gridlist & grid_list, unsigned __int64 time_now);
void PrintClusters(Clusters & clusters);
void PrintGridsByClusters(Clusters & clusters, Gridlist & grid_list);

string CreateCoordsJsonString(float lat, float lng);

int main() {

	float x, y;
	float d_m, d_l; // threshold for dense grid, threshold for sparse grid;
	unsigned __int64 time_now = 0;
	unsigned __int64 gap = 0;
	Gridlist grid_list;
	Clusters clusters;
	Key key;
	// ZeroMQ socket, publisher-subscriber pattern
	float lat = CENTER_Y - 0.03;
	float lng = CENTER_X - 0.1;
	float x_diff;
	float y_diff;
	size_t str_length = 52;
	string str;

	zmq::context_t context(1);
	zmq::socket_t publisher(context, ZMQ_PUB);
	publisher.bind("tcp://*:5556");
	while (1)
	{

		x_diff = ((float)rand() / (float)RAND_MAX)/5; // 0 to 1 -- /5 - 0 to 0.2
		y_diff = ((float)rand() / (float)RAND_MAX)/16; // 0 to 1 -- /16 - 0 to 0.06
		str = CreateCoordsJsonString(lat + y_diff, lng + x_diff);
		str_length = str.length();
		zmq::message_t message(str_length);
		memcpy(message.data(), str.c_str(), str_length);
		publisher.send(message);
		Sleep(50);
	}

	// Start time
	auto start0 = std::chrono::high_resolution_clock::now();
	//	

	CalculateDensityParams(d_m, d_l);
	// fake d_m - just for testing;
	d_m = 5;
	gap = CalculateGapTime();
	// fake gap time - for testing;
	gap = 505;
	srand(0);
	for (int i = 0; i < NO_OF_CYCLES; i++)
	{
		GenerateRandomPair(x, y, i);
		key = Key(x, y);
		if (!grid_list[key])
		{
			grid_list[key] = new CharacteristicVector(time_now);
		}
		else
		{
			grid_list[key]->AddRecord(time_now);
		}
		if (time_now == gap)
		{
			//InitialClustering(grid_list, clusters, time_now, d_m, d_l);
		}
		else if (time_now % gap == 0 && time_now != 0)
		{
			//AdjustClustering(grid_list, clusters, time_now, d_m, d_l);
		}
		time_now++;
	}
	time_now--;

	// Duration of inserts
	auto end0 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff0 = end0 - start0;
	//
	UpdateDensities(grid_list, time_now, d_m, d_l);
	InitialClustering(grid_list, clusters, time_now, d_m, d_l);
	//AdjustClustering(grid_list, clusters, time_now, d_m, d_l);
	// Duration of initial clustering and total time
	auto end1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff1 = end1 - start0;
	std::chrono::duration<double> diff2 = end1 - end0;

	std::cout << "Time for " << NO_OF_CYCLES << " inserts: " << diff0.count() << " s\n";
	std::cout << "Time to cluster after " << NO_OF_CYCLES << " inserts in " << DIMENSIONS << "x" << DIMENSIONS
		<< " area: " << diff2.count() << " s\n";
	std::cout << "Total time: " << diff1.count() << " s\n";
	//

	PrintClusters(clusters);
	PrintGridsByClusters(clusters, grid_list);
	PrintTable(grid_list, time_now);
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
	int gap = 0;
	double dense_to_sparse = log(C_L / C_M) / log(DECAY_FACTOR);// (11) - research paper
	double sparse_to_dense = log((TOTAL_GRIDS - C_M) / (TOTAL_GRIDS - C_L)) / log(DECAY_FACTOR);
	if (dense_to_sparse < sparse_to_dense)
	{
		gap = (int)floor(dense_to_sparse);
	}
	// This param will always be quite small (~10). Because N is big and difference between Cm and Cl is not.
	// Do we need it at all?
	else
	{
		gap = (int)floor(sparse_to_dense);
	}
	return gap;
}

void InitialClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now, float d_m, float d_l)
{
	int x, y; // coordinates of analyzed neighbor;
	unsigned int label;
	GridTuple grid;
	tuple<float, float> neighbors[4];
	CharacteristicVector * vect;
	Key key;

	tuple<float, float> neighbors_join[4];
	CharacteristicVector * vect_neighbor;
	Key key_neighbor;
	int neighbor_score;

	UpdateDensities(grid_list, time_now, d_m, d_l);
	CreateDistinctClusters(grid_list, clusters, d_m);
	//
	//PrintClusters(clusters);
	//
	// Iterate through all clusters
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		// Iterate through each element of cluster (grids might be added during iteration!)
		for (int i = 0; i < it->second->get_size(); i++)
		{
			grid = it->second->GetElement(i);
			vect = grid_list[Key(grid.x, grid.y)];
			if (vect->get_status() > TRANSITIONAL_FROM_DENSE || vect->get_neighbors() <= (OUTSIDE_GRID_BOUNDARY - VALUE_DENSE))
			{
				GetNeighbors(grid.x, grid.y, neighbors);
				for (int j = 0; j < 4; j++)
				{
					tie(x, y) = neighbors[j];
					key = Key(x, y);
					vect_neighbor = grid_list[key];
					if (vect_neighbor)
					{
						label = vect_neighbor->get_label();
						// If analyzed neighbor's and main grid's labels don't match
						if (label != it->first)
						{
							if (label != NO_CLASS)
							{
								if (it->second->get_size() >= clusters[label]->get_size())
								{
									MergeClusters(*it->second, *clusters[label], grid_list);
								}
								else {
									MergeClusters(*clusters[label], *it->second, grid_list);
								}
							}
							// Only TRANSITIONAL grids can be labelled as NO_CLASS at this stage
							else if (vect_neighbor->get_density() >= d_l && vect->get_neighbors() <= VALUE_CAN_ADD)
							{
								// Analyzed neighbor' neighbors. Checking if his neighbors can add it
								GetNeighbors(x, y, neighbors_join);
								neighbor_score = GetNeighborsScore(grid_list, neighbors_join, it->first, true);
								if (neighbor_score <= OUTSIDE_GRID_BOUNDARY)
								{
									vect_neighbor->set_label(it->first);
									vect_neighbor->AddNeighbors(neighbor_score);
									it->second->AddElement(GridTuple(x, y, vect_neighbor->get_status()));
									UpdateNeighborsScores(grid_list, neighbors_join, VALUE_TRANSITIONAL, it->first);
								}
							}
						}
					}
				}
			}
		}
	}
	PrintTable(grid_list, time_now);
	//ResetGridsStatusToUnchanged(grid_list);
	RemoveEmptyClusters(clusters);
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

void CreateDistinctClusters(Gridlist & grid_list, Clusters & clusters, float d_m)
{
	for (auto it = grid_list.begin(); it != grid_list.end(); ++it)
	{
		if (it->second->get_density() >= d_m)
		{
			Key key = Key(it->first.x_, it->first.y_);
			unsigned int label = LabelHash(key);
			it->second->set_label(label);
			clusters[label] = new Cluster(label);
			clusters[label]->AddElement(GridTuple(it->first.x_, it->first.y_, it->second->get_status()));
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

void GetNeighbors(float x, float y, tuple<float, float> neighbors[4]) {
	neighbors[0] = make_tuple(x - STEP, y);
	neighbors[1] = make_tuple(x + STEP, y);
	neighbors[2] = make_tuple(x, y - STEP);
	neighbors[3] = make_tuple(x, y + STEP);
}

unsigned __int8 GetNeighborsScore(Gridlist & grid_list, tuple<float, float>(&neighbors_join)[4], unsigned int label, bool transitional)
{
	unsigned __int8 neighbor_score = 0;
	float x, y;
	Key key_neighbor;
	CharacteristicVector * vect_neighbor;
	unsigned __int8 boundary_value;
	boundary_value = (transitional) ? VALUE_CAN_ADD : (OUTSIDE_GRID_BOUNDARY - VALUE_DENSE);

	for (int i = 0; i < 4; i++)
	{
		tie(x, y) = neighbors_join[i];
		key_neighbor = Key(x, y);
		vect_neighbor = grid_list[key_neighbor];
		if (vect_neighbor && vect_neighbor->get_label() == label)
		{
			if (vect_neighbor->get_status() > TRANSITIONAL_FROM_DENSE)
			{
				neighbor_score += VALUE_DENSE;
			}
			else if (vect_neighbor->get_neighbors() <= boundary_value)
			{
				neighbor_score += VALUE_TRANSITIONAL;
			}
			else
			{
				neighbor_score += OUTSIDE_GRID_BOUNDARY;
				break;
			}
		}
	}
	return neighbor_score;
}

void UpdateNeighborsScores(Gridlist & grid_list, tuple<float, float> neighbors[4], unsigned __int8 score, unsigned int label)
{
	float x, y;
	Key key_neighbor;
	CharacteristicVector * vect_neighbor;
	for (int i = 0; i < 4; i++)
	{
		tie(x, y) = neighbors[i];
		key_neighbor = Key(x, y);
		vect_neighbor = grid_list[key_neighbor];
		if (vect_neighbor && vect_neighbor->get_label() == label)
		{
			vect_neighbor->AddNeighbors(score);
		}
	}
}

void MergeClusters(Cluster & main_cluster, Cluster & secondary_cluster, Gridlist & grid_list)
{
	tuple<float, float> neighbors[4];
	GridTuple grid;
	CharacteristicVector * vect_neighbor;
	CharacteristicVector * vect_itself;
	vector<CharacteristicVector *> grids_to_move;
	float x, y;
	Key key_neighbor;
	Key key_itself;
	unsigned int label;

	unsigned int score_neighbor;
	unsigned int score_itself;

	label = main_cluster.get_label();

	for (int i = 0; i < secondary_cluster.get_size(); i++)
	{
		grid = secondary_cluster.GetElement(i);
		key_itself = Key(grid.x, grid.y);
		vect_itself = grid_list[key_itself];
		grids_to_move.push_back(vect_itself);
		GetNeighbors(grid.x, grid.y, neighbors);
		for (int j = 0; j < 4; j++)
		{
			tie(x, y) = neighbors[j];
			key_neighbor = Key(x, y);
			vect_neighbor = grid_list[key_neighbor];
			if (vect_neighbor && vect_neighbor->get_label() == label)
			{
				score_neighbor = (vect_neighbor->get_status() < DENSE_FROM_SPARSE) ? VALUE_TRANSITIONAL : VALUE_DENSE;
				score_itself = (vect_itself->get_status() < DENSE_FROM_SPARSE) ? VALUE_TRANSITIONAL : VALUE_DENSE;
				vect_neighbor->AddNeighbors(score_itself);
				vect_itself->AddNeighbors(score_neighbor);
			}
		}
	}
	for (int i = 0; i < grids_to_move.size(); i++)
	{
		grids_to_move[i]->set_label(label);
		main_cluster.AddElement(secondary_cluster.PopBack());
	}
}

void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now)
{
	unsigned __int64 time_updated;
	float density;

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
	}
}

void CallClusteringOnGrid(Key grid, Gridlist & grid_list, Clusters & clusters)
{
	CharacteristicVector * vect;
	CharacteristicVector * vect_neighbor;
	Key key;
	Key key_neighbor;
	vector<Key> grids;
	tuple<float, float> neighbors[4];
	tuple<float, float> neighbors_join[4];
	unsigned int label;
	unsigned int label_neighbor;
	float x, y;
	unsigned __int8 neighbor_score = 0;
	unsigned __int8 score = 0;

	grids.push_back(grid);
	for (int i = 0; i < grids.size(); i++)
	{
		key = Key(grids[i].x_, grids[i].y_);
		vect = grid_list[key];
		if (vect->get_status() > TRANSITIONAL_FROM_DENSE || vect->get_neighbors() <= (OUTSIDE_GRID_BOUNDARY - VALUE_DENSE))
		{
			GetNeighbors(key.x_, key.y_, neighbors);
			for (int j = 0; j < 4; j++)
			{
				tie(x, y) = neighbors[j];
				key_neighbor = Key(x, y);
				vect_neighbor = grid_list[key_neighbor];
				if (vect_neighbor && vect_neighbor->get_status() > SPARSE_FROM_DENSE)
				{
					if (vect_neighbor->get_label() != vect->get_label())
					{
						// Another grid that belongs to other cluster - just merge
						if (vect_neighbor->get_label() != NO_CLASS)
						{
							label = vect->get_label();
							label_neighbor = vect_neighbor->get_label();
							if (clusters[label_neighbor]->get_size() >= clusters[label]->get_size())
							{
								MergeClusters(*clusters[label_neighbor], *clusters[label], grid_list);
								delete clusters[label];
								clusters.erase(label);
							}
							else {
								MergeClusters(*clusters[label], *clusters[label_neighbor], grid_list);
								delete clusters[label_neighbor];
								clusters.erase(label_neighbor);
							}
						}
						// DENSE grid cannot be marked NO_CLASS, so Grid must be TRANSITIONAL
						else
						{
							label = vect->get_label();
							GetNeighbors(x, y, neighbors_join);
							neighbor_score = GetNeighborsScore(grid_list, neighbors_join, label, true);
							if (neighbor_score <= OUTSIDE_GRID_BOUNDARY)
							{
								vect_neighbor->set_label(label);
								vect_neighbor->AddNeighbors(neighbor_score);
								clusters[label]->AddElement(GridTuple(x, y, vect_neighbor->get_status()));
								UpdateNeighborsScores(grid_list, neighbors_join, VALUE_TRANSITIONAL, label);
								grids.push_back(key_neighbor);
							}
						}
					}
				}
			}
		}
	}
}

void AdjustClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now, float d_m, float d_l)
{
	unsigned __int8 status;
	tuple<float, float> neighbors[4];
	GridTuple grid;
	vector<CharacteristicVector *> vects_for_recluster; // Newly added, have to look for clusters around them
	CharacteristicVector * vect;
	Key key;
	Key key_neighbor;
	unsigned int label;
	unsigned __int8 score;

	RemoveSporadic(grid_list, time_now);
	UpdateDensities(grid_list, time_now, d_m, d_l);
	for (auto it = grid_list.begin(); it != grid_list.end(); ++it)
	{
		if (it->second->is_changed())
		{
			status = it->second->get_status();
			key = it->first;
			vect = it->second;
			// SPARSE
			if (status < TRANSITIONAL_FROM_SPARSE && status > SPORADIC)
			{
				// IMPLEMENT! Have to check whether have to split cluster or not
				if (vect->get_label() != NO_CLASS) {
					if (vect->get_status() == SPARSE_FROM_DENSE)
					{
						score = -VALUE_DENSE;
					}
					else
					{
						score = -VALUE_TRANSITIONAL;
					}
					GetNeighbors(key.x_, key.y_, neighbors);
					UpdateNeighborsScores(grid_list, neighbors, score, vect->get_label());
					clusters[vect->get_label()]->RemoveElement(key.x_, key.y_);
					vect->set_label(NO_CLASS);
				}
			}
			// TRANSITIONAL
			else if (status < DENSE_FROM_SPARSE)
			{

			}
			// DENSE
			else
			{
				GetNeighbors(key.x_, key.y_, neighbors);
				// If it doesn't belong to cluster, create its distinct cluster
				if (vect->get_label() == NO_CLASS)
				{
					label = LabelHash(key);
					clusters[label] = new Cluster(label);
					clusters[label]->AddElement(GridTuple(key.x_, key.y_, vect->get_status()));
					vect->set_label(label);
				}
				// If it belongs to cluster, update neighbors neighboring scores
				else
				{
					// DENSE_FROM_TRANSITIONAL
					if (vect->get_status() == DENSE_FROM_TRANSITIONAL)
					{
						score = -VALUE_DIFF;
					}
					// DENSE_FROM_SPARSE - does it ever hit this branch?
					// If it was sparse, it probably haven't been in cluster.
					else
					{
						score = VALUE_DENSE;
						cout << "AdjustClustering - DENSE_FROM_SPARSE branch: " << vect->get_label() << endl;
					}
					UpdateNeighborsScores(grid_list, neighbors, score, vect->get_label());
				}
				CallClusteringOnGrid(key, grid_list, clusters);
			}
		}
	}
}

void ResetGridsStatusToUnchanged(Gridlist & grid_list)
{
	int i = 0;
	for (auto it = grid_list.begin(); it != grid_list.end(); ++it)
	{
		if (it->second) {
			i++;
			cout << i << ". " << it->second->get_density() << endl;
			it->second->set_changed();
		}
	}
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

void GenerateRandomPair(float & x, float & y, int counter)
{
	if (counter % 25 == 0)
	{
		x = 4;
		y = 4;
	}
	else if (counter % 40 == 0)
	{
		x = 4;
		y = 5;
	}
	else {
		x = rand() % DIMENSIONS;
		y = rand() % DIMENSIONS;
	}
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
			if (grid_list[key] != nullptr)
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
	int x, y;
	std::cout << endl;
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		std::cout << "Cluster #" << it->first << ":" << endl;
		for (unsigned int i = 0; i < it->second->get_size(); i++)
		{
			it->second->GetPair(i, x, y);
			std::cout << "Grid:     " << x << " : " << y << endl;
		}
	}
	cout << endl;
}

void PrintGridsByClusters(Clusters & clusters, Gridlist & grid_list)
{
	float x, y;
	Key key;
	unsigned int label;
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		cout << "Cluster with label " << it->first << endl;
		for (auto it_vect = it->second->get_begin_iterator(); it_vect < it->second->get_end_iterator(); ++it_vect)
		{
			x = it_vect->x;
			y = it_vect->y;
			key = Key(x, y);
			label = grid_list[key]->get_label();
			cout << "    Grid (" << x << ":" << y << ") label is " << label << endl;
		}
		cout << endl;
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