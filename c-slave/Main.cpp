#include "d-stream\CharacteristicVector.h"
#include "d-stream\Cluster.h"
#include "d-stream\Key.h"

#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <iomanip>

#include <tuple>
#include <unordered_map>

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

typedef unordered_map<Key, CharacteristicVector * > Gridlist;
typedef unordered_map<unsigned int, Cluster *> Clusters;

void UpdateDensities(Gridlist & grid_list, unsigned __int64 time_now, float d_m, float d_l);
void InitialClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now, float d_m, float d_l);
void AdjustClustering(Gridlist & grid_list, unsigned __int64 time_now, float d_m, float d_l);
void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now);
float EstimatedDensitiesSum(unsigned __int64 time_now);
int CalculateGapTime();
float DensityThresholdFunction(unsigned __int64 time_updated, unsigned __int64 time_now);
void CalculateDensityParams(float & d_m, float & d_l);
bool IsSporadic(float density, unsigned __int64 time_updated, unsigned __int64 time_now);


void CreateDistinctClusters(Gridlist & grid_list, Clusters & clusters, float d_m);
void RemoveEmptyClusters(Clusters & clusters);
unsigned int LabelHash(Key & key);
void GetNeighbors(float x, float y, tuple<float, float> neighbors[4]);
void IncreaseNeighborsScore(unsigned int label, unsigned __int8 status);

void GenerateRandomPair(int & x, int & y, int counter);
void PrintTable(Gridlist & grid_list, unsigned __int64 time_now);
void PrintClusters(Clusters & clusters);

int main() {

	int x, y;
	float d_m, d_l; // threshold for dense grid, threshold for sparse grid;
	unsigned __int64 time_now = 0;
	unsigned __int64 gap = 0;
	Gridlist grid_list;
	Clusters clusters;
	CharacteristicVector * vect;
	Key key;

	/*
	auto start0 = std::chrono::high_resolution_clock::now();
	auto end0 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff0 = end0 - start0;
	std::cout << "Time to fill and iterate a vector of " << size << " ints : " << diff0.count() << " s\n";
	*/
	CalculateDensityParams(d_m, d_l);
	// fake d_m - just for testing;
	d_m = 5;
	gap = CalculateGapTime(); 
	// fake gap time - for testing;
	gap = 505;

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
			InitialClustering(grid_list, clusters, time_now, d_m, d_l);
		}
		else if(time_now % gap == 0 && time_now != 0)
		{
			AdjustClustering(grid_list, time_now, d_m, d_l);
		}
		time_now++;
	}
	time_now--;
	UpdateDensities(grid_list, time_now, d_m, d_l);
	InitialClustering(grid_list, clusters, time_now, d_m, d_l);
	PrintClusters(clusters);
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
	int x_t, y_t; // coordinates of grid to transfer;
	unsigned int label;
	GridTuple grid;
	tuple<float, float> neighbors[4];
	tuple<float, float> neighbors_join[4];
	CharacteristicVector * vect;
	Key key;

	Key key_neighbor;
	CharacteristicVector * vect_neighbor;
	int neighbor_score;
	
	UpdateDensities(grid_list, time_now, d_m, d_l);
	CreateDistinctClusters(grid_list, clusters, d_m);
	//
	PrintClusters(clusters);
	//
	// Iterate through all clusters
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		// Iterate through each element of cluster (grids might be added during iteration!)
		for (int i = 0; i < it->second->get_size(); i++)
		{
			grid = it->second->GetElement(i);
			vect = grid_list[Key(grid.x, grid.y)];
			if (vect->get_status() > TRANSITIONAL_FROM_DENSE || vect->get_neighbors() <= VALUE_CAN_ADD)
			{
				GetNeighbors(grid.x, grid.y, neighbors);
				for (int j = 0; j < 4; j++)
				{
					tie(x, y) = neighbors[j];
					key = Key(x, y);
					if (grid_list[key])
					{
						label = grid_list[key]->get_label();
						if (label != it->first)
						{
							if (label != NO_CLASS)
							{
								if (it->second->get_size() >= clusters[label]->get_size())
								{
									for (int k = 0; k > clusters[label]->get_size(); k++)
									{
										clusters[label]->GetPair(k, x_t, y_t);
										grid_list[Key(x_t, y_t)]->set_label(it->first);
									}
									it->second->MergeClusters(clusters[label]);
								}
								else {
									for (int k = 0; k < it->second->get_size(); k++)
									{
										it->second->GetPair(k, x_t, y_t);
										grid_list[Key(x_t, y_t)]->set_label(label);
									}
									clusters[label]->MergeClusters(it->second);
								}
							}
							// Only TRANSITIONAL grids can be labelled as NO_CLASS at this stage
							else if (grid_list[key]->get_density() >= d_l)
							{
								GetNeighbors(x, y, neighbors_join);
								neighbor_score = 0;
								for (int i = 0; i < 4; i++)
								{
									key_neighbor = Key(get<0>(neighbors_join[i]), get<1>(neighbors_join[i]));
									vect_neighbor = grid_list[key_neighbor];
									if (vect_neighbor && vect_neighbor->get_label() == it->first)
									{
										if (vect_neighbor->get_status() > TRANSITIONAL_FROM_DENSE)
										{
											neighbor_score += VALUE_DENSE;
										}
										else if(vect_neighbor->get_neighbors() <= VALUE_CAN_ADD)
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
								if (neighbor_score <= OUTSIDE_GRID_BOUNDARY)
								{
									grid_list[key]->set_label(it->first);
									grid_list[key]->AddNeighbors(neighbor_score);
									it->second->AddElement(GridTuple(x, y, grid_list[key]->get_density()));
									for (int i = 0; i < 4; i++)
									{
										key_neighbor = Key(get<0>(neighbors_join[i]), get<1>(neighbors_join[i]));
										vect_neighbor = grid_list[key_neighbor];
										if (vect_neighbor && vect_neighbor->get_label() == it->first)
										{
											vect_neighbor->AddNeighbors(VALUE_TRANSITIONAL);
										}
									}
								}							
							}
						}
					}
				}
			}
		}
	}
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
			clusters[label]->AddElement(GridTuple(it->first.x_, it->first.y_, it->second->get_density()));
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

void IncreaseNeighborsScore(unsigned int label, unsigned __int8 status)
{

}

void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now)
{
	unsigned __int64 time_updated;
	float density;
	CharacteristicVector * grid;
	for (int i = 0; i < DIMENSIONS; i++)
	{
		for (int j = 0; j < DIMENSIONS; j++)
		{
			if (grid_list[Key(i, j)] != nullptr)
			{
				grid = grid_list[Key(i, j)];
				if (grid->get_status() == SPORADIC)
				{
					if (!grid->is_changed())
					{
						grid_list.erase(Key(i, j));
					}
					else
					{
						time_updated = grid->get_time_updated();
						grid->UpdateDensity(time_now);
						density = grid->get_density();
						if (!IsSporadic(density, time_updated, time_now))
						{
							grid->SetStatus(SPARSE);
						}
					}
				}
				else if (grid->get_status() < TRANSITIONAL_FROM_SPARSE)
				{
					time_updated = grid->get_time_updated();
					grid->UpdateDensity(time_now);
					density = grid->get_density();
					if (IsSporadic(density, time_updated, time_now))
					{
						grid->SetStatus(SPORADIC);
					}
				}
			}
		}
	}
}

void AdjustClustering(Gridlist & grid_list, unsigned __int64 time_now, float d_m, float d_l)
{
	CharacteristicVector * grid;
	unsigned __int8 grid_status;

	//RemoveSporadic(grid_list, time_now);
	//UpdateDensities(grid_list, time_now, d_m, d_l);
	for (int i = 0; i < DIMENSIONS; i++)
	{
		for (int j = 0; j < DIMENSIONS; j++)
		{
			grid = grid_list[Key(i, j)];
			if (grid != nullptr && grid->is_changed())
			{
				grid_status = grid->get_status();
				if (grid_status == SPARSE)
				{
				}
				else if (grid_status == DENSE)
				{

				}
				else if (grid_status == TRANSITIONAL)
				{
				}
			}
		}
	}
}

float EstimatedDensitiesSum(unsigned __int64 time_now)
{
	float estimated_sum = float((1 - pow(DECAY_FACTOR, time_now + 1)) / (1 - DECAY_FACTOR)); // up from (8) - research paper
	return estimated_sum;
}

void CalculateDensityParams(float & d_m, float & d_l)
{
	float denumerator = TOTAL_GRIDS * (1 - DECAY_FACTOR);
	d_m = C_M / denumerator;
	d_l = C_L / denumerator;
}

void GenerateRandomPair(int & x, int & y, int counter)
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
	std::cout << "------------------------------------------------------------------------" << endl;
	std::cout << setprecision(3);
	float sum = 0;
	for (int i = 0; i < DIMENSIONS; i++)
	{
		for (int j = 0; j < DIMENSIONS; j++)
			if (grid_list[Key(i, j)] != nullptr)
			{
				std::cout << "|" << setw(6) << grid_list[Key(i, j)]->get_density();
				sum += grid_list[Key(i, j)]->get_density();
			}
			else
			{
				std::cout << "|" << setw(6) << 0;
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
}