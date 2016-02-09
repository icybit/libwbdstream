#include "d-stream\CharacteristicVector.h"

#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <iomanip>

#include <tuple>
#include <unordered_map>
#include <vector>

using namespace std;

#define DIMENSIONS 10
#define TOTAL_GRIDS DIMENSIONS*DIMENSIONS
#define NO_OF_CYCLES 500

// Parameters controlling the treshold to recognize dense and sparse grids;
#define C_M 3.0f // min threshold for dense grid
#define C_L 0.8f // max threshold for sparse grid

typedef CharacteristicVector * Gridlist[DIMENSIONS][DIMENSIONS];
typedef unordered_map<int, vector<tuple<int, int>>> Clusters;


void UpdateDensities(Gridlist & grid_list, unsigned __int64 time_now);
void InitialClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now, float d_m, float d_l);
void AdjustClustering(Gridlist & grid_list, unsigned __int64 time_now);
void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now);
float EstimatedDensitiesSum(unsigned __int64 time_now);
int CalculateGapTime();
float DensityThreshold(unsigned __int64 time_updated, unsigned __int64 time_now);
void CalculateDensityThresholds(float & d_m, float & d_l);


void CreateDistinctClusters(Gridlist & grid_list, Clusters & clusters, float d_m);
void GetNeighbors(tuple<int, int> & pair, tuple<int, int> neighbors[4]);


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

	CalculateDensityThresholds(d_m, d_l);
	// fake d_m - just for testing;
	d_m = 5;
	//gap = CalculateGapTime(); 
	gap = 505;
	for (int i = 0; i < DIMENSIONS; i++) 
	{
		for (int j = 0; j < DIMENSIONS; j++) 
		{
			grid_list[i][j] = nullptr;
		}
	}
	for (int i = 0; i < NO_OF_CYCLES; i++)
	{
		GenerateRandomPair(x, y, i);
		if (grid_list[x][y] == nullptr)
		{
			grid_list[x][y] = new CharacteristicVector(time_now);
		}
		else
		{
			grid_list[x][y]->AddRecord(time_now);
		}
		if (time_now == gap)
		{
			InitialClustering(grid_list, clusters, time_now, d_m, d_l);
		}
		else if(time_now % gap == 0)
		{
			RemoveSporadic(grid_list, time_now);
			AdjustClustering(grid_list, time_now);
		}
		time_now++;
	}
	time_now--;
	
	UpdateDensities(grid_list, time_now);
	InitialClustering(grid_list, clusters, time_now, d_m, d_l);
	PrintClusters(clusters);
	PrintTable(grid_list, time_now);
	return 0;
}

float DensityThreshold(unsigned __int64 time_updated, unsigned __int64 time_now)
{
	float threshold = 0.0f;
	// (27) - research paper
	float numerator = float(C_L * (1 - pow(DECAY_FACTOR, time_now - time_updated + 1)));
	float denumerator = TOTAL_GRIDS * (1 - DECAY_FACTOR);
	threshold = numerator / denumerator;
	return threshold;
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
	tuple<int, int> grid;
	tuple<int, int> neighbors[4];
	tuple<int, int> grid_to_transfer;
	UpdateDensities(grid_list, time_now);
	CreateDistinctClusters(grid_list, clusters, d_m);
	//
	PrintClusters(clusters);
	//
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		for (int i = 0; i < it->second.size(); i++)
		{
			int smth = it->second.size();
			grid = it->second.at(i);
			GetNeighbors(grid, neighbors);
			for (int j = 0; j < 4; j++)
			{
				tie(x, y) = neighbors[j];
				if (x >= 0 && x < DIMENSIONS && y >= 0 && y < DIMENSIONS)
				{
					label = grid_list[x][y]->get_label();
					if (label != it->first)
					{
						if (label != NO_CLASS)
						{
							if (it->second.size() >= clusters[label].size())
							{
								for (int k = clusters[label].size() - 1; k >= 0; k--)
								{
									grid_to_transfer = clusters[label].at(k);
									clusters[label].pop_back();
									it->second.push_back(grid_to_transfer);
									tie(x_t, y_t) = grid_to_transfer;
									grid_list[x_t][y_t]->set_label(it->first);
								}
							}
							else {
								for (int k = it->second.size() - 1; k >= 0; k--)
								{
									grid_to_transfer = it->second.at(k);
									it->second.pop_back();
									clusters[label].push_back(grid_to_transfer);
									tie(x_t, y_t) = grid_to_transfer;
									grid_list[x_t][y_t]->set_label(label);
								}
							}
						}
						else if (grid_list[x][y]->get_density() >= d_l)
						{
							grid_list[x][y]->set_label(it->first);
							it->second.push_back(make_tuple(x, y));
						}
					}
				}
			}
		}
	}
}


void UpdateDensities(Gridlist & grid_list, unsigned __int64 time_now)
{
	for (int i = 0; i < DIMENSIONS; i++)
		for (int j = 0; j < DIMENSIONS; j++)
		{
			if (grid_list[i][j] != nullptr)
			{
				grid_list[i][j]->UpdateDensity(time_now);
			}
		}
}

void CreateDistinctClusters(Gridlist & grid_list, Clusters & clusters, float d_m)
{
	int n = 1;
	for (int i = 0; i < DIMENSIONS; i++)
		for (int j = 0; j < DIMENSIONS; j++)
		{
			if (grid_list[i][j] != nullptr && grid_list[i][j]->get_density() >= d_m)
			{
				grid_list[i][j]->set_label(n);
				clusters[n].push_back(make_tuple(i, j));
				int pf = get<0>(clusters[n].at(0));
				n++;
			}
		}
}

void GetNeighbors(tuple<int, int> & pair, tuple<int, int> neighbors[4]) {
	neighbors[0] = make_tuple(get<0>(pair) - 1, get<1>(pair));
	neighbors[1] = make_tuple(get<0>(pair) + 1, get<1>(pair));
	neighbors[2] = make_tuple(get<0>(pair), get<1>(pair) - 1);
	neighbors[3] = make_tuple(get<0>(pair), get<1>(pair) + 1);
}

void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now)
{

}

void AdjustClustering(Gridlist & grid_list, unsigned __int64 time_now)
{

}

float EstimatedDensitiesSum(unsigned __int64 time_now)
{
	float estimated_sum = float((1 - pow(DECAY_FACTOR, time_now + 1)) / (1 - DECAY_FACTOR)); // up from (8) - research paper
	return estimated_sum;
}

void CalculateDensityThresholds(float & d_m, float & d_l)
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
	cout << "------------------------------------------------------------------------" << endl;
	cout << setprecision(3);
	float sum = 0;
	for (int i = 0; i < DIMENSIONS; i++)
	{
		for (int j = 0; j < DIMENSIONS; j++)
			if (grid_list[i][j] != nullptr)
			{
				cout << "|" << setw(6) << grid_list[i][j]->get_density();
				sum += grid_list[i][j]->get_density();
			}
			else
			{
				cout << "|" << setw(6) << 0;
			}
		cout << "|" << endl;
		cout << "-----------------------------------------------------------------------" << endl;
	}
	cout << setprecision(7) << endl;
	cout << "Time now: " << time_now << " (+1)." << endl;
	cout << "Cycles: " << NO_OF_CYCLES << endl;
	cout << "Sum: " << sum << endl;
	cout << "Estimated sum: " << EstimatedDensitiesSum(time_now);
}

void PrintClusters(Clusters & clusters)
{
	cout << endl;
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		cout << "Cluster #" << it->first << ":" << endl;
		for (unsigned int i = 0; i < it->second.size(); i++)
		{
			cout << "Grid:     " << get<0>(it->second.at(i)) << " : " << get<1>(it->second.at(i)) << endl;
		}
	}
}