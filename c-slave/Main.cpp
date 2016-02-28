#include "d-stream\CharacteristicVector.h"
#include "d-stream\Cluster.h"

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

// Parameters controlling the treshold to recognize dense and sparse grids;
#define C_M 3.0f // min threshold for dense grid
#define C_L 0.8f // max threshold for sparse grid

#define OUTSIDE_GRID_BOUNDARY 3 // =<3 grid is OUTSIDE grid

typedef CharacteristicVector * Gridlist[DIMENSIONS][DIMENSIONS];
typedef unordered_map<int, Cluster> Clusters;


void UpdateDensities(Gridlist & grid_list, unsigned __int64 time_now, float d_m, float d_l);
void InitialClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now, float d_m, float d_l);
void AdjustClustering(Gridlist & grid_list, unsigned __int64 time_now, float d_m, float d_l);
void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now);
float EstimatedDensitiesSum(unsigned __int64 time_now);
int CalculateGapTime();
float DensityThreshold(unsigned __int64 time_updated, unsigned __int64 time_now);
void CalculateDensityParams(float & d_m, float & d_l);
bool IsSporadic(float density, unsigned __int64 time_updated, unsigned __int64 time_now);


void CreateDistinctClusters(Gridlist & grid_list, Clusters & clusters, float d_m);
void GetNeighbors(GridTuple & pair, tuple<int, int> neighbors[4]);


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
	CalculateDensityParams(d_m, d_l);
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
			AdjustClustering(grid_list, time_now, d_m, d_l);
		}
		time_now++;
	}
	time_now--;
	
	UpdateDensities(grid_list, time_now, d_m, d_l);
	InitialClustering(grid_list, clusters, time_now, d_m, d_l);
	PrintClusters(clusters);
	PrintTable(grid_list, time_now);
	cout << "Press enter to exit." << endl;
	cin.ignore();
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

bool IsSporadic(float density, unsigned __int64 time_updated, unsigned __int64 time_now)
{
	bool is_sporadic = false;
	float density_threshold;
	density_threshold = DensityThreshold(time_updated, time_now);
	if (density < density) {
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
	tuple<int, int> neighbors[4];
	
	UpdateDensities(grid_list, time_now, d_m, d_l);
	CreateDistinctClusters(grid_list, clusters, d_m);
	//
	PrintClusters(clusters);
	//
	// Iterate through all clusters
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		// Iterate through each element of cluster (grids might be added during iteration)
		for (int i = 0; i < it->second.get_size(); i++)
		{
			grid = it->second.GetElement(i);
			if (grid_list[grid.x][grid.y]->get_status() > TRANSITIONAL_FROM_DENSE || grid.neighbors < OUTSIDE_GRID_BOUNDARY)
			{
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
								if (it->second.get_size() >= clusters[label].get_size())
								{
									for (int k = 0; k > clusters[label].get_size(); k++)
									{
										clusters[label].GetPair(k, x_t, y_t);
										grid_list[x_t][y_t]->set_label(it->first);
									}
									it->second.MergeClusters(&clusters[label]);
								}
								else {
									for (int k = 0; k < it->second.get_size(); k++)
									{
										it->second.GetPair(k, x_t, y_t);
										grid_list[x_t][y_t]->set_label(label);
									}
									clusters[label].MergeClusters(&it->second);
								}
							}
							else if (grid_list[x][y]->get_density() >= d_l)
							{
								grid_list[x][y]->set_label(it->first);
								it->second.AddElement(GridTuple(x, y, grid_list[x][y]->get_density(), 1));
								it->second.increase_neighbors(i);
							}
						}
					}
				}
			}
		}
	}
}


void UpdateDensities(Gridlist & grid_list, unsigned __int64 time_now, float d_m, float d_l)
{
	float density;
	for (int i = 0; i < DIMENSIONS; i++)
		for (int j = 0; j < DIMENSIONS; j++)
		{
			if (grid_list[i][j] != nullptr)
			{
				grid_list[i][j]->UpdateDensity(time_now);
				density = grid_list[i][j]->get_density();
				if (density >= d_m) 
				{
					grid_list[i][j]->SetStatus(DENSE);
				}
				else if (density >= d_l)
				{
					grid_list[i][j]->SetStatus(TRANSITIONAL);
				}
				else if (grid_list[i][j]->get_status() != SPORADIC)
				{
					grid_list[i][j]->SetStatus(SPARSE);
				}
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
				clusters[n].AddElement(GridTuple(i, j, grid_list[i][j]->get_density(), 0));
				n++;
			}
		}
}

void GetNeighbors(GridTuple & pair, tuple<int, int> neighbors[4]) {
	neighbors[0] = make_tuple(pair.x - 1, pair.y);
	neighbors[1] = make_tuple(pair.x + 1, pair.y);
	neighbors[2] = make_tuple(pair.x, pair.y - 1);
	neighbors[3] = make_tuple(pair.x, pair.y + 1);
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
			if (grid_list[i][j] != nullptr)
			{
				grid = grid_list[i][j];
				if (grid->get_status() == SPORADIC)
				{
					if (!grid->is_changed())
					{
						delete grid_list[i][j];
						grid_list[i][j] = nullptr;
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
				else if (grid->get_status() == SPARSE)
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

	RemoveSporadic(grid_list, time_now);
	UpdateDensities(grid_list, time_now, d_m, d_l);
	for (int i = 0; i < DIMENSIONS; i++)
	{
		for (int j = 0; j < DIMENSIONS; j++)
		{
			grid = grid_list[i][j];
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
	cout << "Estimated sum: " << EstimatedDensitiesSum(time_now) << endl;
}

void PrintClusters(Clusters & clusters)
{
	int x, y;
	cout << endl;
	for (auto it = clusters.begin(); it != clusters.end(); ++it)
	{
		cout << "Cluster #" << it->first << ":" << endl;
		for (unsigned int i = 0; i < it->second.get_size(); i++)
		{	
			it->second.GetPair(i, x, y);
			cout << "Grid:     " << x << " : " << y << endl;
		}
	}
}