#include "d-stream\CharacteristicVector.h"

#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <iomanip>

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
typedef unordered_map<int, vector<CharacteristicVector * >> Clusters;





void UpdateDensities(Gridlist & grid_list, unsigned __int64 time_now);
void InitialClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now);
void AdjustClustering(Gridlist & grid_list, unsigned __int64 time_now);
void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now);
float EstimatedDensitiesSum(unsigned __int64 time_now);
int CalculateGapTime();
float DensityThreshold(unsigned __int64 time_updated, unsigned __int64 time_now);
void CalculateDensityThresholds(float & d_m, float & d_l);


void IsDense(Gridlist & grid_list, float d_m);


void GenerateRandomPair(int & x, int & y);
void PrintTable(Gridlist & grid_list, unsigned __int64 time_now);

int main() {
	
	int x, y;

	float d_m, d_l; // threshold for dense grid, threshold for sparse grid;
	unsigned __int64 time_now = 0;
	unsigned __int64 gap = 0;

	Gridlist grid_list;
	Clusters clusters;

	CalculateDensityThresholds(d_m, d_l);
	//gap = CalculateGapTime(); 
	gap = 499;
	for (int i = 0; i < DIMENSIONS; i++) 
	{
		for (int j = 0; j < DIMENSIONS; j++) 
		{
			grid_list[i][j] = nullptr;
		}
	}
	for (int i = 0; i < NO_OF_CYCLES; i++)
	{
		GenerateRandomPair(x, y);
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
			InitialClustering(grid_list, clusters, time_now);
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
	else
	{
		gap = (int)floor(sparse_to_dense);
	}
	return gap;
}



void IsDense(Gridlist & grid_list, float d_m)
{
	int n = 0;
	for (int i = 0; i < DIMENSIONS; i++)
		for (int j = 0; j < DIMENSIONS; j++)
		{
			if (grid_list[i][j] != nullptr && grid_list[i][j]->get_density() >= d_m)
			{
				n++;
			}
		}
}


void InitialClustering(Gridlist & grid_list, Clusters & clusters, unsigned __int64 time_now)
{
	UpdateDensities(grid_list, time_now);
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

void RemoveSporadic(Gridlist & grid_list, unsigned __int64 time_now)
{

}

void AdjustClustering(Gridlist & grid_list, unsigned __int64 time_now)
{

}

void GenerateRandomPair(int & x, int & y)
{
	x = rand() % DIMENSIONS;
	y = rand() % DIMENSIONS;
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