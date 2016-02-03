#include "d-stream\CharacteristicVector.h"

#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <iomanip>

#define DIMENSIONS 10
#define TOTAL_GRIDS DIMENSIONS*DIMENSIONS
#define NO_OF_CYCLES 50

// Parameters controlling the treshold to recognize dense and sparse grids;
#define DENSE_PARAM 3
#define SPARSE_PARAM 0.8

using namespace std;


void UpdateDensities(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS], unsigned __int64 time_now);
void InitialClustering(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS], unsigned __int64 time_now);
void AdjustClustering(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS], unsigned __int64 time_now);
void RemoveSporadic(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS]);
float EstimatedDensitiesSum(unsigned __int64 time_now);
float DensityThreshold(unsigned __int64 time_updated, unsigned __int64 time_now);


void GenerateRandomPair(int & x, int & y);
void PrintTable(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS], unsigned __int64 time_now);

int main() {
	int x, y;
	unsigned __int64 time_now = 0;
	unsigned __int64 gap = 0;
	CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS];
	
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
		/*if (time_now == gap)
		{
			InitialClustering(grid_list);
		}
		else
		{
			RemoveSporadic(grid_list);
			AdjustClustering(grid_list);
		}*/
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
	float numerator = SPARSE_PARAM * (1 - pow(DECAY_FACTOR, time_now - time_updated + 1));
	float denumerator = TOTAL_GRIDS * (1 - DECAY_FACTOR);
	threshold = numerator / denumerator;
	return threshold;
}

void UpdateDensities(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS], unsigned __int64 time_now)
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

void InitialClustering(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS])
{

}

void RemoveSporadic(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS])
{

}

void AdjustClustering(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS])
{

}

void GenerateRandomPair(int & x, int & y)
{
	x = rand() % DIMENSIONS;
	y = rand() % DIMENSIONS;
}

float EstimatedDensitiesSum(unsigned __int64 time_now)
{
	float estimated_sum = (1 - pow(DECAY_FACTOR, time_now + 1)) / (1 - DECAY_FACTOR); // up from (8) - research paper
	return estimated_sum;
}

void PrintTable(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS], unsigned __int64 time_now)
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
	cout << "Time now: " << time_now << endl;
	cout << "Cycles: " << NO_OF_CYCLES << endl;
	cout << "Sum: " << sum << endl;
	cout << "Estimated sum: " << EstimatedDensitiesSum(time_now);
}