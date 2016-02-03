#include "d-stream\CharacteristicVector.h"

#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <iomanip>

#define DIMENSIONS 10
#define NO_OF_CYCLES 50

using namespace std;

void GenerateRandomPair(int & x, int & y);
void PrintTable(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS]);
void UpdateDensities(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS], unsigned __int64 time_now);
void InitialClustering(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS], unsigned __int64 time_now);
void AdjustClustering(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS], unsigned __int64 time_now);
void RemoveSporadic(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS]);

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
		if (time_now == gap)
		{
			InitialClustering(grid_list);
		}
		else
		{
			RemoveSporadic(grid_list);
			AdjustClustering(grid_list);
		}

	}
	PrintTable(grid_list);
	return 0;
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

void PrintTable(CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS])
{
	cout << "-----------------------------------------" << endl;
	float sum = 0;
	for (int i = 0; i < DIMENSIONS; i++)
	{
		for (int j = 0; j < DIMENSIONS; j++)
			if (grid_list[i][j] != nullptr)
			{
				cout << "|" << setw(3) << grid_list[i][j]->get_density();
				sum += grid_list[i][j]->get_density();
			}
			else
			{
				cout << "|" << setw(3) << 0;
			}
		cout << "|" << endl;
		cout << "-----------------------------------------" << endl;
	}
	cout << "Cycles: " << NO_OF_CYCLES << endl;
	cout << "Sum: " << sum << endl;
}