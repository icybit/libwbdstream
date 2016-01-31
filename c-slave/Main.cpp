#include "d-stream\CharacteristicVector.h"

#include <stdlib.h>
#include <time.h>

#define DIMENSIONS 100
#define NO_OF_CYCLES 10

using namespace System;

void GenerateRandomPair(int & x, int & y);

int main() {
	int x, y;
	unsigned __int64 time_now = 0;
	CharacteristicVector * grid_list[DIMENSIONS][DIMENSIONS];
	grid_list[0][0] = new CharacteristicVector();
	delete grid_list[0][0];
	if (grid_list[0][0] == nullptr) {
		Console::Write("Tuscias");
		
	}
	else {
		Console::Write("Nepaejo");
	}
	for (int i = 0; i < NO_OF_CYCLES; i++)
	{
		GenerateRandomPair(x, y);
	}
}

void GenerateRandomPair(int & x, int & y)
{
	x = rand() % DIMENSIONS + 1;
	y = rand() % DIMENSIONS + 1;
}