#include "Parameters.h"

#include <stdint.h>

#include "Common.h"

float Parameters::c_l = C_L;
float Parameters::c_m = C_M;
float Parameters::decay_factor = DECAY_FACTOR;
uint64_t Parameters::total_grids = TOTAL_GRIDS;

void Parameters::InitParams(float c_m, float c_l, float decay_factor, uint64_t total_grids)
{
	Parameters::c_l = c_l;
	Parameters::c_m = c_m;
	Parameters::decay_factor = decay_factor;
	Parameters::total_grids = total_grids;
}