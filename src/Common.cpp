#include "Common.h"

float Common::c_l = C_L;
float Common::c_m = C_M;
float Common::decay_factor = DECAY_FACTOR;
int Common::total_grids = TOTAL_GRIDS;

void Common::InitParams(float c_l, float c_m, float decay_factor, int total_grids)
{
	Common::c_l = c_l;
	Common::c_m = c_m;
	Common::decay_factor = decay_factor;
	Common::total_grids = total_grids;
}

void Common::InitTotalGrids(int total_grids)
{
	Common::total_grids = total_grids;
}

extern "C"  DSTREAM_PUBLIC void dstream_init_params(float c_l, float c_m, float decay_factor, int total_grids)
{
	Common::InitParams(c_l, c_m, decay_factor, total_grids);
}

extern "C"  DSTREAM_PUBLIC void dstream_init_total_grids(int total_grids)
{
	Common::InitTotalGrids(total_grids);
}