#ifndef __DSTREAM_PARAMETERS_H__
#define __DSTREAM_PARAMETERS_H__

#include <stdint.h>

class Parameters
{
public:
	static float c_m;
	static float c_l;
	static float decay_factor;
	static uint32_t distance;
	static uint64_t total_grids;

	static void InitParams(float c_m, float c_l, float decay_factor, uint64_t total_grids, uint32_t distance);
};

#endif /* __DSTREAM_PARAMETERS_H__ */

