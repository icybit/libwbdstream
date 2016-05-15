#ifndef __DSTREAM_COMMON_H__
#define __DSTREAM_COMMON_H__

#if __GNUC__ >= 4
#define DSTREAM_PUBLIC __attribute__ ((visibility ("default")))
#define DSTREAM_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define DSTREAM_PUBLIC
#define DSTREAM_LOCAL
#endif

#define C_M 3.0f /* min threshold for dense grid */
#define C_L 0.8f /* max threshold for sparse grid */

#define COORDS "coords"
#define CLUSTERS "clusters"
#define CHAR_VECT_SIZE 18
#define DECAY_FACTOR 0.998f
#define DIMENSIONS 40
#define EPSILON 0.0000001f
#define GRIDS "grids"
#define NO_OF_CYCLES 20000000
#define STEP 1
#define TOTAL_GRIDS 400
#define R_EARTH 6371000 

#define NO_CLASS 0
#define SPORADIC 0
#define SPARSE_FROM_TRANSITIONAL 1
#define SPARSE 2
#define SPARSE_FROM_DENSE 3
#define TRANSITIONAL_FROM_SPARSE 4
#define TRANSITIONAL 5
#define TRANSITIONAL_FROM_DENSE 6
#define DENSE_FROM_SPARSE 7
#define DENSE 8
#define DENSE_FROM_TRANSITIONAL 9
#define INITIAL 10

#endif /* __DSTREAM_COMMON_H__ */