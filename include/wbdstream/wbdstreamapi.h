#ifndef __WB_DSTREAM_API_H__
#define __WB_DSTREAM_API_H__

#if !defined (__WB_DSTREAM_H_INSIDE__) && !defined (WB_DSTREAM_COMPILATION)
	#error "Only <wbdstream.h> can be included directly."
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
	#ifdef __cplusplus
		#include <unordered_map>
		typedef std::unordered_map<Key, CharacteristicVector * > Gridlist;
		typedef std::unordered_map<unsigned int, Cluster *> Clusters;
	#endif


	void dstream_clusterize(unsigned char * buffer, uint32_t buffer_size);
	int dstream_calculate_gap_time();
	void dstream_calculate_xy_coords(double dx, double dy, double & x, double & y);
	void dstream_calculate_xy_distance(double x, double y, double & dx, double & dy);

#ifdef __cplusplus
}
#endif

#endif /* __WB_DSTREAM_API_H__ */