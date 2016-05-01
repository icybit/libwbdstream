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


	void Clusterize(unsigned char * buffer, uint32_t buffer_size);

#ifdef __cplusplus
}
#endif

#endif /* __WB_DSTREAM_API_H__ */