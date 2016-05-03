#ifndef __WB_CHAR_VECT_API_H__
#define __WB_CHAR_VECT_API_H__

#if !defined (__WB_DSTREAM_H_INSIDE__) && !defined (WB_DSTREAM_COMPILATION)
	#error "Only <wbdstream.h> can be included directly."
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

	struct dstream_char_vect_s;
	typedef struct dstream_char_vect_s dstream_char_vect_t;

	dstream_char_vect_t * dstream_char_vect_new(uint64_t time_now);
	void dstream_char_vect_free(void ** vect);

	void dstream_char_vect_add_record(dstream_char_vect_t * vect, uint64_t time_now);
	void dstream_char_vect_print(dstream_char_vect_t * vect);
	void dstream_char_vect_serialize(dstream_char_vect_t * vect, uint8_t * buffer);
	void dstream_char_vect_update_density(dstream_char_vect_t * vect, uint64_t time_now);
	
#ifdef __cplusplus
}
#endif

#endif /* __WB_CHAR_VECT_API_H__ */