#ifndef __WB_DSTREAM_H__
#define __WB_DSTREAM_H__

#define __WB_DSTREAM_H_INSIDE__

#define DSTREAM_VERSION_MAJOR 0
#define DSTREAM_VERSION_MINOR 1
#define DSTREAM_VERSION_PATCH 0

#define DSTREAM_MAKE_VERSION(major, minor, patch) \
    ((major) * 10000 + (minor) * 100 + (patch))
#define DSTREAM_VERSION \
	DSTREAM_MAKE_VERSION(DSTREAM_VERSION_MAJOR, DSTREAM_VERSION_MINOR, DSTREAM_VERSION_PATCH)

#include "wbdstream/wbdstreamapi.h"
#include "wbdstream/wbcharvectapi.h"

#undef __WB_DSTREAM_H_INSIDE__

#endif /* __WB_DSTREAM_H_INSIDE__ */