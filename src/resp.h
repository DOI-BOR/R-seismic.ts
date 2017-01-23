#pragma once

#include "common.h"

/* enum defining implemented methods to compute response spectrum */
typedef enum {
		RM_IIR = 0,
		RM_BZ,
		RM_LS,
		RM_CONV,
#ifdef HAVE_FFT
		RM_FFT
#endif
} RespMethod;
