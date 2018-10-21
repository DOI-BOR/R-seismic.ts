#pragma once

#include "common.h"
#ifdef HAVE_FFTW3
# include "fftw3.h"
#endif
#include "ts.h"

/* enum defining implemented methods to compute response spectrum */
typedef enum {
		RM_IIR = 0,
		RM_BZ,
		RM_LS,
		RM_CONV,
#ifdef HAVE_FFTW3
		RM_FFT
#endif
} RespMethod;

#ifdef HAVE_FFTW3
typedef struct dcmplx {
	double r;	/* real part */
	double i;	/* imaginary part */
} DCMPLX;
#endif

#ifdef __RESP_SRC
# ifdef DllImport
#  undef DllImport
#  define DllImport DllExport
# endif
#endif
#include "resp_proto.h"
