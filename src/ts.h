#pragma once

#include "common.h"

/* define complex types */
typedef struct {
	double a;	/* amplitude */
	double p;	/* phase */
} CMPLX;

/* enum defining implemented orders for finite-difference derivatives */
typedef enum { FD_ORDER_2 = 0, FD_ORDER_4, FD_ORDER_6, FD_ORDER_8 } FdOrder;
/* enum defining implemented taper types */
typedef enum { TW_HANNING = 0, TW_BARTLETT, TW_PARZEN, TW_BLACKMANN_HARRIS, TW_EXACT_BLACKMANN } TaperType;
/* enum defining implemented filter types */
typedef enum { FT_BUTTERWORTH = 0, FT_BESSEL, FT_CHEBYSHEV_I, FT_CHEBYSHEV_II } FilterType;
/* enum defining implemented filter pass-band types */
typedef enum { PT_BAND_PASS = 0, PT_BAND_REJECT, PT_LO_PASS, PT_HI_PASS,  } PassBandType;
/* enum defining implemented filter direction types */
typedef enum { DT_FORWARD = 0, DT_REVERSE, DT_ZERO_PHASE } DirectionType;

#ifdef __TS_SRC
# ifdef DllImport
#  undef DllImport
#  define DllImport DllExport
# endif
#endif
#include "ts_proto.h"
