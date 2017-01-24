#include "ts.h"
#include "ts_proto.h"

/* be sure to include these last so they can override some windows defs */
#include <R.h>
#include <Rinternals.h>

DllExport SEXP CALLfd_deriv(SEXP ts_d, SEXP dt_d, SEXP nd_i, SEXP order_i)
{
	int ii, pcnt = 0;
	double *dts = NULL, *dtp;
	SEXP dts_d;

	double *ts = REAL(ts_d);
	int len = length(ts_d);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	int nd = asInteger(nd_i) == NA_INTEGER ? 1 : INTEGER(nd_i)[0];
	FdOrder order = asInteger(order_i) == NA_INTEGER ? FD_ORDER_8 : INTEGER(order_i)[0];

	if ( len < 3 )
		error("time series must have at least 3 points");
	/* if ( nd < 1 || nd > 2 )
	error("only first and second derivatives supported"); */
	if ( dt <= 0 )
		error("dt must be positive");

	dts = fd_deriv(ts, len, dt, nd, order);

	/* allocate space for R structures for derivative, and copy fd_deriv output */
	dts_d = PROTECT(allocVector(REALSXP, len)); pcnt++;
	dtp = REAL(dts_d);
	for ( ii = 0 ; ii < len ; ii++ )
		dtp[ii] = dts[ii];

	free(dts);

	UNPROTECT(pcnt);

	return dts_d;
}

DllExport SEXP CALLwindow_ts(SEXP ts_d, SEXP dt_d, SEXP t0_d, SEXP tw_d, SEXP demean_i, SEXP ptap_d, SEXP type_s, SEXP do_norm_i)
{
	int ii, pcnt = 0, win_len;
	double *wts = NULL, *wtp;
	SEXP wts_d;

	double *ts = REAL(ts_d);
	int len = length(ts_d);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	double t0 = ! R_FINITE(asReal(dt_d)) ? 0 : REAL(t0_d)[0];
	double tw = ! R_FINITE(asReal(dt_d)) ? len * dt - t0 : REAL(tw_d)[0];
	BOOL demean = asLogical(demean_i) == NA_LOGICAL ? FALSE : LOGICAL(demean_i)[0];
	double ptap = ! R_FINITE(asReal(ptap_d)) ? 0 : REAL(ptap_d)[0];
	TaperType type = asChar(type_s) == NA_STRING ? TW_HANNING : /* if not set, use Hanning */
			strncasecmp(CHAR(STRING_ELT(type_s,0)), "ba", 2) == 0 ? TW_BARTLETT : /* Bartlett */
			strncasecmp(CHAR(STRING_ELT(type_s,0)), "p", 1) == 0 ? TW_PARZEN : /* Parzen */
			strncasecmp(CHAR(STRING_ELT(type_s,0)), "h", 1) == 0 ? TW_HANNING : /* Hanning */
			strncasecmp(CHAR(STRING_ELT(type_s,0)), "bl", 2) == 0 ? TW_BLACKMANN_HARRIS : /* Blackmann-Harris */
			strncasecmp(CHAR(STRING_ELT(type_s,0)), "e", 1) == 0 ? TW_EXACT_BLACKMANN : /* Exact Blackmann */
			TW_HANNING; /* silently ignore anything else, and use Hanning */
	BOOL do_norm = asLogical(do_norm_i) == NA_LOGICAL ? FALSE : LOGICAL(do_norm_i)[0];

	if ( len < 3 )
		error("time series must have at least 3 points");
	if ( dt <= 0 )
		error("dt must be positive");

	wts = window_ts(ts, len, FALSE, dt, t0, tw, demean, ptap, type, do_norm, &win_len);

	/* allocate space for R structures for windowed time series, and copy output */
	wts_d = PROTECT(allocVector(REALSXP, win_len)); pcnt++;
	wtp = REAL(wts_d);
	for ( ii = 0 ; ii < win_len ; ii++ )
		wtp[ii] = wts[ii];

	free(wts);

	UNPROTECT(pcnt);

	return wts_d;
}

DllExport SEXP CALLhilbertr_fir(SEXP ts_d, SEXP dt_d)
{
	int ii, pcnt = 0;
	double *hts = NULL, *htp;
	SEXP hts_d;

	double *ts = REAL(ts_d);
	int len = length(ts_d);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];

	if ( len < 3 )
		error("time series must have at least 3 points");
	if ( dt <= 0 )
		error("dt must be positive");

	hts = hilbertr_fir(ts, len);

	/* allocate space for R structures for derivative, and copy fd_deriv output */
	hts_d = PROTECT(allocVector(REALSXP, len)); pcnt++;
	htp = REAL(hts_d);
	for ( ii = 0 ; ii < len ; ii++ )
		htp[ii] = hts[ii];

	free(hts);

	UNPROTECT(pcnt);

	return hts_d;
}

DllExport SEXP CALLfilter_ts(SEXP ts_d, SEXP dt_d, SEXP order_i, SEXP pb_type_s, SEXP filt_type_s, SEXP flo_d, SEXP fhi_d,
															SEXP dir_s, SEXP cheb_attn_d, SEXP cheb_tr_bw_d)
{
	int ii, pcnt = 0;
	double *fts = NULL, *ftp;
	SEXP fts_d;

	double *ts = REAL(ts_d);
	int len = length(ts_d);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	int order = asInteger(order_i) == NA_INTEGER ? 4 : INTEGER(order_i)[0];
	PassBandType pb_type = asChar(pb_type_s) == NA_STRING ? PT_BAND_PASS : /* if not set, use band pass */
			strncasecmp(CHAR(STRING_ELT(pb_type_s,0)),"bp",2) == 0 ? PT_BAND_PASS : /* "bp" - band pass */
			strncasecmp(CHAR(STRING_ELT(pb_type_s,0)),"bandp",5) == 0 ? PT_BAND_PASS : /* "bp" - band pass */
			strncasecmp(CHAR(STRING_ELT(pb_type_s,0)),"band-p",6) == 0 ? PT_BAND_PASS : /* "bp" - band pass */
			strncasecmp(CHAR(STRING_ELT(pb_type_s,0)),"br",2) == 0 ? PT_BAND_REJECT : /* "br" - band reject (notch) */
			strncasecmp(CHAR(STRING_ELT(pb_type_s,0)),"bandr",5) == 0 ? PT_BAND_REJECT : /* "br" - band reject (notch) */
			strncasecmp(CHAR(STRING_ELT(pb_type_s,0)),"band-r",6) == 0 ? PT_BAND_REJECT : /* "br" - band reject (notch) */
			strncasecmp(CHAR(STRING_ELT(pb_type_s,0)),"n",1) == 0 ? PT_BAND_REJECT : /* "br" - band reject (notch) */
			strncasecmp(CHAR(STRING_ELT(pb_type_s,0)),"l",1) == 0 ? PT_LO_PASS :  /* "lp" - low pass */
			strncasecmp(CHAR(STRING_ELT(pb_type_s,0)),"h",1) == 0 ? PT_HI_PASS : /* "hp" - high pass */
			PT_BAND_PASS; /* silently ignore anything else, and use band pass */
	FilterType filt_type = asChar(filt_type_s) == NA_STRING ? FT_BUTTERWORTH : /* if not set, use Butterworth */
			strncasecmp(CHAR(STRING_ELT(filt_type_s,0)),"bu",2) == 0 ? FT_BUTTERWORTH : /* "bu" - Butterworth */
			strncasecmp(CHAR(STRING_ELT(filt_type_s,0)),"be",2) == 0 ? FT_BESSEL : /* "be" - Bessel */
			strncasecmp(CHAR(STRING_ELT(filt_type_s,0)),"c1",2) == 0 ? FT_CHEBYSHEV_I :  /* "c1" - Chebyshev Type I */
			strcasecmp(CHAR(STRING_ELT(filt_type_s,0)),"chebyshev1") == 0 ? FT_CHEBYSHEV_I :  /* "c1" - Chebyshev Type I */
			strcasecmp(CHAR(STRING_ELT(filt_type_s,0)),"chebyshev-type-i") == 0 ? FT_CHEBYSHEV_I :  /* "c1" - Chebyshev Type I */
			strncasecmp(CHAR(STRING_ELT(filt_type_s,0)),"c2",2) == 0 ? FT_CHEBYSHEV_II : /* "c2" - Chebyshev Type II */
			strcasecmp(CHAR(STRING_ELT(filt_type_s,0)),"chebyshev2") == 0 ? FT_CHEBYSHEV_II :  /* "c1" - Chebyshev Type II */
			strcasecmp(CHAR(STRING_ELT(filt_type_s,0)),"chebyshev-type-ii") == 0 ? FT_CHEBYSHEV_II :  /* "c1" - Chebyshev Type II */
			FT_BUTTERWORTH; /* silently ignore anything else, and use Butterworth */
	double flo = ! R_FINITE(asReal(flo_d)) ? 2 / (len * dt) : REAL(flo_d)[0];
	double fhi = ! R_FINITE(asReal(fhi_d)) ? 1 / (3 * dt) : REAL(fhi_d)[0];
	double cheb_attn = ! R_FINITE(asReal(cheb_attn_d)) ? 30 : REAL(cheb_attn_d)[0];
	DirectionType dir = asChar(dir_s) == NA_STRING ? DT_ZERO_PHASE : /* if not set, use zero-phase (forward and backwards) */
			strncasecmp(CHAR(STRING_ELT(dir_s,0)),"f",1) == 0 ? DT_FORWARD : /* "f*" - Forward */
			strncasecmp(CHAR(STRING_ELT(dir_s,0)),"r",1) == 0 ? DT_REVERSE : /* "r*" - reverse */
			strncasecmp(CHAR(STRING_ELT(dir_s,0)),"z",1) == 0 ? DT_ZERO_PHASE :  /* "z*" - zero-phase */
			DT_ZERO_PHASE; /* silently ignore anything else, and use zero-phase */
	double cheb_tr_bw = ! R_FINITE(asReal(cheb_tr_bw_d)) ? 0.3 : REAL(cheb_tr_bw_d)[0];

	if ( len < 3 )
		error("time series must have at least 3 points");
	if ( dt <= 0 )
		error("dt must be positive");
	if ( order < 1 || order > 8 )
		error("filter order must be between 1 and 8");

	if ( (fts = calloc(len, sizeof(*fts))) == NULL )
		oops("CALLfilter_ts", "Can't get space");
	design(order, pb_type, filt_type, cheb_attn, cheb_tr_bw, flo, fhi, dt);
	memcpy(fts, ts, len * sizeof(*ts));
	apply(fts, len, dir);

	/* allocate space for R structures for windowed time series, and copy output */
	fts_d = PROTECT(allocVector(REALSXP, len)); pcnt++;
	ftp = REAL(fts_d);
	for ( ii = 0 ; ii < len ; ii++ )
		ftp[ii] = fts[ii];

	free(fts);

	UNPROTECT(pcnt);

	return fts_d;
}
