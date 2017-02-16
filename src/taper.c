/*
	Written by	Chris Wood
	U.S. Bureau of Reclamation
	Denver, CO
*/

#include "ts.h"

#ifdef DllImport
#  undef DllImport
#  define DllImport DllExport
# endif
#include "ts_proto.h"

/* Mk_taper - construct various time-domain tapers
	Arguments:
		thlen - returned half-length of the taper
		len - input length of the taper
		type - taper type. Available taper types are:
			Bartlett Taper
			Parzen Taper
			Hanning Taper
			Blackmann-Harris
			Exact Blackmann
		ptap - Percentage taper; ranges between 0 and 50
		do_norm - Boolean flag: 1 - normalize taper to preserve power, 0 - no normalization
		norm_wt - normalization weight to apply to the un-tapered points, if norm = TRUE
	Returns a pointer to double with the taper - caller must free
*/
static double *mk_taper(int *thlen, int len, TaperType type, double ptap,
                        BOOL do_norm, double *norm_wt)
{
	double ww, *wt;
	double a0, a1, a2, a3, sum_wt2;
	int ii, thl;

	if ( thlen == NULL || (do_norm == TRUE && norm_wt == NULL) )
	  oops("mk_taper","Argument exception: thlen or norm_wt are NULL");

	/* taper 1/2-length in samples */
	thl = MIN(ROUND(len * ptap / 100), len / 2);
	*thlen = thl;

	if ( (wt = calloc(thl, sizeof(*wt))) == NULL )
		oops("taper","can't get space for taper");

	/* unnormalized weight outside of taper is just 1. */
	sum_wt2 = len - 2 * thl;

	/* construct the taper */
	switch ( type ) {
		case TW_BARTLETT:  /* bartlett taper */
			for ( ii = 0 ; ii < thl ; ii++ ) {
				ww = (double)ii / thl;
				wt[ii] = ww;
				sum_wt2 += 2. * ww * ww;
			}
			break;
		case TW_PARZEN: /* parzen taper */
			for ( ii = 0 ; ii < thl ; ii++ ) {
				ww = (double)ii / (thl + 1.);
				if ( ww - 0.5 <= 0 )
					ww = 6. * ww * ww * ( 1. - ww );
				else
					ww = 1. - 2. * pow( 1. - ww , 3.);
				wt[ii] = ww;
				sum_wt2 += 2. * ww * ww;
			}
			break;
		case TW_HANNING: /* hanning taper */
			for ( ii = 0 ; ii < thl ; ii++ ) {
				ww = 0.5 - 0.5 * cos(M_PI*ii/(thl+1.));
				wt[ii] = ww;
				sum_wt2 += 2. * ww * ww;
			}
			break;
		case TW_BLACKMANN_HARRIS: /* Blackmann-Harris 4 term -92 dB taper, Elliot
							and Rao, p. 218.  a0 = 0.35825 ? */
			a0 = 0.35875;	a1 = 0.48829;
			a2 = 0.14128;	a3 = 0.01168;
			for ( ii = 0 ; ii < thl ; ii++ ) {
				ww = a0 - a1 * cos(M_PI*ii/(thl+1.))
					+ a2 * cos(2.*M_PI*ii/(thl+1.))
					- a3 * cos(3.*M_PI*ii/(thl+1.));
				wt[ii] = ww;
				sum_wt2 += 2. * ww * ww;
			}
			break;
		case TW_EXACT_BLACKMANN: /* exact Blackman.  Elliot and Rao, p. 215-216 */
			a0 = 7938. / 18608.;
			a1 = 9240. / 18608.;
			a2 = 1430. / 18608.;
			for ( ii = 0 ; ii < thl ; ii++ ) {
				ww = a0 - a1 * cos(M_PI*ii/(thl+1.))
					+ a2 * cos(2.*M_PI*ii/(thl+1.));
				wt[ii] = ww;
				sum_wt2 += 2. * ww * ww;
			}
			break;
	}

	/* normalize so that sum(wt^2) = 1 */
	if ( do_norm ) {
		double norm = sqrt(len / sum_wt2);
		for ( ii = 0 ; ii < thl ; ii++ )
			wt[ii] *= norm;
		*norm_wt = norm; /* the un-tapered points will need to be multiplied by the same normalization factor */
	}

	return(wt);
}


/* Window_ts - window, demean, and taper a time series.
	Arguments:
		buf - pointer to a 1D array of double or CMPLX
		len - number of points in the array
		is_cmplx - TRUE if input buffer is COMPLX, FALSE if double
		dt - sample interval, in seconds
		demean - TRUE if the windowed data should be demeaned (after windowing, but before tapering)
		t0 - start time of the window, in seconds, relative to the start of the data. Default is 0
		tw - window length, in seconds. Note that t0 + tw < = len * dt. Default is len * dt - t0
		ptap - Percentage taper; ranges between 0 and 50. 0 for no tapering
		type - taper type for ptrap > 0. Available taper types are:
			Bartlett Taper
			Parzen Taper
			Hanning Taper
			Blackmann-Harris
			Exact Blackmann
		do_norm - Boolean flag: 1 - normalize taper to preserve power, 0 - no normalization
	Returns a pointer to the windowed data - caller must free
	*/
DllExport void *window_ts(void *buf, int len, BOOL is_cmplx, double dt,
									double t0, double tw, BOOL demean,
									double ptap, TaperType type, BOOL do_norm,
									int *wlen)
{
	int winlen, start;
	double *wt, *xr;
	CMPLX *xc;

	/* define window start and length  */
	if ( dt <= 0. )
		dt = 1.;	/* dt is positive */
	t0 = MAX(t0,0.);	/* t0 is non-negative */
	start = MIN(ROUND(t0 / dt),len-1);
	winlen = ROUND(tw / dt);
	if ( start + winlen > len ) {
		winlen = len - start;
		tw = winlen * dt;
	}

	/* get space for the window, and copy data from the input buffer */
	if ( is_cmplx ) {
		if ( (xc = calloc(winlen, sizeof(*xc))) == NULL )
			oops("window_ts","can't get space for time series");
		memcpy(xc, (CMPLX *)buf + start, winlen * sizeof(*xc));
	} else {
		if ( (xr = calloc(winlen, sizeof(*xr))) == NULL )
			oops("window_ts","can't get space for time series");
		memcpy(xr, (double *)buf + start, winlen * sizeof(*xr));
	}
	*wlen = winlen;

	if ( demean ) {
		/* de-mean points within the (untapered) window (time series is assumed to be real) */
		int ii;
		double mean;
		for ( mean = 0., ii = 0 ; ii < winlen ; ii++ )
			if ( is_cmplx )
				mean += xc[ii].a * cos( xc[ii].p );
			else
				mean += xr[ii];
		mean /= winlen;
		for ( ii = 0 ; ii < winlen ; ii++ ) {
			if ( is_cmplx ) {
				double tmp;
				tmp = xc[ii].a * cos( xc[ii].p ) - mean;
				xc[ii].a = ABS( tmp );
				xc[ii].p = (tmp < 0. ? M_PI : 0.);
			} else
				xr[ii] -= mean;
		}
	}

	if ( ptap > 0 ) {
		/* apply a taper to points within window */
		int ii, jj, kk;
		int whlen;
		double norm_wt;
		/* get taper weights */
		wt = mk_taper(&whlen, winlen, type, ptap, do_norm, &norm_wt);

		/* apply taper to window */
		for ( ii = 0 ; ii < whlen ; ii++ ) {
			kk = winlen - 1 - ii;
			if ( is_cmplx ) {
				xc[ii].a *= wt[ii];
				xc[kk].a *= wt[ii];
			} else {
				xr[ii] *= wt[ii];
				xr[kk] *= wt[ii];
			}
		}
		free(wt);
	}
# ifndef LOG2FFT
	else {
		/* make endpoint equal to average of discontinuity */
		if ( is_cmplx ) {
			double tmp;
			tmp = xc[0].a * cos( xc[0].p );
			tmp += xc[(len - 1)].a * cos( xc[(len - 1)].p );
			tmp /= 2.;
			xc[len - 1].a = ABS( tmp );
			xc[len - 1].p = (tmp < 0. ? M_PI : 0.);
		} else {
			double tmp;
			tmp = (xr[0] + xr[len - 1]) / 2.;
			xr[len - 1] = tmp;
		}
	}
# endif

	if ( is_cmplx )
		return(xc);
	else
		return(xr);
}
