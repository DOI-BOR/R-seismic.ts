/* RESP - calculate a response spectrum from an input time series.
	Written by	Chris Wood
			U.S. Bureau of Reclamation
			Denver, CO
*/

#include "resp.h"

#ifdef DllImport
#  undef DllImport
#  define DllImport DllExport
# endif
#include "resp_proto.h"

#define DEF_DAMPING	0.05

/* GETPER - get nicely spaced periods between taulo and tauhi */
#define DEF_TAU_LO 0
#define DEF_TAU_HI 10
/* static double dtau[] = 		{ .01, .025, .05, 0.1, .25, 0.5, 1.0, 2.5, -1. };
static double taumax[] = 	{ .05,  0.1, 0.3, 0.5, 1.0, 3.0, 5.0, 10., -1. }; */
static double dtau[] = 		{ .005, .005, .01, .02, .05, 0.1, 0.2, 0.5, -1. };
static double taumax[] = 	{ .05,  0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10., -1. };
#define SP_EPS	1.1e-16	/* infimum of single-precision delta */
static int getper(double **periods, double taulo, double tauhi)
{
	int ii, jj, kk = 0, nn;
	double *tau = NULL, dt, t, otaumax = 0;

	if ( periods == NULL )
		oops("getper","periods pointer is NULL");
	if ( taulo < 0 )
		taulo = DEF_TAU_LO;
	if ( tauhi < 0 )
		tauhi = DEF_TAU_HI;
	if ( taulo >= 0 && tauhi >= 0 && tauhi < taulo )
		oops("getper","high period < low period");

	t = otaumax;
	if ( t >= taulo ) {
		if ( (tau = realloc(tau, ++kk * sizeof(*tau))) == NULL )
			oops("getper","can't get space for periods");
		tau[kk-1] = t;
	}
	for ( ii = 0 ; taumax[ii] > 0. ; ii++ ) {
		dt = dtau[ii];
		nn = ROUND((taumax[ii] - otaumax) / dt);
		for ( t = otaumax, jj = 0 ; jj < nn ; jj++ ) {
			t += dt;
			if ( t > tauhi + SP_EPS )
				break;
			if ( t >= taulo - SP_EPS ) {
				if ( (tau = realloc(tau, ++kk * sizeof(*tau))) == NULL )
					oops("getper","can't get space for periods");
				tau[kk-1] = t;
			}
		}
		if ( t > tauhi + SP_EPS )
			break;
		otaumax = taumax[ii];
	}

	*periods = tau;
	return(kk);
}

/* get spectral intensity.  rectangular integration */
static double getsi(double *rspect, double *per, int nper, double tau1, double tau2) {
	int ii;
	double si, prev, next;

	si = 0.; prev = next = per[0];
	for ( ii = 0 ; ii < nper ; ii++ ) {
		double dper;
		if ( per[ii] < tau1 ) {
			prev = per[ii];
			continue;
		}
		if ( per[ii] > tau2 )
			break;
		if ( ii < nper - 1 )
			next = per[ii+1];
		if ( prev < tau1 || ii == 0 )
			dper = 0.5 * (per[ii] + next) - tau1;
		else if ( next > tau2 || ii == nper - 1 )
			dper = tau2 - 0.5 * (prev + per[ii]);
		else
			dper = 0.5 * (next - prev);
		prev = per[ii];
		si += rspect[ii] * dper;
	}

	return(si);
}

/* GETPSR - calculate a pseudo-relative response spectrum */
static double *getpsr(double *ts, int len, double dt, double *per, int nper,
		double lam, RespMethod rm, char otype)
{
	int ii;
	double *rspect;

	if ( nper <= 0 || per == NULL )
		oops("getpsr","empty or null response-period array");

	if ( (rspect = calloc(nper, sizeof(*rspect))) == NULL )
		oops("getpsr","can't get space for response spectrum");

	/* for each oscillator period, get peak response */
	for ( ii = 0 ; ii < nper ; ii++ ) {
		double tau = per[ii];
		if ( rm == RM_IIR )
			rspect[ii] = rpeak(tau,lam,ts,len,dt,otype);
		else if ( rm == RM_BZ )
			rspect[ii] = bzpeak(tau,lam,ts,len,dt,otype);
		else if ( rm == RM_LS )
			rspect[ii] = lspeak(tau,lam,ts,len,dt,otype);
		else if ( rm == RM_CONV )
			rspect[ii] = cpeak(tau,lam,ts,len,dt,otype);
#ifdef HAVE_FFTW3
		else if ( rm == RM_FFT )
			rspect[ii] = fpeak(tau,lam,ts,len,dt,otype);
#endif
		else
			oops("getpsr","unknown response method\n");
	}

	return (rspect);
}

DllExport void resp(char in_type, char out_type,
	double *ts, int len, double dt, double ptap,
	double **lambda, int *nlam,
	double tau_lo, double tau_hi, double **per, int *nper,
	RespMethod rm, double ***rspectra, char *rstype, char *units,
	double **si, double *tausi_lo, double *tausi_hi, char *siunits)
{
	int ii, nn, win_len;
	double tausi1, tausi2, t0, tw, *mts, **rspect;
	BOOL fftflg = FALSE;

	/* sanity checking */
	if ( ts == NULL || len < 3 )
		oops("do_resp","null time series, and/or length < 3");
	else if ( dt <= 0 )
		oops("do_resp","dt <= 0");
	else if ( per == NULL || nper == NULL )
		oops("do_resp","null period and/or count pointer");
	else if ( lambda == NULL || nlam == NULL )
		oops("do_resp","null damping and/or count pointer");
	else if ( strchr("avd",in_type) == NULL || strchr("avd",out_type) == NULL )
		oops("do_resp","unknown input or output type");
	else if ( ptap > 50 )
		oops("do_resp","ptap > 50");
	else if ( rspectra == NULL || si == NULL )
		oops("do_resp","null response spectrum and/or SI pointer");

	/* de-mean and window the input time series with a Hanning taper */
	t0 = 0; tw = len * dt;	/* use full time series length */
	mts = window_ts(ts,len,FALSE,dt,t0,tw,FALSE,ptap,TW_HANNING,FALSE,&win_len);	/* returns windowed time series; ts is unchanged */

	nn = getpow(in_type,'a');
#ifdef HAVE_FFTW3
	if ( fftflg && nn < 0 ) {
		fft_int(mts, dt, len, nn);
	} else
#endif
	{
		if ( nn > 0 ) {
			double *dts = fd_deriv(mts, win_len, dt, nn, FD_ORDER_8);
			free(mts);
			mts = dts;
		}
	}

	/* get damping, if not specified on input */
	if ( *nlam <= 0 || *lambda == NULL ) {
		if ( (*lambda = calloc(1, sizeof(**lambda))) == NULL )
			oops("do_resp","can't get space for lambda");
		(*lambda)[0] = DEF_DAMPING;	/* set default value */
		*nlam = 1;
	}

	/* get oscillator periods if not specified on input */
	if ( *nper <= 0 || *per == NULL ) {
		/* if *per != NULL, perhaps we should first free? no. */
		*nper = getper(per,tau_lo,tau_hi);
	}

	/* decide upon proper limits of integration for SI, and get output units */
	switch ( out_type ) {
		case 'a': /* acceleration */
			tausi1 = *tausi_lo < 0 ? 0.1 : *tausi_lo;
			tausi2 = *tausi_hi < 0 ? 0.5 : *tausi_hi;
			strcpy(rstype,"pseudo-absolute acceleration");
			strcpy(units,"cm/s^2");
			strcpy(siunits,"cm/s");
			break;
		case 'v': /* velocity */
			tausi1 = *tausi_lo < 0 ? 0.1 : *tausi_lo;
			tausi2 = *tausi_hi < 0 ? 2.5 : *tausi_hi;
			strcpy(rstype,"pseudo-relative velocity");
			strcpy(units,"cm/s");
			strcpy(siunits,"cm");
			break;
		case 'd': /* displacement */
			tausi1 = *tausi_lo < 0 ? 0.1 : *tausi_lo;
			tausi2 = *tausi_hi < 0 ? (*per)[*nper-1] : *tausi_hi;
			strcpy(rstype,"relative displacement");
			strcpy(units,"cm");
			strcpy(siunits,"cm-s");
			break;
	}
	*tausi_lo = tausi1;
	*tausi_hi = tausi2;

	/* get space for response spectra pointer array and SI array */
	if ( (rspect = calloc(*nlam, sizeof(*rspect))) == NULL )
		oops("do_resp","can't get space for rspect pointer array");
	if ( (*si = calloc(*nlam, sizeof(**si))) == NULL )
		oops("do_resp","can't get space for SI array");

	/* get the response spectrum and SI for each damping value */
	for ( ii = 0 ; ii < *nlam ; ii++ ) {
		/* get pseudo relative response spectrum */
		rspect[ii] = getpsr(mts, win_len, dt, *per, *nper, (*lambda)[ii], rm, out_type);
		(*si)[ii] = getsi(rspect[ii], *per, *nper, tausi1, tausi2);
	}
	free(mts);	/* done with modified time series */

	*rspectra = rspect;	/* note: rspect is an array of pointers, each element of which points to an array of double */
	return;
}
