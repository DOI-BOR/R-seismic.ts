#pragma once

extern void oops ( char *, char * );
extern void smsg ( char *, char * );
extern double bzpeak ( double, double, double *, int, double, char );
extern double lspeak ( double, double, double *, int, double, char );
extern double rpeak ( double, double, double *, int, double, char );
extern double cpeak ( double, double, double *, int, double, char );
extern double fpeak ( double, double, int, double, char );
extern int getpow ( char, char );
#ifdef HAVE_FFT
extern void omn ( CMPLX *, double, int, int );
#endif

DllImport void resp(char in_type, char out_type,
	double *ts, int len, double dt, double ptap,
	double **lambda, int *nlam,
	double tau_lo, double tau_hi, double **per, int *nper,
	RespMethod rm, double ***rspectra, char *rstype, char *units,
	double **si, double *tausi_lo, double *tausi_hi, char *siunits);
