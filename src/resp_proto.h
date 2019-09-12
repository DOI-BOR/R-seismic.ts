#pragma once

extern void oops ( char *, char * );
extern void smsg ( char *, char * );
extern double bzpeak(double tau, double lambda, double *ts, int len, double dt, char type);
extern double lspeak(double tau, double lambda, double *ts, int len, double dt, char type);
extern double rpeak(double tau, double lambda, double *ts, int len, double dt, char type);
extern double cpeak(double tau, double lambda, double *ts, int len, double dt, char type);
extern double getmax(double *ts, int len, int *maxi);
extern int getpow(char inflg, char outflg);
#ifdef HAVE_FFTW3
 extern double fpeak(double tau, double lambda, double *ts, int len, double dt, char type);
 extern void fft_int(double *ts, double dt, int len, int nn);
#endif

DllImport void resp(char in_type, char out_type,
	double *ts, int len, double dt, double ptap,
	double **lambda, int *nlam,
	double tau_lo, double tau_hi, double **per, int *nper,
	RespMethod rm, double ***rspectra, char *rstype, char *units,
	double **si, double *tausi_lo, double *tausi_hi, char *siunits);
