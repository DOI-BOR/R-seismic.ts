#pragma once

extern void oops ( char *, char * );
extern void smsg ( char *, char * );

DllExport double *fd_deriv( double *ts, int len, double dt, int nd, FdOrder order );
DllExport void *window_ts( void *buf, int len, BOOL is_cmplx, double dt,
													 double t0, double tw, BOOL demean,
													 double ptap, TaperType type, BOOL do_norm,
													 int *wlen );
DllExport double *hilbertr_fir( double *ts, int len );
DllExport void apply( double *data, int nsamps, DirectionType dir );
DllExport void design( int iord, PassBandType type, FilterType aproto, double a, double trbndw, double fl, double fh, double ts );
