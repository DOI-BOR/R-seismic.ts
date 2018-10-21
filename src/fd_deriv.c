/* Finite difference derivatives */

#include "common.h"
DllImport void oops ( char *, char * );
DllImport void smsg ( char *, char * );
DllImport char msgbuf[];

#define __TS_SRC
#include "ts.h"

static void fd2(double *ts, int len, double dt, int nd, double *dts) {
	double fd2_1_1 = 0.5 / dt;
	double fd2_2_0 = -2 / (dt * dt), fd2_2_1 = 1 / (dt * dt);
	int ii;
	dts[0] = 0;
	if ( nd == 1 ) {
		for ( ii = 1 ; ii < len - 2 ; ii++ )
			dts[ii] = fd2_1_1 * (ts[ii+1] - ts[ii-1]);
	} else if ( nd == 2 ) {
		for ( ii = 1 ; ii < len - 2 ; ii++ )
			dts[ii] = fd2_2_0 * ts[ii] + fd2_2_1 * (ts[ii+1] + ts[ii-1]);
	} else
		oops("fd2","only first and second derivative currently supported");
	dts[len - 1] = 0;
}

static void fd4(double *ts, int len, double dt, int nd, double *dts) {
	double fd2_1_1 = 0.5 / dt;
	double fd2_2_0 = -2 / (dt * dt), fd2_2_1 = 1 / (dt * dt);
	double fd4_1_1 = 2 / (3 * dt), fd4_1_2 = -1 / (12 * dt);
	double fd4_2_0 = -5 / (2 * dt * dt), fd4_2_1 = 4 / (3 * dt * dt), fd4_2_2 = -1 / (12 * dt * dt);
	int ii;
	dts[0] = 0;
	if ( nd == 1 ) {
		dts[1] = fd2_1_1 * (ts[2] - ts[0]);
		for ( ii = 2 ; ii < len - 2 ; ii++ )
			dts[ii] = fd4_1_1 * (ts[ii+1] - ts[ii-1]) + fd4_1_2 * (ts[ii+2] - ts[ii-2]);
		dts[len - 2] = fd2_1_1 * (ts[len-1] - ts[len-3]);
	} else if ( nd == 2 ) {
		dts[1] = fd2_2_0 * ts[1] + fd2_2_1 * (ts[2] + ts[0]);
		for ( ii = 2 ; ii < len - 2 ; ii++ )
			dts[ii] = fd4_2_0 * ts[ii] + fd4_2_1 * (ts[ii+1] + ts[ii-1]) + fd4_2_2 * (ts[ii+2] + ts[ii-2]);
		dts[len - 2] = fd2_2_0 * ts[len-2] + fd2_2_1 * (ts[len-1] + ts[len-3]);
	} else
		oops("fd4","only first and second derivative currently supported");
	dts[len - 1] = 0;
}

static void fd6(double *ts, int len, double dt, int nd, double *dts) {
	double fd2_1_1 = 0.5 / dt;
	double fd2_2_0 = -2 / (dt * dt), fd2_2_1 = 1 / (dt * dt);
	double fd4_1_1 = 2 / (3 * dt), fd4_1_2 = -1 / (12 * dt);
	double fd4_2_0 = -5 / (2 * dt * dt), fd4_2_1 = 4 / (3 * dt * dt), fd4_2_2 = -1 / (12 * dt * dt);
	double fd6_1_1 = 3 / (4 * dt), fd6_1_2 = -3 / (20 * dt), fd6_1_3 = 1 / (60 * dt);
	double fd6_2_0 = -49 / (18 * dt * dt), fd6_2_1 = 3 / (2 * dt * dt), fd6_2_2 = -3 / (20 * dt * dt), fd6_2_3 = 1 / (90 * dt * dt);
	int ii;
	dts[0] = 0;
	if ( nd == 1 ) {
		dts[1] = fd2_1_1 * (ts[2] - ts[0]);
		dts[2] = fd4_1_1 * (ts[3] - ts[1]) + fd4_1_2 * (ts[4] - ts[0]);
		for ( ii = 3 ; ii < len - 3 ; ii++ )
			dts[ii] = fd6_1_1 * (ts[ii+1] - ts[ii-1]) + fd6_1_2 * (ts[ii+2] - ts[ii-2]) + fd6_1_3 * (ts[ii+3] - ts[ii-3]);
		dts[len - 3] = fd4_1_1 * (ts[len-2] - ts[len-4]) + fd4_1_2 * (ts[len-1] - ts[len-5]);
		dts[len - 2] = fd2_1_1 * (ts[len-1] - ts[len-3]);
	} else if ( nd == 2 ) {
		dts[1] = fd2_2_0 * ts[1] + fd2_2_1 * (ts[2] + ts[0]);
		dts[2] = fd4_2_0 * ts[2] + fd4_2_1 * (ts[3] + ts[1]) + fd4_2_2 * (ts[4] + ts[0]);
		for ( ii = 3 ; ii < len - 3 ; ii++ )
			dts[ii] = fd6_2_0 * ts[ii] + fd6_2_1 * (ts[ii+1] + ts[ii-1]) + fd6_2_2 * (ts[ii+2] + ts[ii-2]) + fd6_2_3 * (ts[ii+3] + ts[ii-3]);
		dts[len - 3] = fd4_2_0 * ts[len-3] + fd4_2_1 * (ts[len-2] + ts[len-4]) + fd4_2_2 * (ts[len-1] + ts[len-5]);
		dts[len - 2] = fd2_2_0 * ts[len-2] + fd2_2_1 * (ts[len-1] + ts[len-3]);
	} else
		oops("fd6","only first and second derivative currently supported");
	dts[len - 1] = 0;
}

static void fd8(double *ts, int len, double dt, int nd, double *dts) {
	double fd2_1_1 = 0.5 / dt;
	double fd2_2_0 = -2 / (dt * dt), fd2_2_1 = 1 / (dt * dt);
	double fd4_1_1 = 2 / (3 * dt), fd4_1_2 = -1 / (12 * dt);
	double fd4_2_0 = -5 / (2 * dt * dt), fd4_2_1 = 4 / (3 * dt * dt), fd4_2_2 = -1 / (12 * dt * dt);
	double fd6_1_1 = 3 / (4 * dt), fd6_1_2 = -3 / (20 * dt), fd6_1_3 = 1 / (60 * dt);
	double fd6_2_0 = -49 / (18 * dt * dt), fd6_2_1 = 3 / (2 * dt * dt), fd6_2_2 = -3 / (20 * dt * dt), fd6_2_3 = 1 / (90 * dt * dt);
	double fd8_1_1 = 4 / (5 * dt), fd8_1_2 = -1 / (5 * dt), fd8_1_3 = 4 / (105 * dt), fd8_1_4 = -1 / (280 * dt);
	double fd8_2_0 = -205 / (72 * dt * dt), fd8_2_1 = 8 / (5 * dt * dt), fd8_2_2 = -1 / (5 * dt * dt),
		fd8_2_3 = 8 / (315 * dt * dt), fd8_2_4 = -1 / (560 * dt * dt);
	int ii;
	dts[0] = 0;
	if ( nd == 1 ) {
		dts[1] = fd2_1_1 * (ts[2] - ts[0]);
		dts[2] = fd4_1_1 * (ts[3] - ts[1]) + fd4_1_2 * (ts[4] - ts[0]);
		dts[3] = fd6_1_1 * (ts[4] - ts[2]) + fd6_1_2 * (ts[5] - ts[1]) + fd6_1_3 * (ts[6] - ts[0]);
		for ( ii = 4 ; ii < len - 4 ; ii++ )
			dts[ii] = fd8_1_1 * (ts[ii+1] - ts[ii-1]) + fd8_1_2 * (ts[ii+2] - ts[ii-2]) + fd8_1_3 * (ts[ii+3] - ts[ii-3]) + fd8_1_4 * (ts[ii+4] - ts[ii-4]);
		dts[len - 4] = fd6_1_1 * (ts[len-3] - ts[len-5]) + fd6_1_2 * (ts[len-2] - ts[len-6]) + fd6_1_3 * (ts[len-1] - ts[len-7]);
		dts[len - 3] = fd4_1_1 * (ts[len-2] - ts[len-4]) + fd4_1_2 * (ts[len-1] - ts[len-5]);
		dts[len - 2] = fd2_1_1 * (ts[len-1] - ts[len-3]);
	} else if ( nd == 2 ) {
		dts[1] = fd2_2_0 * ts[1] + fd2_2_1 * (ts[2] + ts[0]);
		dts[2] = fd4_2_0 * ts[2] + fd4_2_1 * (ts[3] + ts[1]) + fd4_2_2 * (ts[4] + ts[0]);
		dts[3] = fd6_2_0 * ts[3] + fd6_2_1 * (ts[4] + ts[2]) + fd6_2_2 * (ts[5] + ts[1]) + fd6_2_3 * (ts[6] + ts[0]);
		for ( ii = 4 ; ii < len - 4 ; ii++ )
			dts[ii] = fd8_2_0 * ts[ii] + fd8_2_1 * (ts[ii+1] + ts[ii-1]) + fd8_2_2 * (ts[ii+2] + ts[ii-2]) + fd8_2_3 * (ts[ii+3] + ts[ii-3]) + fd8_2_4 * (ts[ii+4] + ts[ii-4]);
		dts[len - 4] = fd6_2_0 * ts[len-4] + fd6_2_1 * (ts[len-3] + ts[len-5]) + fd6_2_2 * (ts[len-2] + ts[len-6]) + fd6_2_3 * (ts[len-1] + ts[len-7]);
		dts[len - 3] = fd4_2_0 * ts[len-3] + fd4_2_1 * (ts[len-2] + ts[len-4]) + fd4_2_2 * (ts[len-1] + ts[len-5]);
		dts[len - 2] = fd2_2_0 * ts[len-2] + fd2_2_1 * (ts[len-1] + ts[len-3]);
	} else
		oops("fd8","only first and second derivative currently supported");
	dts[len - 1] = 0;
}

DllExport double *fd_deriv(double *ts, int len, double dt, int nd, FdOrder order) {
	double *dts;

	if ( len < 3 || ts == NULL )
		oops("fd_deriv","NULL input data, or length < 3 points");
	if ( (dts = calloc(len, sizeof(*ts))) == NULL )
		oops("fd_deriv","can't get space for derivative");

	switch ( order ) {
		case FD_ORDER_2:
			fd2(ts, len, dt, nd, dts);
			break;
		case FD_ORDER_4:
			fd4(ts, len, dt, nd, dts);
			break;
		case FD_ORDER_6:
			fd6(ts, len, dt, nd, dts);
			break;
		case FD_ORDER_8:
			fd8(ts, len, dt, nd, dts);
			break;
		default:
			oops("fd_deriv","order not implemented");
			break;
	}

	return(dts);
}
