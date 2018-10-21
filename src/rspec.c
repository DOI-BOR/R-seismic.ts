#include "resp.h"
#include "resp_proto.h"

/* CHKMAX - check for a max */
static void chkmax(double xx, double *xmax, int ii, int *imax)
{
	double tst;

	tst = ABS(xx);
	if ( tst > *xmax ) {
		*xmax = tst;
		*imax = ii;
	}
}

/* NP - get corresponding power of omega according to flag */
static int np(char flg)
{
	int nn = 0;

	switch ( flg ) {
		case 'a':	/* acceleration */
			nn = 2;
			break;
		case 'v':	/* velocity */
			nn = 1;
			break;
		case 'd':	/* displacement */
			nn = 0;
			break;
		default:
			oops("np","bad flag");
	}

	return(nn);
}

/* GETPOW - get power of i * omega to go from inout to output type */
int getpow(char inflg, char outflg)
{
	return( np(outflg) - np(inflg) );
}

/* multiply relative displacement response by
	appropriate power of f0 to get desired response.
	If this is a zero-period call, assume we have
	the max of the acceleration */
static double scale(double tau, double max_disp, char type) {
	double max_val = max_disp;
	double c = 2. * M_PI / tau;

	switch ( type ) {
		case 'a': /* pseudo absolute acceleration */
			if ( tau > 0.00001 )
				max_val *= c * c;
			break;
		case 'v': /* pseudo relative velocity */
			if ( tau > 0.00001 )
				max_val *=  c;
			else
				max_val = 0.;
			break;
		case 'd': /* relative displacement */
			if ( tau < 0.00001 )
				max_val = 0;
			break;
	}
	return (max_val);
}

/* RPEAK - get peak response for a particular natural frequency using
	recursive method by IIR filter method */
double rpeak(double tau, double lambda, double *ts, int len, double dt, char type)
{
	int ii, imax;
	double omegas, alpha, beta, aa[3], bb[1];
	double yy, yy1, yy2;
	double max;

	/* if its the zero-period case, just find max of time series */
	if ( tau < 0.00001 ) { /* zero-period case */
		/* find max - skip first few points */
		max = 0.; imax = 0;
		for ( ii = 6 ; ii < len ; ii++ )
			chkmax(ts[ii],&max,ii,&imax);
		return(scale(tau, max, type));
	}

	if ( lambda > 1 )
		oops("rpeak","fractional damping greater than 1.0");

	/* get ARMA coefficients based on 2-pole IIR filter with poles
		at omega = i*lambda*omegas +/- sqrt(1 - lambda^2)*omegas */
	omegas = 2. * M_PI / tau;				/* omega sub s */
	alpha = lambda * omegas * dt;
	beta = lambda; beta *= beta; beta = 1. - beta;	/* 1 - lambda^2 */
	beta = sqrt(beta) * omegas * dt;
	aa[0] = 1.0;
	aa[1] = -2. * exp(-alpha) * cos(beta);
	aa[2] = exp(-2. * alpha);
	bb[0] = -(aa[0] + aa[1] + aa[2]) / omegas / omegas;

	/* recursively generate output series and look for max */
	max = 0.;  imax = 0;
	yy2 = bb[0] * ts[0];
	yy1 = bb[0] * ts[1] - aa[1] * yy2;
	for ( ii = 2 ; ii < len ; ii++ ) {
		yy = bb[0] * ts[ii] - aa[1] * yy1 - aa[2] * yy2;
		yy2 = yy1;
		yy1 = yy;
		chkmax(yy,&max,ii,&imax);
	}

	/* now have relative displacement - scale and return */
	return (scale(tau, max, type));
}

/* BZPEAK - get peak response for a particular natural frequency using
	bilinear Z transform of equation of motion */
double bzpeak(double tau, double lambda, double *ts, int len, double dt, char type)
{
	int ii, imax;
	double omegas, cw, sw, denom, aa[3], bb[3];
	double yy, yy1, yy2;
	double max;

	/* if its the zero-period case, just find max of time series */
	if ( tau < 0.00001 ) { /* zero-period case */
		/* find max - skip first few points */
		max = 0.; imax = 0;
		for ( ii = 6 ; ii < len ; ii++ )
			chkmax(ts[ii],&max,ii,&imax);
		return(scale(tau, max, type));
	}

	if ( lambda > 1 )
		oops("bzpeak","fractional damping greater than 1.0");

	/* get ARMA coefficients based on bilinear Z transform */
	omegas = 2. * M_PI / tau;				/* omega sub s */
	cw = cos(omegas * dt);
	sw = sin(omegas * dt);
	denom = 1. + lambda * sw;
	aa[0] = 1.0;
	aa[1] = -2. * cw / denom;
	aa[2] = (1. - lambda * sw) / denom;
	bb[0] = -dt * dt * (1. + cw) / (8. * denom);
	bb[1] = 2. * bb[0];
	bb[2] = bb[0];

	/* recursively generate output series and look for max */
	max = 0.;  imax = 0;
	yy2 = bb[0] * ts[0];
	yy1 = bb[0] * ts[1] + bb[1] * ts[0] - aa[1] * yy2;
	for ( ii = 2 ; ii < len ; ii++ ) {
		yy = bb[0] * ts[ii] + bb[1] * ts[ii-1] + bb[2] * ts[ii-2];
		yy -= aa[1] * yy1 + aa[2] * yy2;
		yy2 = yy1;
		yy1 = yy;
		chkmax(yy,&max,ii,&imax);
	}

	/* now have relative displacement - scale and return */
	return (scale(tau, max, type));
}

/* LSPEAK - get peak response for a particular natural frequency using
	method of Nigam and Jennings, 1969, wherein they do an exact
	integration of the equation of motion by approximating the recorded
	acceleration by sucessive linear segments (time series analysis
	apparently unknown).  This implementation assumes constant dt */
double lspeak(double tau, double lambda, double *ts, int len, double dt, char type)
{
	int ii;
	double omegas, omegad, sqrl, swd, cwd, elw;
	double t1, t2, t3, t4, t5, t6, t7;
	double a11, a12, a21, a22, b11, b12, b21, b22;
	double xx, xx1, xxdot, xxdot1, zz;
	double xxmax, xxdotmax, zzmax;
	int   ixmax, ixdotmax, izmax;

	/* if its the zero-period case, just find max of time series */
	if ( tau < 0.00001 ) {  /* zero-period case */
		if ( type == 'v' || type == 'd' )
			return( 0. );
		/* find max - skip first few points */
		zzmax = 0.; izmax = 0;
		for ( ii = 6 ; ii < len ; ii++ )
			chkmax(ts[ii],&zzmax,ii,&izmax);
		return(zzmax);
	}

	if ( lambda > 1 )
		oops("lspeak","fractional damping greater than 1.0");

	omegas = 2. * M_PI / tau;				/* omega sub s */
	elw = exp(-lambda * omegas * dt);
	sqrl = sqrt(1 - lambda * lambda);
	omegad =  omegas * sqrl;
	swd = sin(omegad * dt);
	cwd = cos(omegad * dt);
		/* get a bunch of terms */
	t1 = 1. / (omegas * omegas * dt);
	t2 = (2. * lambda * lambda - 1.) * t1;
	t3 = t2 + lambda / omegas;
	t4 = 2. * lambda * t1 / omegas;
	t5 = t4 + 1. / (omegas * omegas);
	t6 = cwd - lambda * swd / sqrl;
	t7 = omegad * swd + lambda * omegas * cwd;
		/* get the A and B matrix elements */
	a11 = elw * (cwd + lambda * swd / sqrl);
	a12 = elw * swd / omegas / sqrl;
	a21 = -elw * omegas * swd / sqrl;
	a22 = elw * (cwd - lambda * swd / sqrl);
	b11 = elw * (t3 * swd / omegad + t5 * cwd) - t4;
	b12 = -elw * (t2 * swd / omegad + t4 * cwd);
	b12 += - 1. / (omegas * omegas) + t4;
	b21 = elw * (t3 * t6 - t5 * t7) + t1;
	b22 = -elw * (t2 * t6 - t4 * t7) - t1;

	/* recursively generate output series and look for max */
	t1 = -2. * lambda * omegas; t2 = -omegas * omegas;
	xxmax = xxdotmax = zzmax = 0.;
	ixmax = ixdotmax = izmax = 0;
	xx = b12 * ts[0]; xxdot = b22 * ts[0];
	for ( ii = 1 ; ii < len ; ii++ ) {
		xx1 = a11 * xx + a12 * xxdot + b11 * ts[ii-1] + b12 * ts[ii];
		xxdot1 = a21 * xx + a22 * xxdot + b21 * ts[ii-1] + b22 * ts[ii];
		xx = xx1;
		xxdot = xxdot1;
		zz = t1 * xxdot + t2 * xx;
		chkmax(xx,&xxmax,ii,&ixmax);
		chkmax(xxdot,&xxdotmax,ii,&ixdotmax);
		chkmax(zz,&zzmax,ii,&izmax);
	}

	switch ( type ) {
		case 'a': /* pseudo absolute acceleration */
			return(omegas * omegas * xxmax);
		case 'v': /* pseudo relative velocity */
			return(omegas * xxmax);
		case 'd': /* relative displacement */
			return(xxmax);
	}
	return(-1.);
}

/* multiply convolution output by
	appropriate factor to get desired response.
	If this is a zero-period call, assume we have
	the max of the acceleration */
static double scalec(double tau, double lambda, double dt, double max_disp, char type) {
	double max_val = max_disp;
	double tmp = 2. * M_PI * sqrt(1. - lambda * lambda);
	double c = tau * dt / tmp;

	if ( type == 'a' ) {
		if ( tau > 0.00001 ) {
			max_val *= c;
			max_val *= 4. * M_PI * M_PI / ( tau * tau );
		}
	} else if ( type == 'v' ) {
		if ( tau > 0.00001 ) {
			max_val *= 2. * M_PI * c / tau;
		} else {
			max_val = 0.;
		}
	} else {
		if ( tau > 0.00001 ) {
			max_val *= c;
		} else {
			max_val = 0.;
		}
	}

	return(max_val);
}

/* CPEAK - calculate a response spectrum using convolution */
double cpeak(double tau, double lambda, double *ts, int len, double dt, char type)
{
	int ii, kk, jj, jlo, jhi, imax, lenosc;
	double kw, kwp, *osc;
	double tmp, max;
	double sum;

	if ( tau < 0.00001 ) { /* zero-period case */
		/* find max - skip first few points */
		max = 0.; imax = 0;
		for ( ii = 6 ; ii < len ; ii++ ) {
			tmp = ABS(ts[ii]);
			if ( tmp > max ) {
				max = tmp;
				imax = ii;
			}
		}
		return(scalec(tau, lambda, dt, max, type));
	}

	if ( lambda > 1 )
		oops("cpeak","fractional damping greater than 1.0");

	/* get constants for current tau */
	kw = 2. * M_PI * dt * lambda / tau;
	kwp = 2. * M_PI * dt * sqrt(1. - lambda * lambda);
	kwp /= tau;

	/* get oscillator length.  For small tau a length shorter than
	   the time series may be found.  Criterion is point where
	   exp is down by factor of 10^3.  Otherwise, just use len */
	lenosc = MIN(ROUND(3 * log(10.) / kw + 1), len);

	if ( (osc = calloc(lenosc, sizeof(*osc))) == NULL )
		oops("cpeak","can's get space for oscillaror");

	/* fill oscillator buffer */
	for ( ii = 0 ; ii < lenosc ; ii++ )
		osc[ii] = exp(-kw * ii) * sin(kwp * ii);

	/* Do convolution only out to time len * dt since oscillator
	   and time series are decreasing in amplitude and max will
	   occur before then.  For large tau, this won't work */
	max = 0.; imax = 0;
	for ( kk = 0 ; kk < len ; kk++ ) {
		sum = 0.;
		/* the following scheme assumes lenosc <= len */
		jlo = ( kk > lenosc ? kk - lenosc : 0 );
		jhi = ( kk < len ? kk : len - 1 );
		for ( jj = jlo ; jj <= jhi ; jj++ ) {
			ii = kk - jj;
			sum += ts[jj] * osc[ii];
		}
		/* check for max - skip first few points */
		if ( kk > 6 )
			chkmax(sum,&max,kk,&imax);
		/* rube exit */
		if ( kk > len / 2 && kk > 2 * imax )
			break;
	}

	return(scalec(tau, lambda, dt, max, type));
}


#ifdef HAVE_FFTW3
// return the modulus of a DCMPLX value
static double modulus( DCMPLX *x ) {
	return sqrt( x->r*x->r + x->i*x->i );
}

// return the phase of a DCMPLX value
static double phase( DCMPLX *x ) {
	return atan2( x->i, x->r );
}

// fftw3 wrapper
static CMPLX* fft(DCMPLX *x, int N) {
	int ii;
	fftw_plan p = fftw_plan_dft_1d(N, (fftw_complex *)x, (fftw_complex *)x, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	// convert in-place to CMPLX
	for ( ii = 0 ; ii < N ; ii++ ) {
		CMPLX ctmp;
		ctmp.a = modulus(x + ii);
		ctmp.p = phase(x + ii);
		memcpy( x + ii, &ctmp, sizeof( *x ) );
	}
	return (CMPLX *)x;
}

// fftw3 wrapper
static DCMPLX* ifft(CMPLX *x, int N) {
	int ii;
	// convert in-place to CMPLX, and scale by N for inverse transform
	for ( ii = 0 ; ii < N ; ii++ ) {
		DCMPLX dtmp;
		dtmp.r = x[ii].a * cos( x[ii].p ) / N;
		dtmp.i = x[ii].a * sin( x[ii].p ) / N;
		memcpy( x + ii, &dtmp, sizeof( *x ) );
	}
	fftw_plan p = fftw_plan_dft_1d(N, (fftw_complex *)x, (fftw_complex *)x, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	return (DCMPLX *)x;
}

/* fft_int - take the fft of a real time series, multiply/divide by (i*omega)^n,
	and take the inverse fft */
void fft_int(double *ts, double dt, int len, int nn)
{
	int ii;
	double cnst, df;
	DCMPLX *sbuf, *cts;
	CMPLX* s;

	if ( nn == 0 )
		return;

	df = 1 / (len * dt);

	sbuf = calloc( len, sizeof( *sbuf ) );
	for ( ii = 0 ; ii < len ; ii++ )
		sbuf[ii].r = ts[ii];
	s = fft( sbuf, len );

	/* multiply spectrum by (i * omega)^nn */
	cnst = 2. * M_PI * df;
	for ( ii = 0 ; ii <= len/2 ; ii++ ) {
		int jj;
		double amp, phase, om;
		om = cnst * ii;
		if ( ii == 0 && nn < 0 )
			continue;
		amp = ( nn < 0 ? 1. / om : om);
		phase = ( nn < 0 ? -1. : 1. ) * M_PI / 2.;
		for ( jj = 0 ; jj < ABS(nn) ; jj++ ) {
				/* positive frequencies */
			s[ii].a *= amp;
			s[ii].p += phase;
				/* negative frequencies */
			if ( ii == 0 || (ii == len/2 && len % 2 == 0) )
				continue;
			s[len-ii].a *= amp;
			s[len-ii].p -= phase;
		}
	}

	// take the inverse FFT, and replace input with real part
	cts = ifft( s, len );
	for ( ii = 0 ; ii < len ; ii++ )
		ts[ii] = cts[ii].r;

}

/* multiply FFT output by
	appropriate factor to get desired response.
	If this is a zero-period call, assume we have
	the max of the acceleration */
static double scalef(double tau, double max_disp, double f0, char type) {
	double max_val = max_disp;
	double c = 2. * M_PI / tau;
	c = 2. * M_PI * f0;
	switch ( type ) {
		case 'a': /* pseudo absolute acceleration */
			break;
		case 'v': /* pseudo relative velocity */
			if ( tau > 0. )
				max_val /= c;
			else
				max_val = 0.;
			break;
		case 'd': /* relative displacement */
			if ( tau > 0. )
				max_val /= c * c;
			else
				max_val = 0;
			break;
	}

	return(max_val);
}


/* FPEAK - get peak response for a particular natural frequency using
	FFT method */
double fpeak(double tau, double lambda, double *ts, int len, double dt, char type)
{
	int ii, imax;
	double df, f0, max;
	DCMPLX *sbuf, *cts;
	CMPLX *s;

	df = 1./ ( len * dt );

	sbuf = calloc( len, sizeof( *sbuf ) );
	for ( ii = 0 ; ii < len ; ii++ )
		sbuf[ii].r = ts[ii];

	s = fft(sbuf, len);	// in-place fft, so just a pointer to sbuf

	/* multiply previously transformed acceleration time series
		by oscillator response to get displacement spectrum */

	if ( tau > 0. ) {
		f0 = 1. / tau;
		for ( ii = 0 ; ii <= len/2 ; ii++) {
			int jj;
			double amp, phase, freq, f1, f2;
			freq = ii * df;
			f1 = freq / f0;
			f2 = f1 * f1;
			amp = f2 - 1.; amp *= amp;
			amp += 4. * lambda * lambda * f2;
			amp = 1. / sqrt(amp);
			phase = atan2(2. * lambda * f1,f2 - 1.);
			/* positive frequencies */
			s[ii].a *= amp;
			s[ii].p += phase;
			/* negative frequencies */
			if ( ii == 0 || ii == len/2 )
				continue;
			jj = len - ii;
			s[jj].a *= amp;
			s[jj].p -= phase;
		}
	} else {
		/* zero-period case */
		f0 = -1.;
	}

	/* inverse transform */
	cts = ifft(s, len);	// in-place fft, so just another pointer to sbuf

	/* find the peak absolute value of the time series, which is
		real. Skip the first few points */
	max = 0; imax = 0;
	for ( ii = 6 ; ii < len ; ++ii)
		chkmax(ABS(cts[ii].r),&max,ii,&imax);

	free( sbuf );

	/* now have relative displacement - scale and return */
	return(scalef(tau, max, f0, type));
}
#endif /* HAVE_FFTW3 */
