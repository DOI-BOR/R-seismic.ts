# Augment a time series (i.e., sample at a higher rate without
# changing the frequency content)
# Required Arguments:
#		x.data:
#			series, signalSeries, or vector data
# Optional Arguments
#		factor:
#			augmentation factor (>= 1). If factor = N, then the
#			sample rate dt goes to dt/N. The frequency content
#			is unchanged; the time series is simply interpolated to
#			a higher sample rate.

augment <- function(x.data, factor = NA) {
	if ( is.na(factor) )
		factor = 1.;
	if ( factor < 1. ) {
		# decimate, rather than augment. Data must already be suitably low-pass filtered
		x.aug <- ilsd::decimate(x.data, 1./factor)
	} else if ( factor > 1. ) {
		# augment by zero-padding the FFT at frequencies greater than the Nyquist, and then taking the inverse transform
		if ( is(x.data, "series") || is(x.data, "signalSeries") )
			x <- x.data
		else
			x <- splus2R::signalSeries(as.vector(x.data))
		x.len <- length(x)
		x.fft <- fft(x)
		x.fft.len <- length(x.fft)
		x.fft <- x.fft / x.fft.len
		z.len <- as.integer(round((factor - 1.)*x.fft.len))
		if ( x.fft.len == 2 * floor(x.fft.len/2) ) {
			# n is even
			x.fft.aug <- as.vector(c( x.fft[1:(x.fft.len/2 + 1)], rep(0+0i,z.len), x.fft[(x.fft.len/2 + 2):x.fft.len] ))
		} else {
			# n is odd
			x.fft.aug <- as.vector(c( x.fft[1:((x.fft.len + 1)/2)], rep(0+0i,z.len), x.fft[((x.fft.len + 1)/2 + 1):x.fft.len] ))
		}
		x.aug <- splus2R::signalSeries(Re(fft(x.fft.aug,inverse=T)), from=x@positions@from, by=x@positions@by/factor)
	} else {
		# factor is 1, so no need to do anything
		x.aug <- x.data
	}
	return ( x.aug )
}
