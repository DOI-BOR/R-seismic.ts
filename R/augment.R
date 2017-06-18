#' Augment a Time Series
#'
#' @description
#' \code{augment} augments a real, univariate time series. Augmentation is
#' defined as resampling a time series at a higher sample rate, without
#' changing the frequency content.
#'
#' @param xt Equally-sampled input series. Must convert to a numeric vector.
#' @param factor Augmentation factor > 0. If factor = N > 1, then the
#' sample rate dt goes to dt/N. The frequency content is unchanged;
#' the time series is simply interpolated to a higher sample rate. If
#' factor = 1, then the time series is unchanged. If factor < 0, then
#' the time series is decimated by 1/factor. Default is 1.
#' @return The augmented time series.
#' @details For factor > 1, this code augments by zero-padding the Fourier
#' transform of the input series at frequencies greater than the Nyquist,
#' and then taking the inverse transform.
#' @keywords ts

augment <- function(xt, factor=NA) {
	if ( is.na(factor) )
		factor = 1.;
	if ( factor < 1. ) {
		# decimate, rather than augment. Data must already be suitably low-pass filtered
		x.aug <- decimate(xt, 1./factor)
	} else if ( factor > 1. ) {
		# augment by zero-padding the FFT at frequencies greater than the Nyquist, and
	  # then taking the inverse transform
		if ( is(xt, "series") || is(xt, "signalSeries") )
			x <- xt
		else
			x <- splus2R::signalSeries(as.vector(xt))
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
		x.aug <- xt
	}
	return ( x.aug )
}
