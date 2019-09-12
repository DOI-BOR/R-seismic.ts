#' Augment a Time Series
#'
#' @description
#' \code{augment} augments a real, univariate time series. Augmentation is
#' defined as resampling a time series at a higher sample rate, without
#' changing the frequency content.
#'
#' @param xt Equally-sampled input series
#' represented as a numeric \code{\link{vector}}, \code{\link{matrix}},
#' \code{\link{data.frame}}, \code{\link{ts}}, \code{\link{mts}}, or
#' \code{\link{signalSeries}} object.
#' @param factor Augmentation factor > 0. If factor = N > 1, then the
#' sample rate dt goes to dt/N. The frequency content is unchanged;
#' the time series is simply interpolated to a higher sample rate. If
#' factor = 1, then the time series is unchanged. If factor < 0, then
#' the time series is decimated by 1/factor. Default is 1.
#' @return The augmented time series, in the same representation as the input.
#' @details For factor > 1, this code augments by zero-padding the Fourier
#' transform of the input series at frequencies greater than the Nyquist,
#' and then taking the inverse transform.
#' @keywords ts

augment <- function(xt, factor=NA) {
	if ( is.na(factor) )
		factor = 1.;
	xt.aug <- NULL
	if ( factor < 1. ) {
		# decimate, rather than augment. Data must already be suitably low-pass filtered
		xt.aug <- decimate(xt, 1./factor)
	} else if ( factor > 1. ) {
		# augment by zero-padding the FFT at frequencies greater than the Nyquist, and
	  # then taking the inverse transform

    multi.trace <- is.mts(xt) ||
	    ( ! is.null(dim(xt)) && length(dim(xt)) > 1 && dim(xt)[2] > 1 )

    zt.len <- NULL
    if ( multi.trace ) {
      zt.len <- dim(xt)[1]
      PLAN <- fftw::planFFT(zt.len)
      aug <- NULL
      cnames <- NULL
	    for ( ii in 1:dim(xt)[2] ) {
	      if ( is(xt, "signalSeries") ) {
	        zt <- xt@data[,ii]
	        cnames <- colnames(xt@data)
	      } else {
	        zt <- xt[,ii]
	        cnames <- colnames(xt)
	      }

        # take FFT and scale by length
        zt.fft <- FFT(zt, plan=PLAN)
        zt.fft.len <- length(zt.fft)
        zt.fft <- zt.fft / zt.fft.len

        # zero-pad FFT
        zero.len <- as.integer(round((factor - 1.)*zt.fft.len))
        if ( zt.fft.len == 2 * floor(zt.fft.len/2) ) {
          # FFT length is even
          zt.fft.aug <- as.vector(c( zt.fft[1:(zt.fft.len/2 + 1)],
                                     rep(0+0i,zero.len),
                                     zt.fft[(zt.fft.len/2 + 2):zt.fft.len] ))
        } else {
          # FFT length is odd
          zt.fft.aug <- as.vector(c( zt.fft[1:((zt.fft.len + 1)/2)],
                                     rep(0+0i,zero.len),
                                     zt.fft[((zt.fft.len + 1)/2 + 1):zt.fft.len] ))
        }
        zt.aug.len <- length(zt.fft.aug)

        # add real part of inverse transform to data frame
        PLAN.aug <- fftw::planFFT(zt.aug.len)
        zt.aug <- Re(IFFT(zt.fft.aug, PLAN.aug, scale=FALSE))
	      if ( is.null(aug) )
	        aug <- zt.aug
	      else
	        aug <- data.frame(aug, zt.aug)
	    }
	    colnames(aug) <- cnames
	  } else {
	    if ( is(xt, "signalSeries") ) {
	      if ( ! is.null(dim(xt)) && length(dim(xt)) > 1 && dim(xt)[2] == 1 )
	        zt <- xt@data[,1]
	      else
	        zt <- xt@data
	    } else
	      zt <- xt
      zt.len <- length(zt)
      PLAN <- fftw::planFFT(zt.len)

      # take FFT and scale by length
      zt.fft <- FFT(zt, plan=PLAN)
      zt.fft.len <- length(zt.fft)
      zt.fft <- zt.fft / zt.fft.len

      # zero-pad
      zero.len <- as.integer(round((factor - 1.)*zt.fft.len))
      if ( zt.fft.len == 2 * floor(zt.fft.len/2) ) {
        # FFT length is even
        zt.fft.aug <- as.vector(c( zt.fft[1:(zt.fft.len/2 + 1)],
                                   rep(0+0i,zero.len),
                                   zt.fft[(zt.fft.len/2 + 2):zt.fft.len] ))
      } else {
        # FFT length is odd
        zt.fft.aug <- as.vector(c( zt.fft[1:((zt.fft.len + 1)/2)],
                                   rep(0+0i,zero.len),
                                   zt.fft[((zt.fft.len + 1)/2 + 1):zt.fft.len] ))
      }
      zt.aug.len <- length(zt.fft.aug)

      # add real part of inverse transform to data frame
      PLAN.aug <- fftw::planFFT(zt.aug.len)
      zt.aug <- Re(IFFT(zt.fft.aug, PLAN.aug, scale=FALSE))
      aug <- zt.aug
	  }
    # get augmented sample increment
    if ( is(xt, "signalSeries") || is(xt, "ts") )
      dt <- deltat(xt)
    if ( ! is.finite(dt) )
      dt <- 1
    dt.aug <- dt * zt.len / zt.aug.len
	  # create a new object with same type as input
	  if ( is(xt, "signalSeries") ) {
	    xt.aug <- splus2R::signalSeries(aug, from=xt@positions@from, by=dt.aug,
	                                    units=xt@units,
	                                    units.position=xt@units.position)
	  } else if ( is(xt, "ts") ) {
	    xt.aug <- ts(aug, start=start(xt)[1], deltat=dt)
	  } else {
	    xt.aug <- aug
	  }
	} else {
		# factor is 1, so no need to do anything
		xt.aug <- xt
	}
	return ( xt.aug )
}
