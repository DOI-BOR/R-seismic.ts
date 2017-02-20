#' Compute the Analytic Signal
#'
#' @description
#' \code{analytic.ts} computes the analytic signal for a real,
#' univariate or multivariate time series.
#'
#' @param x.data Equally-sampled univariate or multivariate series. Supported
#' types include a numeric \code{\link{vector}}, \code{\link{matrix}},
#' \code{\link{data.frame}}, \code{\link{ts}}, or \code{\link{signalSeries}}.
#' @details The analytic signal s(t) of a function f(t) is defined
#' by s(t) = f(t) - i * Hilbert(f(t)), where H(f(t)) is the Hilbert
#' transform of f(t). This code uses the Fourier transform to compute
#' the anatic signal, using the property that s(t) <=> 2 * H(omega) * F(omega),
#' where H() is the Heaviside function, and f(t) <=> F(omega).
#' @return The analytic signal of the input time series.
#' @seealso \href{https://en.wikipedia.org/wiki/Analytic_signal}{Wikipedia}
#' @keywords ts
#'
analytic.ts <- function(xt) {
  # check for multi-trace input
  multi.trace <- is.mts(xt) || is.matrix(xt) || length(dim(xt)) > 1 ||
    ( is(xt, "signalSeries") && ! is.null(dim(xt)) )

  if ( multi.trace == TRUE ) {
    if ( is(xt, "signalSeries") )
      x.len <- dim(xt)[1]
    else
      x.len <- dim(xt)[1]
  } else
    x.len <- length(xt)

	# make the Heaviside function times 2. For fft(), note that at discontinuities, h(f) -> 0.5 * [h(f-) + h(f+)],
	#	and that the indices are defined as follows:
	#	n even
	#		1 = zero frequency (discontinuity)
	#		2:n/2 = positive frequencies
	#		n/2+1 = nyquist (discontinuity)
	#		n/2+2:n = negative frequencies
	#	n odd
	#		1 = zero frequency (discontinuity)
	#		2:(n+1)/2 = positive frequencies
	#		(n+1)/2+1:n = negative frequencies
	n <- x.len
	Hf <- rep(0,n) # initialize all values to zero
	if ( n == 2 * floor(n/2) ) {
		# n is even
		Hf[c(1, n/2 + 1)] = 1	# zero and Nyquist frequencies - disconinuity value = 0.5*(0+2) = 1
		Hf[2:n/2] = 2	# positive frequencies
	} else {
		# n is odd
		Hf[1] = 1 # zero frequency - discontinuity value = 1
		Hf[2:(n + 1)/2] = 2	# positive frequencies
	}

	if ( multi.trace == TRUE ) {
		xa <- xt
		for ( cn in 1:dim(xt)[2] ) {
			x <- as.vector(xt[,cn])
			if ( is(xt, "signalSeries") )
			  x <- as.vector(xt[,cn]@data)
			else
			  x <- as.vector(xt[,cn])
			Xf <- fft(x) / n   # if length(x) is prime, the fft is just the dft
			# multiply dft of input by Heaviside function, and take inverse transform
			xa[,cn] <- fft(Xf * Hf, inverse=TRUE)
		}
	} else {
		# get dft of the input data
	  if ( is(xt, "signalSeries") )
	    x <- as.vector(xt@data)
	  else
	    x <- as.vector(xt)
	  Xf <- fft(x) / n   # if length(x) is prime, the fft is just the dft
		# multiply dft of input by Heaviside function, and take inverse transform
		xa <- fft(Xf * Hf, inverse=TRUE)
	}

	return (xa)
}
