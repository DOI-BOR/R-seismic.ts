#' Minimum-phase representation of an input time series
#'
#' \code{minimum.phase} returns the minimum-phase representation of an
#' input univariate time series
#'
#' @param x.data Equally-sampled input series. Must convert to numeric vector.
#' @param zero.pad Zero-pad the input using \code{\link{zero.pad}} before
#' transforming. Default is TRUE.
#'
#' @details Uses a clever method...
#' @return The minimum-phase representation.
#' @seealso \code{\link{zero.pad}}
#'
#' @keywords ts

minimum.phase <- function(x.data, zero.pad=TRUE) {
	if ( zero.pad )
		x <- zero.pad(x.data)
	else
		x <- as.vector(x.data)
	x.fft <- fft(x) / length(x)
	x.aspec <- abs(x.fft)
	x.min.phase <- hilbert(log(x.aspec), "hilbert", zero.pad)
	x.min.phase.fft <- x.aspec * exp((0+1i) * x.min.phase)
	# tsplot(x.min.phase)
	x.min.phase.ifft <- fft(x.min.phase.fft, inverse=TRUE)
	x.min.phase <- Re(x.min.phase.ifft)
	return (x.min.phase)
}
