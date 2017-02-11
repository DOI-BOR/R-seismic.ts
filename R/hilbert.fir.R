#' Hilbert Transform of a time series using FIR filter method
#'
#' \code{hilbert.fir} returns the Hilbert transform of a univariate
#' time series using using a FIR filter approach.
#'
#' @param xt Equally-sampled input series. Must convert to numeric vector.
#' @param dt Sample interval, in seconds. Default is 0.01.
#'
#' @return the Hilbert transform of the data
#' @seealso \code{\link{hilbert}}
#' @keywords ts

hilbert.fir <- function(xt, zero.pad = TRUE, dt = 0.01) {

	xt <- as.double(xt[!is.na(xt)])
	# note: length(x) does not have to be a power of 2
	if ( zero.pad )
		x <- zero.pad(xt)
	else
		x <- as.vector(xt)
	len <- length(xt)
	if ( len < 3 )
		stop("input time series must have at least 3 valid points")

	# call C function
	out <- .Call("CALLhilbertr_fir", as.double(x), as.double(dt))
	return(out)
}
