#' Compute the derivative of a time series using the finite difference method.
#'
#' \code{deriv.fd} is used to compute the derivative of a numerical vector
#' using the finite difference method.
#'
#' @param xt required (actual) equally-sampled input series. Must convert to numeric vector.
#' @param dt (optional) sample interval, in seconds. Default is 0.01.
#' @param nd (optional) integer order of the derivative. Currently only first and
#' second derivatives are supported. Default is 1.
#' @param order (optional) integer order of the finite difference. Currently only
#' orders of 2, 4, 6, and 8 are supported. Default is 8.
#' @param pct (optional) percentage of data window to apply a Hanning taper. Must be
#' between 0 and 50. Default is 0.
#' @return the derivative of the windowed data
deriv.fd <- function(xt, dt=0.01, nd=1, order=8, pct=NA) {

	xt <- as.double(xt[!is.na(xt)])
	len <- length(xt)
	if ( len < 3 )
		stop("input time series must have at least 3 valid points")

	order <- if ( order == 8 ) 3 else if ( order == 6 ) 2 else if ( order == 4 ) 1 else 0

	# call C function
	out <- .Call("CALLfd_deriv",
							 as.double(xt), as.double(dt), as.integer(nd),
							 as.integer(order), as.double(pct))
	return(out)
}
