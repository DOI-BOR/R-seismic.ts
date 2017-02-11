#' Zero-pad a univariate time series
#'
#' \code{zero.pad} adds zeros to the end of an input timte series. The number
#' of added zeros can be specified, or, by default, enough zeros are added to
#' make the length a power of 2.
#'
#' @param x.data Equally-sampled input series. Must convert to numeric vector.
#' @param padded.len Length of input series after zero-padding. Default is to pad
#' the input until the length is the next power of 2.
#' @param extra If padded.len is not set, then pad the input data to the
#' next-next power of 2.
#' @return the zero-padded data.

#' @keywords ts

zero.pad <- function(x.data, padded.len = NA, extra = FALSE) {
	x <- as.vector(x.data)
	n = length(x)
	if ( is.na(padded.len) ) {
		n2 = 2
		while ( n2 < n )
			n2 <- 2 * n2
		if ( extra )
			n2 <- 2 * n2
		x.aug = c(x, rep(0, n2 - n))
	} else {
		npad = as.integer(round(padded.len)) - n
		if ( npad > 0 )
			x.aug = c(x, rep(0, npad))
		else
			x.aug = x
	}
	return (x.aug)
}
