#' Direct inverse of a time series
#'
#' \code{inverseWavelet.} finds the inverse of a
#' minimum-phase univariate time series
#'
#' @param wavelet Equally-sampled input series. Must convert to numeric vector.
#' @param inv.len Length of inverse (default is the input length)
#'
#' @details Directly determines the inverse through division, but the
#' inverse will diverge unless the input is minimum phase. For signals
#' that are not minimum-phase, use \code{\link{inverseWavelet.ls}}.
#' @return List with the inverse and various measures.
#' @seealso \code{\link{inverseWavelet.ls}}, \code{\link{minimum.phase}}
#'
#' @keywords ts

inverseWavelet <- function(wavelet, inv.len = NA) {
	bt <- as.vector(wavelet)
	bt.len <- length(bt)
	if ( is.na(inv.len) )
		inv.len <- bt.len
	at.len <- inv.len
	at <- 1. / bt[1]
	for ( tt in seq(from=2, by=1, length=(at.len-1)) ) {
		if ( tt > bt.len ) {
			btmax <- bt.len
		} else {
			btmax <- tt
		}
		at <- c(at, -sum(at[(tt-1):(tt-btmax+1)]*bt[2:btmax])/bt[1])
	}
	at
}
