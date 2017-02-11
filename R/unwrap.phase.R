#' Unwrap a phase time series.
#'
#' \code{unwrap.phase} replaces jumps of greater than pi with a jump of
#' phase - 2 * pi or pase + 2 * pi.
#'
#' @param x.data required (first) equally-spaced phases Must convert to
#' numeric vector.
#' @return vector with unwrapped phase
#' @keywords ts
unwrap.phase <- function(x.data) {
	x <- as.vector(x.data)
	max.jump <- pi
	x.last <- x[1]
	phase <- 0
	xu <- NULL
	for ( x.cur in x ) {
		x.diff <- x.cur - x.last
		x.last <- x.cur
		if ( x.diff > max.jump ) {
			phase <- phase - 2 * pi
		} else if ( x.diff < -max.jump ) {
			phase <- phase + 2 * pi
		}
		xu <- c(xu, x.cur + phase)
	}
	return(xu)
}
