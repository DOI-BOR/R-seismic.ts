#' Hilbert Transform of a time series using analytic signal method
#'
#' \code{hilbert} returns the Hilbert transform of a univariate
#' time series using the analytic signal approach.
#'
#' @param x.data Equally-sampled input series. Must convert to numeric vector.
#' @param op Operation to perform. One of
#' \code{c("hilbert", "envelope", "instphase", "instfreq")}. Default is \code{"hilbert"}
#' @param zero.pad Zero-pad the input using \code{\link{zero.pad}} before
#' transforming. Default is TRUE.
#' @param dt Sample interval, in seconds. Default is 0.01.
#' @details This function may be used to determine the Hilbert transform, envelope
#' function, instantaneous phase, or instantaneous frequency of an input
#' time series.
#' @return the transform of the data.
#' @seealso \code{\link{analytic.ts}}, \code{\link{unwrap.phase}},
#' \code{\link{zero.pad}}, \code{\link{hilbert}},
#' \href{https://en.wikipedia.org/wiki/Analytic_signal}{Wikipedia}
#' @keywords ts

hilbert <- function(x.data, op = "hilbert", zero.pad = TRUE, dt = 0.01) {
	# note: length(x) does not have to be a power of 2
	if ( zero.pad )
		x <- zero.pad(x.data)
	else
		x <- as.vector(x.data)
	xa <- analytic.ts(x)
	result <- NULL
	if ( op == "hilbert" ) {
		result <- Im(xa)
	} else if ( op == "envelope" ) {
		result <- Mod(xa)
	} else if ( op == "instphase" ) {
		inst.phase <- atan2(Im(xa),Re(xa))
		result <- unwrap.phase(inst.phase)
	} else if ( op == "instfreq" ) {
		inst.phase <- atan2(Im(xa),Re(xa))
		result <- diff(unwrap.phase(inst.phase),1) / (2 * pi * dt)
	}
	return (result)
}
