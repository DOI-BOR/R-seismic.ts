#' Decimate an Augmented Time Series
#'
#' @description
#' \code{decimate} down-samples a real, univariate time series that has
#' been augmented, or suitably low-pass filtered.
#'
#' @param x.data Equally-sampled input series, signalSeries, or vector.
#' Must convert to a numeric vector.
#' @param factor Decimation factor (default is 1). If factor = N,
#' then dt -> dt * N
#' @return The decimated time series.
#' @details low-pass filtering is assumed to not be needed because the input
#' time series was augmented, or has been suitably filtered. This function
#' should not be used to decimate a time series that has non-zero spectral
#' amplitudes for frequencies above 1 / (dt * factor).
#' @seealso \code{\link{signalSeries}}
#' @keywords ts
decimate <- function(x.data, factor = NA) {
	if ( is.na(factor) )
		factor = 1.;
	skip <- as.integer(round(factor))
	if ( factor < 1. ) {
		x.dec <- ilsd::augment(x.data, 1./factor)
	} else if ( skip > 1 ) {
		if ( is(x.data, "series") || is(x.data, "signalSeries") )
			x <- x.data
		else
			x <- splus2R::signalSeries(as.vector(x.data))
		x.len <- length(x)
		x.dec <- splus2R::signalSeries(x@data[seq(1,x.len,skip)],from=x@positions@from,by=(x@positions@by*skip))
	} else {
		x.dec <- x.data
	}
	return ( x.dec )
}
