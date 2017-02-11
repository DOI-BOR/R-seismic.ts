#' Determine the prediction-error filter
#'
#' \code{predictionFilter} finds the prediction-error filter for a specified
#' univariate input time series.
#'
#' @param x.data Equally-sampled input series. Must convert to numeric vector.
#' @param length.filt Length of filter (default is 0.2 the input length)
#' @param lag Lag used for determing the filter (default is 1)
#' @param length.wavelet Lag used for determining the inverse of the prediction-
#' error filter
#'
#' @details What this all means....
#' @return List with the prediction-error filter, innovations, wavelet, and
#' various measures.
#' @seealso \code{\link{acf}}, \code{\link{predictionFilter}},
#' \code{\link{wienerFilter}}
#'
#' @keywords ts

predictionErrorFilter <- function(x.data, length.filt = NULL, lag = NULL, length.wavelet = NA) {
	pred.filt <- ilsd::predictionFilter(x.data, length.filt, lag)

	# get prediction error filter
	if ( pred.filt$lag > 1 )
		perr.filt <- c(1, rep(0, lag - 1), -pred.filt$filter)
	else
		perr.filt <- c(1, -pred.filt$filter)

	# get the inovations
	ino <- ilsd::my.filter(pred.filt$x, perr.filt)

	# if lag = 1, get the wavelet, which is the inverse of the (assumed minimum delay)
	# prediction error filter
	wavelet <- NA
	if ( pred.filt$lag == 1 )
		wavelet <- inverseWavelet(perr.filt, length.wavelet)

	results = list( filter = perr.filt, lag = pred.filt$lag,
									min.error = pred.filt$min.error, est.min.error = pred.filt$est.min.error,
									ino = ino, wavelet = wavelet, x.mean = pred.filt$x.mean )
	return ( results )
}
