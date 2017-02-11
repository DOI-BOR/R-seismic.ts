#' Determine a prediction filter
#'
#' \code{predictionFilter} finds the optimal filter, in the least-squares sense,
#' that transforms a specified univariate input time series into a lagged
#' version of itself.
#'
#' @param x.data Equally-sampled input series. Must convert to numeric vector.
#' @param length.filt Length of filter (default is 0.2 the input length)
#' @param lag Lag used for determing the filter (default is 0)
#'
#' @details Uses \code{\link{acf}} to form the normal equations, which by default
#' de-means the series, so if the desired output z has non-zero mean, the output
#' filter will be for de-meaned z. The returned value of estimated z will add the
#' mean back in. For long time series, we could just correct for the means, but
#' this won't work for short wavelets since the time series are not assumed to
#' repeat.
#' @return List with the Wiener filter and various measures.
#' @seealso \code{\link{acf}}, \code{\link{wienerFilter}}
#'
#' @keywords ts

predictionFilter <- function(x.data, length.filt = NULL, lag = NULL) {
	x <- as.vector( x.data )

	# get filter length; default is 10% of input series length
	x.data.len <- length(x)
	f.len = length.filt
	if ( is.null(length.filt) )
		f.len <- as.integer(max( min(2, x.data.len), 0.1 * x.data.len))
	if ( is.null(lag) )
		lag <- 1

	# de-mean the input series and put into a signalSeries
	x.mean <- mean(x)
	x <- x - x.mean
	xs <- splus2R::signalSeries(x)

	# get the auto covariance
	x.acvf <- acf(xs, lag.max = f.len + lag - 1, type="covariance", plot=F)
	axx <- x.acvf$acf[1:f.len,1,1] # + x.mean * x.mean

	# create the autocovariance matrix from the acvf
	axx.m <- toeplitz(axx)

	# get inverse using generic solver
	#axx.m.inv <- solve(axx.m)

	# get inverse using toeplitz solver
	axx.m.inv <- ltsa::TrenchInverse(axx.m)

	# get the lagged autocovariance
	axx.lag <- x.acvf$acf[(lag+1):(lag+f.len),1,1] # + x.mean * x.mean

	# multiply the inverse of the autocovariance matrix with the lagged autocovariance to get the prediction filter coefficients
	ls.filter <- axx.m.inv %*% axx.lag

	# get the normalized error
	min.error <- ilsd::minNMSE(axx, axx.lag, ls.filter)

	# get predicted output (first lag points are zero)
	x.est <- c(rep(0,lag), ilsd::my.filter(x, ls.filter) + x.mean)

	# directly calculate the normalized error (skip lag)
	est.min.error <- ilsd::estMinNMSE(x[(1+lag):x.data.len] + x.mean, x.est[(1+lag):x.data.len], x.acvf$acf[1,1,1])

	# return a list of the results
	results = list( filter = ls.filter, lag = lag, x = x, x.est = x.est,
									min.error = min.error, est.min.error = est.min.error,
									x.acvf = x.acvf, axx = axx, axx.lag = axx.lag, x.mean = x.mean )
	return ( results )
}
