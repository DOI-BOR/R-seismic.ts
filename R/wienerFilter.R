#' Determine a Wiener (least-squares) filter
#'
#' \code{wienerFilter} finds the optimal filter, in the least-squares sense,
#' that transforms a specified univariate input time series to a specified
#' univariate target output time series.
#'
#' @param x.data Equally-sampled input series. Must convert to numeric vector.
#' @param z.data Equally-sampled target output series. Must convert to numeric vector.
#' @param length.filt Length of filter (default is 0.2 the input length)
#' @param lag Lag used for determing the filter (default is 0)
#' @param demean Should windowed data be demeaned? Default is FALSE.
#'
#' @details Uses \code{\link{acf}}, which by default de-means the series, so if
#' the desired output z has non-zero mean, the output filter will be for
#' de-meaned z. The returned value of estimated z will add the mean back in.
#' For long time series, we could just correct for the means, but this won't
#' work for short wavelets since the time series are not assumed to repeat.
#' @return List with the Wiener filter and various measures.
#' @seealso \code{\link{acf}}
#'
#' @keywords ts

wienerFilter <- function(x.data, z.data, length.filt = NULL, lag = NULL,
                         demean = FALSE) {

	# check the lag value
	if ( is.null(lag) )
		lag <- 0
	if ( lag < 0 )
		stop ("negative lags are not supported")

	x <- as.vector(x.data)
	z <- as.vector(z.data)

	# get filter length; default is 20% of input series length
	x.data.len <- length(x)
	f.len = length.filt
	if ( is.null(length.filt) )
		f.len <- max( min(2, x.data.len), 0.2 * x.data.len)

	# handle case where lag exceeds the filter length
	diff <- lag - f.len
	max.acf.lag <- f.len - 1 + max(0, diff)

	# acf needs equal-length x and z, so get amount to truncate input (x) series, if necessary
	z.data.len <- length(z)
	x.data.trim <- 0
	if ( z.data.len < x.data.len + max.acf.lag ) {
		x.data.trim <- x.data.len + max.acf.lag - z.data.len
		if ( is.null(length.filt) )
			f.len <- max( min(2, x.data.len - x.data.trim), 0.2 * (x.data.len - x.data.trim))
	}

	# acf needs equal-length x and z, so get amount to pad input series so as to match output length
	x.data.pad <- z.data.len - x.data.len + x.data.trim

	# trim, de-mean, and then pad the input series (pad only after de-meaning)
	x <- x[1:(x.data.len - x.data.trim)]
	x.mean <- mean(x)
	if ( demean )
		x <- x - x.mean
	if ( x.data.pad > 0 )
		x <- as.vector( c(x, rep(0, x.data.pad)) )
	z.mean <- mean(z)
	if ( demean )
		z <- z - z.mean

	# bind the input (x) and output (z) series into a multivariate time series
	#		of length equal to the output length
	xz <- ts(cbind(x,z))

	# get the auto and cross covariances. normally, we only need the first f.len - 1 values of acf
	xz.acvf <- acf(xz, lag.max = max.acf.lag, type="covariance", plot=FALSE)

	axx <- xz.acvf$acf[1:f.len,1,1] # + x.mean * x.mean

	azz <- xz.acvf$acf[1:f.len,2,2] # + z.mean * z.mean
	if ( lag == 0 ) {
		# all points are from zx
		azx <- xz.acvf$acf[1:f.len,2,1] # + x.mean * z.mean
	} else if ( lag < f.len ) {
		# prepend points from xz. note: xz[1] = zx[1], so skip over first xz point
		azx <- c(xz.acvf$acf[(lag+1):(1+1),1,2], xz.acvf$acf[1:(f.len-lag),2,1]) # + x.mean * z.mean
	} else {
		# all points are from xz
		azx <- xz.acvf$acf[(lag+1+diff):(1+1+diff),1,2] # + x.mean * z.mean
	}

	# create the autocovariance matrix from the acvf
	axx.m <- toeplitz(axx)

	# get inverse using generic solver
	#axx.m.inv <- solve(axx.m)

	# get inverse using toeplitz solver
	axx.m.inv <- ltsa::TrenchInverse(axx.m)

	# solve f * x = z for filter f, given input data x and desired output z

	# multiply the inverse of the autocovariance matrix with the crosscovariance to get the Wiener filter coefficients
	ls.filter <- axx.m.inv %*% azx

	# get the normalized error
	min.error <- getNMSE(azz, azx, ls.filter)

	# get predicted output
	z.est <- my.filter(x, ls.filter)
	if ( demean )
		z.est <- z.est + z.mean

	# directly calculate the normalized error
	if ( demean )
		est.min.error <- getEstNMSE(z + z.mean, z.est, xz.acvf$acf[1,2,2], lag)
	else
		est.min.error <- getEstNMSE(z, z.est, xz.acvf$acf[1,2,2], lag)

	# return a list of the results
	results = list( filter = ls.filter, z.est = z.est,
									min.error = min.error, est.min.error = est.min.error,
									xz.acvf = xz.acvf, axx = axx, azz = azz, azx = azx,
									len = f.len, lag = lag,
									x.mean = x.mean, z.mean = z.mean, xz = xz )
	return ( results )
}
