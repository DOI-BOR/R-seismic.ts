# get prediction filter
predictionErrorFilter <- function(x.data, length.filt = NULL, lag = NULL, length.wavelet = NA) {
	pred.filt <- ilsd::predictionFilter(x.data, length.filt, lag)

	# get prediction error filter
	if ( pred.filt$lag > 1 )
		perr.filt <- c(1, rep(0, lag - 1), -pred.filt$filter)
	else
		perr.filt <- c(1, -pred.filt$filter)

	# get the inovations
	ino <- ilsd::my.filter(pred.filt$x, perr.filt)

	# if lag = 1, get the wavelet, which is the inverse of the (assumed minimum delay) prediction error filter
	wavelet <- NA
	if ( pred.filt$lag == 1 )
		wavelet <- ilsd::inverseWavelet(perr.filt, length.wavelet)

	results = list( filter = perr.filt, lag = pred.filt$lag,
									min.error = pred.filt$min.error, est.min.error = pred.filt$est.min.error,
									ino = ino, wavelet = wavelet, x.mean = pred.filt$x.mean )
	return ( results )
}
