# least-squares inverse f = x^-1 from directly solving
#		f * x = {0,...,1,0,...}
# finds optimal lag for a given length inverse,
# unless fixed.lag is set to a specified lag
inverseWavelet.ls <- function(wavelet, inv.len = NA, fixed.lag = NA) {

	x <- as.vector(wavelet)
	x.len <- length(x)
	# get inverse filter length
	#	note: the best-performing filter will be infinite length, but as a compromise, use input length (or twice input length)
	if ( is.na(inv.len) )
		inv.len <- 2 * x.len
	f.len <- inv.len
	spike.len = x.len + f.len - 1
	minE <- 1
	minLag <- NA
	minFilter <- NA
	if ( is.na(fixed.lag) ) {
		delta = max(1,length=spike.len / 100)
		for ( lag in seq(from=0, to=(spike.len-1), by=delta) ) {
			if ( lag == 0 )
				spike <- c(1, rep(0, spike.len - 1))
			else if ( lag == spike.len - 1 )
				spike <- c(rep(0, spike.len - 1), 1)
			else
				spike <- c(rep(0, lag), 1, rep(0, spike.len - 1 - lag))
			ls.filt <- ilsd::wienerFilter(x, spike, f.len)
			if ( ls.filt$est.min.error < minE ) {
				minFilter <- ls.filt
				minLag <- lag
				minE <- ls.filt$est.min.error
			}
			print(sprintf("--- %d: %.5f %.5f",lag,ls.filt$min.error,ls.filt$est.min.error))
			# print(spike)
			# print(ls.filt$z.est)
		}
	} else {
		# note: for sufficiently long filters, the optimum lag for minimum delay input is 0,
		#	the optimum lag for maximum delay input is spike.length, and for mixed delay input
		#	the optimum values is in between those limits. Assuming an average mixed delay
		#	input, a good guess for the optimum lag is 0.5 * spike.len
		if ( fixed.lag < 0 ) # magic value indicating function should pick a good value
			fixed.lag = as.integer(0.5 * spike.len)
		if ( fixed.lag == 0 )
			spike <- c(1, rep(0, spike.len - 1))
		else if ( fixed.lag >= spike.len - 1 )
			spike <- c(rep(0, spike.len - 1), 1)
		else
			spike <- c(rep(0, fixed.lag), 1, rep(0, spike.len - 1 - fixed.lag))
		ls.filt <- ilsd::wienerFilter(x, spike, f.len)
		minFilter <- ls.filt
		minLag <- fixed.lag
		minE <- ls.filt$est.min.error
		print(sprintf("--- %d: %.5f %.5f",fixed.lag,ls.filt$min.error,ls.filt$est.min.error))
	}
	results = list( min.lag = minLag, filter = minFilter$filter, z.est = minFilter$z.est,
									min.error = minFilter$min.error, est.min.error = minFilter$est.min.error )
	return (results)
}
