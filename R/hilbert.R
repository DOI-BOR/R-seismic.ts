hilbert <- function(x.data, op = "hilbert", zero.pad = TRUE, dt = 0.01) {
	# note: length(x) does not have to be a power of 2
	if ( zero.pad )
		x <- zero.pad(x.data)
	else
		x <- as.vector(x.data)
	#xa <- ilsd::analytic.ts(x)
	xa <- analytic.ts(x)
	result <- NULL
	if ( op == "hilbert" ) {
		result <- -Im(xa)
	} else if ( op == "envelope" ) {
		result <- abs(xa)
	} else if ( op == "instphase" ) {
		inst.phase <- atan2(Im(xa),Re(xa))
		# result <- ilsd::unwrap.phase(inst.phase)
		result <- unwrap.phase(inst.phase)
	} else if ( op == "instfreq" ) {
		inst.phase <- atan2(Im(xa),Re(xa))
		# result <- diff(ilsd::unwrap.phase(inst.phase),1) / (2 * pi * dt)
		result <- diff(unwrap.phase(inst.phase),1) / (2 * pi * dt)
	}
	return (result)
}
