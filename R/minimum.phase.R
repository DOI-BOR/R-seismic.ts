minimum.phase <- function(x.data, zero.pad=TRUE) {
	if ( zero.pad )
		x <- zero.pad(x.data)
	else
		x <- as.vector(x.data)
	x.fft <- fft(x) / length(x)
	x.aspec <- abs(x.fft)
	#x.min.phase <- ilsd::hilbert(log(x.aspec), "hilbert", zero.pad=FALSE)
	x.min.phase <- hilbert(log(x.aspec), "hilbert", zero.pad=FALSE)
	x.min.phase.fft <- x.aspec * exp((0+1i) * x.min.phase)
	# tsplot(x.min.phase)
	x.min.phase.ifft <- fft(x.min.phase.fft, inverse=TRUE)
	x.min.phase <- Re(x.min.phase.ifft)
	return (x.min.phase)
}
