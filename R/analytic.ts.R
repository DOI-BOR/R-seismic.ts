# returns f(t) - i * hilbert(f(t))
#	uses property that f(t) - i * hilbert(f(t)) <=> 2 * H(omega) * F(omega)
#	where H() is the Heaviside function, and f(t) <=> F(omega)

analytic.ts <- function(x.data) {
	# check for multi-trace input
	multi.trace <- is.mts(x.data) || is.matrix(x.data) || length(dim(x.data)) > 1 ||
		( is(x.data, "signalSeries") && ! is.null(dim(x.data)) )
	if ( multi.trace == TRUE ) {
		x.len <- dim(x.data)[1]
	} else {
		x.len <- length(x.data)
	}

	# make the Heaviside function times 2. For fft(), note that at discontinuities, h(f) -> 0.5 * [h(f-) + h(f+)],
	#	and that the indices are defined as follows:
	#	n even
	#		1 = zero frequency (discontinuity)
	#		2:n/2 = positive frequencies
	#		n/2+1 = nyquist (discontinuity)
	#		n/2+2:n = negative frequencies
	#	n odd
	#		1 = zero frequency (discontinuity)
	#		2:(n+1)/2 = positive frequencies
	#		(n+1)/2+1:n = negative frequencies
	n <- x.len
	Hf <- rep(0,n) # initialize all values to zero
	if ( n == 2 * floor(n/2) ) {
		# n is even
		Hf[c(1, n/2 + 1)] = 1	# zero and Nyquist frequencies - disconinuity value = 0.5*(0+2) = 1
		Hf[2:n/2] = 2	# positive frequencies
	} else {
		# n is odd
		Hf[1] = 1 # zero frequency - discontinuity value = 1
		Hf[2:(n + 1)/2] = 2	# positive frequencies
	}

	if ( multi.trace == TRUE ) {
		xa <- x.data
		for ( cn in 1:dim(x.data)[2] ) {
			x <- as.vector(x.data[,cn])
			Xf <- fft(x) / n   # if length(x) is prime, the fft is just the dft
			# multiply dft of input by Heaviside function, and take inverse transform
			xa[,cn] <- fft(Xf * Hf, inverse=TRUE)
		}
	} else {
		# get dft of the input data
		x <- as.vector(x.data)
		Xf <- fft(x) / n   # if length(x) is prime, the fft is just the dft
		# multiply dft of input by Heaviside function, and take inverse transform
		xa <- fft(Xf * Hf, inverse=TRUE)
	}

	return (xa)
}
