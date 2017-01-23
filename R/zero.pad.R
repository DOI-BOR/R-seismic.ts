# zero pad input series to padded.len (if defined), or to the
# next power of two
# Required Arguments:
#		x.data:
#			vector data
# Optional Arguments
#		padded.len:
#			if padded.len = NA, then pad the current lenght to the
#			next power of two.
#		extra: If extra = TRUE, and padded.len=NA, then pad
#			to the next next power of two
zero.pad <- function(x.data, padded.len = NA, extra = FALSE) {
	x <- as.vector(x.data)
	n = length(x)
	if ( is.na(padded.len) ) {
		n2 = 2
		while ( n2 < n )
			n2 <- 2 * n2
		if ( extra )
			n2 <- 2 * n2
		x.aug = c(x, rep(0, n2 - n))
	} else {
		npad = as.integer(round(padded.len)) - n
		if ( npad > 0 )
			x.aug = c(x, rep(0, npad))
		else
			x.aug = x
	}
	return (x.aug)
}
