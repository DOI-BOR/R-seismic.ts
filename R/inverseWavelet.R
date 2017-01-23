# compute wavlet inverse by directly solving at * bt = 1
inverseWavelet <- function(wavelet, inv.len = NA) {
	bt <- as.vector(wavelet)
	bt.len <- length(bt)
	if ( is.na(inv.len) )
		inv.len <- bt.len
	at.len <- inv.len
	at <- 1. / bt[1]
	for ( tt in seq(from=2, by=1, length=(at.len-1)) ) {
		if ( tt > bt.len ) {
			btmax <- bt.len
		} else {
			btmax <- tt
		}
		at <- c(at, -sum(at[(tt-1):(tt-btmax+1)]*bt[2:btmax])/bt[1])
	}
	at
}
