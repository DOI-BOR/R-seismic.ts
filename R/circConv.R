# direct circular convolution
circConv <- function(x, filt) {
	xt <- as.vector(x)
	nx <- length(xt)
	ft <- as.vector(filt)
	nf <-length(ft)
	cc <- ft[1] * xt[1]
	for ( tt in seq(from=2, by=1, length=(nx + nf - 2)) ) {
		if ( tt > nf ) {
			ftmax <- nf
			xtmin <- tt - nf + 1
		} else {
			ftmax <- tt
			xtmin <- 1
		}
		if ( tt > nx ) {
			ftmin <- tt - nx + 1
			xtmax <- nx
		} else {
			ftmin <- 1
			xtmax <- tt
		}
		cc <- c(cc, sum(ft[ftmin:ftmax]*xt[xtmax:xtmin]))
	}
	cc
}
