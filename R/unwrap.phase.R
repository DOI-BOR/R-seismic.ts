unwrap.phase <- function(x.data) {
	x <- as.vector(x.data)
	max.jump <- pi
	x.last <- x[1]
	phase <- 0
	xu <- NULL
	for ( x.cur in x ) {
		x.diff <- x.cur - x.last
		x.last <- x.cur
		if ( x.diff > max.jump ) {
			phase <- phase - 2 * pi
		} else if ( x.diff < -max.jump ) {
			phase <- phase + 2 * pi
		}
		xu <- c(xu, x.cur + phase)
	}
	return(xu)
}
