# Decimate an augmented time series by a factor not greater
# than 1/aug_factor
#   note: low-pass filtering is assumed to not be needed
#		because the data was augmented, or has been suitably
#		filtered. This method should not be used to decimate
# 	a time series that has non-zero spectral amplitudes
#		for frequencies above 1 / (dt * factor)
# Required Arguments:
#		x.data:
#			series, signalSeries, or vector data
# Optional Arguments
#		factor:
#			Decimation factor. If factor = N, then dt -> dt * N
decimate <- function(x.data, factor = NA) {
	if ( is.na(factor) )
		factor = 1.;
	skip <- as.integer(round(factor))
	if ( factor < 1. ) {
		x.dec <- ilsd::augment(x.data, 1./factor)
	} else if ( skip > 1 ) {
		if ( is(x.data, "series") || is(x.data, "signalSeries") )
			x <- x.data
		else
			x <- splus2R::signalSeries(as.vector(x.data))
		x.len <- length(x)
		x.dec <- splus2R::signalSeries(x@data[seq(1,x.len,skip)],from=x@positions@from,by=(x@positions@by*skip))
	} else {
		x.dec <- x.data
	}
	return ( x.dec )
}
