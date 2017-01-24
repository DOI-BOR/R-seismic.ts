signalSeries2ts <- function(ss) {
	ts(ss@data,deltat=ss@positions@by,start=ss@positions@from)
}