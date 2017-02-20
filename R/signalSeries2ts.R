setAs("signalSeries", "ts",
function(from) {
  warning("ss2ts: start=",from@positions@from," by=",from@positions@by)
	ts(from@data, start=from@positions@from, deltat=from@positions@by)
}
)

as.ts.signalSeries <- function(x, ...) {
  ts(x@data, start=x@positions@from, deltat=x@positions@by, ...)
}
setMethod("as.ts","signalSeries",as.ts.signalSeries)
