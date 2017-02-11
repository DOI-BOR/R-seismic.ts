setAs("signalSeries", "ts",
function(from) {
	ts(from@data, deltat=from@positions@by, start=from@positions@from)
}
)
