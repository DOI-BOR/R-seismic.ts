#' Compute the Arias Intensity of a Univariate Time Series
#'
#' @description
#' \code{ai} computes the Arias Intensity and duration of a univariate
#' acceleration time series.
#'
#' @param at Equally-sampled acceleration series (units of cm/s^2 or m/s^2).
#' Must convert to a numeric \code{vector}, \code{ts}, or \code{signalSeries}.
#' @param dt Sample interval, in seconds. Default is 0.01, if input is not a
#' \code{ts}, or \code{signalSeries} object, in which case it is set from the object.
#' @param in.units Units of input vector. One of \code{c("mks", "m/s/s", "m/s^2",
#' "cgs", "cm/s/s", "cm/s^2", "g")} for m/s^2, cm/s^2, or g. Default is "cgs".
#' @return A list with AI, AI.units, dur, dur.units, and ai.integral (which
#' will be the same object type as the input).
#' @examples
#' amp <- 1500 # cm/s^2
#' dt <- .01 # sec
#' f0 <- 0.2 # Hz
#' dur <- 20 # sec
#' syn.ber <- ts(c(rep(0, 0.05 * dur/dt), berlage(amp,f0,dur,dt)),
#'               start=0, deltat=dt)
#' tsplot(syn.ber)
#' ai.ber <- ai(syn.ber,in.units="cgs")
#' sprintf("Arias Intensity = %.5f %s, Arias Duration = %.3f %s",
#'          ai.ber$AI, ai.ber$AI.units, ai.ber$dur, ai.ber$dur.units)
#' tsplot(ai.ber$ai.integral) # ai.integral is a ts object because input was
#' @seealso \href{https://en.wikipedia.org/wiki/Arias_Intensity}{Wikipedia}
#' @keywords ts
#'
ai <- function(at, dt=NA, in.units=NA) {

	if ( is(at, "signalSeries") || is(at, "ts") )
	  dt <- deltat(at)
	if ( is.na(dt) )
	  dt <- 0.01

	start <- 0
	if ( is(at, "signalSeries") ) {
	  in.units <- at@units
	  start <- at@from
	}
	if ( is(at, "ts") )
	  start <- start(at)[1]

	if ( is.na(in.units) )
	  in.units <- "cgs"

	len <- length(at)
	g <- 9.81
	coeff <- 0.5 * pi / g
	if ( in.units == "cgs" || in.units == "cm/s/s" || in.units == "cm/s^2" )
		coeff <- coeff * 1e-4
	else if ( in.units == "g" )
	  coeff <- coeff * g
	ai.integral <- coeff * cumsum(at^2 * dt)

	AI = ai.integral[len]
	ai.integral <- ai.integral / AI
	if ( is(at, "signalSeries") )
	  ai.integral <- signalSeries(ai.integral, from=start, by=dt, units=AI.units)
	else if ( is(at, "ts") )
	  ai.integral <- ts(ai.integral, delta=dt, start=start)

	ai.dur <- length(which(ai.integral >= .05 & ai.integral <= 0.95)) * dt

	ai <- list(AI=AI, AI.units="m/s", dur=ai.dur, dur.units="seconds",
						 ai.integral=ai.integral)

	return(ai)
}


