#' Compute the Arias Intensity for a time series.
#'
#' \code{ai} is used to compute the Arias Intensity and duration of an input
#' acceleration time series (units of cm/s^2 or m/s^2). An optional Hanning
#' taper may be applied before computing the values.
#'
#' @param at required (actual) equally-sampled acceleration series. Must convert to numeric vector.
#' @param dt (optional) sample interval, in seconds
#' @param in.units (optional). Units of input vector."mks" or "cgs"
#' @return a list with the AI, duration, and units
ai <- function(at, dt=0.01, in.units="cgs") {

	len <- length(at)

	g <- 9.81
	coeff <- 0.5 * pi / g
	if ( in.units == "cgs" )
		coeff <- coeff * 1e-4
	ai.integral <- coeff * cumsum(at^2 * dt)

	AI = ai.integral[len]
	ai.integral <- ai.integral / AI

	ai.dur <- length(which(ai.integral >= .05 & ai.integral <= 0.95)) * dt

	ai <- list(AI=AI, AI.units="m/s", dur=ai.dur, dur.units="seconds",
						 ai.integral=ai.integral)

	return(ai)
}


