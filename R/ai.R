#' Compute the Arias Intensity for a time series
#'
#' @description
#' \code{ai} computes the Arias Intensity and duration of a univariate
#' acceleration time series.
#'
#' @param at Equally-sampled acceleration series (units of cm/s^2 or m/s^2).
#' Must convert to a numeric vector.
#' @param dt Sample interval, in seconds.
#' @param in.units Units of input vector. One of c("mks", "cgs", "g") for
#' m/s^2, cm/s^2, or g
#' @return A list with AI, AI.units, dur, dur.units, and ai.integral
#' @examples
#' amp <- 1500 # cm/s^2
#' dt <- .01 # sec
#' f0 <- 0.2 # Hz
#' dur <- 20 # sec
#' syn.ber <- ts(c(rep(0, 0.05 * dur/dt), berlage(amp,f0,dur,dt)),
#'               start=0, deltat=dt)
#' tsplot(syn.ber)
#' ai.ber <- ai(syn.ber,dt,in.units="cgs")
#' sprintf("Arias Intensity = %.5f %s, Arias Duration = %.3f %s",
#'          ai.ber$AI, ai.ber$AI.units, ai.ber$dur, ai.ber$dur.units)
#' tsplot(ts(ai.ber$ai.integral,start=0,deltat=dt))
#' @seealso \href{https://en.wikipedia.org/wiki/Arias_Intensity}{Wikipedia}
#' @keywords ts
#'
ai <- function(at, dt=0.01, in.units="cgs") {

	len <- length(at)

	g <- 9.81
	coeff <- 0.5 * pi / g
	if ( in.units == "cgs" )
		coeff <- coeff * 1e-4
	else if ( in.units == "g" )
	  coeff <- coeff * g
	ai.integral <- coeff * cumsum(at^2 * dt)

	AI = ai.integral[len]
	ai.integral <- ai.integral / AI

	ai.dur <- length(which(ai.integral >= .05 & ai.integral <= 0.95)) * dt

	ai <- list(AI=AI, AI.units="m/s", dur=ai.dur, dur.units="seconds",
						 ai.integral=ai.integral)

	return(ai)
}


