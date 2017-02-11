#' Compute the energy integral for a time series.
#'
#' \code{ei} is used to compute the energy integral of an input velocity
#' time series (units of cm/s or m/s). An optional Hanning
#' taper may be applied before computing the values.
#'
#' @param at required (actual) equally-sampled acceleration series. Must convert to numeric vector.
#' @param dt (optional) sample interval, in seconds
#' @param in.units (optional). Units of input vector, density, and mu. "mks" or "cgs"
#' @return a list with the the total energy and units
#' @keywords ts
ei <- function(vt, dt=0.01, dens=2.7, mu=3e11, in.units="cgs") {

	len <- length(vt)

	beta <- sqrt(mu / dens)
	coeff <- dens * beta
	if ( in.units == "cgs" )
		coeff <- coeff * 1e3 * 1e-2 * 1e-4
	ei.integral <- coeff * cumsum(vt^2 * dt)

	EI = ei.integral[len]
	ei.integral <- ei.integral / EI

	ei.dur <- length(which(ei.integral >= .05 & ei.integral <= 0.95)) * dt

	ei <- list(EI=EI, EI.units="J/m^2", dur=ei.dur, dur.units="seconds",
						 ei.integral=ei.integral)

	return(ei)
}