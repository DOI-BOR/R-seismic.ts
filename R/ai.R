#' Compute Arias intensity and duration
#'
#' @description
#' \code{ai} computes the Arias Intensity and duration of a univariate
#' or multivariate acceleration time series.
#'
#' @param at Equally-sampled acceleration series (units of cm/s^2, m/s^2, or g)
#' represented as a numeric \code{\link{vector}}, \code{link{matrix}},
#' \code{\link{data.frame}}, \code{\link{ts}}, \code{\link{mts}}, or
#' \code{\link{signalSeries}} object.
#' @param dt Sample interval, in seconds. Normally only used for \code{\link{vector}},
#' \code{\link{matrix}}, or \code{\link{data.frame}} input objects. If input is a
#' \code{\link{ts}}, \code{\link{mts}}, or \code{\link{signalSeries}}, then the sample
#' interval will be taken from the object attributes. However, if \code{\link{deltat}}
#' is not set, then \code{dt} will be used. Default is 0.01 seconds.
#' @param in.units Units of the input acceleration data. One of \code{c("mks", "m/s/s",
#' "m/s^2", "cgs", "cm/s/s", "cm/s^2", "g")} for m/s^2, cm/s^2, or g. If the input
#' is a \code{\link{signalSeries}} object, and the \code{units} attribute is set,
#' then that will be used instead. Default is \code{"mks"}.
#' @return A list with \code{AI, AI.units, dur, dur.units,} and \code{ai.integral}
#' (which will be the same object type as \code{at}). If \code{at} is a multivariate
#' time series, then the list elements (e.g., \code{AI}) will be vectors with one
#' element for each time series, and \code{ai.integral} will be a multivariate
#' time series.
#' @examples
#' amp <- 1500 # cm/s^2
#' dt <- .01 # sec
#' f0 <- 0.2 # Hz
#' dur <- 20 # sec
#' syn.ber <- ts(c(rep(0, 0.05 * dur/dt), berlage(amp,f0,dur,dt)),
#'               start=0, deltat=dt)
#' tsplot(syn.ber)
#' ai.ber <- ai(syn.ber, in.units="cgs")
#' sprintf("Arias Intensity = %.5f %s, Arias Duration = %.3f %s",
#'          ai.ber$AI, ai.ber$AI.units, ai.ber$dur, ai.ber$dur.units)
#' tsplot(ai.ber$ai.integral) # ai.integral is a ts object because input was
#' @seealso \code{\link{ei}}, \code{\link{ts}}, \code{\link{mts}}, \code{\link{signalSeries}}
#' \itemize{
#' \item \href{https://en.wikipedia.org/wiki/Arias_Intensity}{Wikipedia}
#' \item \href{http://www.iitk.ac.in/nicee/wcee/article/13_243.pdf}{Anderson (2004)}
#' Quantitative measure of the goodness-of-fit of synthetic seismograms.
#' \item \href{http://link.aip.org/link/?EQS/22/985/1}{Kempton and Stewart (2006)}
#'  Equations for Significant Duration of Earthquake Ground Motions Considering Site and Near-Source Effects.
#' }
#' @keywords ts
#'
ai <- function(at, dt=NA, in.units=NA) {

  if ( is(at, "signalSeries") || is(at, "ts") )
	  dt <- deltat(at)
  if ( ! is.finite(dt) )
    dt <- 0.01

	start <- 0
	if ( is(at, "signalSeries") ) {
	  in.units <- at@units
	  units.position <- at@units.position
	  if ( is.na(units.position) || is.null(units.position) )
	    units.position <- "seconds"
	  start <- at@positions@from
	}
	if ( is.na(in.units) || is.null(in.units) )
	  in.units <- "mks"
	if ( is(at, "ts") )
	  start <- start(at)[1]

	g <- 9.81
	coeff <- 0.5 * pi / g
	if ( in.units == "cgs" || in.units == "cm/s/s" || in.units == "cm/s^2" )
	  coeff <- coeff * 1e-4
	else if ( in.units == "g" )
	  coeff <- coeff * g

	multi.trace <- is.mts(at) || is.matrix(at) || length(dim(at)) > 1 ||
    ( is(at, "signalSeries") && ! is.null(dim(at)) )
  if ( multi.trace ) {
    len <- dim(at)[1]
    ai.integral <- NULL
    AI <- NULL
    AI.dur <- NULL
    cnames <- NULL
    for ( ii in 1:dim(at)[2] ) {
      if ( is(at, "signalSeries") ) {
        xt <- at@data[,ii]
        cnames <- colnames(at@data)
      } else {
        xt <- at[,ii]
        cnames <- colnames(at)
      }
      ai.int <- coeff * cumsum(xt^2 * dt)
      ai.val = ai.int[len]
      ai.int <- ai.int / ai.val
      ai.dur <- length(which(ai.int >= .05 & ai.int <= 0.95)) * dt
      if ( is(at, "signalSeries") ) {
        if ( is.null(ai.integral) )
          ai.integral <- signalSeries(ai.int, from=start, by=dt,
                                   units.position=units.position)
        else
          ai.integral <- signalSeries(data.frame(ai.integral@data, ai.int),
                                   from=start, by=dt, units.position=units.position)
      } else if ( is(at, "ts") ) {
        if ( is.null(ai.integral) )
          ai.integral <- ts(ai.int, start=start, deltat=dt)
        else
          ai.integral <- ts(data.frame(ai.integral, ai.int), start=start, deltat=dt)
      } else {
        if ( is.null(ai.integral) )
          ai.integral <- data.frame(ai.int)
        else
          ai.integral <- data.frame(ai.integral, ai.int)
      }
      if ( is.null(AI) ) {
        AI <- c(ai.val)
        AI.dur <- c(ai.dur)
      } else {
        AI <- c(AI, ai.val)
        AI.dur <- c(AI.dur, ai.dur)
      }
    }
    names(AI) <- cnames
    names(AI.dur) <- cnames
    if ( is(at, "signalSeries") )
      colnames(ai.integral@data) <- cnames
    else
      colnames(ai.integral) <- cnames

  } else {
	  len <- length(at)
	  if ( is(at, "signalSeries") ) {
	    xt <- at@data
	  } else {
	    xt <- at
	  }
	  ai.int <- coeff * cumsum(xt^2 * dt)
	  AI = ai.int[len]
	  ai.int <- ai.int / AI
	  AI.dur <- length(which(ai.int >= .05 & ai.int <= 0.95)) * dt
  	if ( is(at, "signalSeries") )
  	  ai.integral <- signalSeries(ai.int, from=start, by=dt,
  	                              units.position=units.position)
  	else if ( is(at, "ts") )
  	  ai.integral <- ts(ai.int, delta=dt, start=start)
	  else
	    ai.integral <- ai.int
  }

	ret <- list(AI=AI, AI.units="m/s", dur=AI.dur, dur.units="seconds",
						 ai.integral=ai.integral)

	return(ret)
}


