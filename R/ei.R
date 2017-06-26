#' Compute energy-integral intensity and duration
#'
#' \code{ei} computes the intensity and duration of the energy integral
#' of a univariate or multivariate velocity time series.
#'
#' @param vt Equally-sampled velocity series (units of cm/s or m/s)
#' represented as a numeric \code{\link{vector}}, \code{link{matrix}},
#' \code{\link{data.frame}}, \code{\link{ts}}, , \code{\link{mts}}, or
#' \code{\link{signalSeries}}.
#' @param dt Sample interval, in seconds. Normally only used for \code{\link{vector}},
#' \code{\link{matrix}}, or \code{\link{data.frame}} input objects. If input is a
#' \code{\link{ts}}, \code{\link{mts}}, or \code{\link{signalSeries}}, then the sample
#' interval will be taken from the object attributes. However, if \code{\link{deltat}}
#' is not set, then \code{dt} will be used. Default is 0.01 seconds.
#' @param dens Assumed density for the energy integral, in cgs or mks units. Default
#' is 2400 kg/m^3.
#' @param vs Assumed shear-wave velocity for the energy integral, in cgs or mks units.
#' Default 760 m/s.
#' @param in.units units system for the velocity data \code{vt}, \code{dens},
#' and \code{vs}. Each input must use the same units system. Supported units
#' systems are \code{c("mks","cgs"}). If \code{vt} is a \code{\link{signalSeries}},
#' and the \code{units} attribute is set, then the system for those units will be
#' used instead (in which case, it must match that used for \code{dens} and
#' and \code{mu}). Default is \code{"mks"}.
#' @param unscaled If \code{TRUE}, then don't scale the energy integral, which will
#' result in units of m^2/s. If \code{FALSE}, then scale the energy integral by
#' \code{dens * vs}, which will result in units of J/m^2. Default is \code{TRUE}.
#' @return A list with \code{EI, EI.units, dur, dur.units,} and \code{ei.integral}
#' (which will be the same object type as \code{vt}). If \code{vt} is a multivariate
#' time series, then the list elements (e.g., \code{EI}) will be vectors with one
#' element for each time series, and \code{ei.integral} will be a multivariate
#' time series.
#' @seealso \code{\link{ai}}, \code{\link{ts}}, \code{\link{mts}}, \code{\link{signalSeries}}
#' \itemize{
#' \item \href{http://www.iitk.ac.in/nicee/wcee/article/13_243.pdf}{Anderson (2004)}
#' Quantitative measure of the goodness-of-fit of synthetic seismograms.
#' \item \href{http://link.aip.org/link/?EQS/22/985/1}{Kempton and Stewart (2006)}
#'  Equations for Significant Duration of Earthquake Ground Motions Considering Site and Near-Source Effects.
#' \item \href{http://dx.doi.org/10.1016/j.soildyn.2013.04.010}{Trifunac and Todorovska (2013)}
#'  A note on the power of strong ground motion during the January 17, 1994 earthquake in Northridge, California.
#' }
#' @keywords ts

ei <- function(vt, dt=0.01, dens=2400, vs=760, in.units="mks", unscaled=TRUE) {
	if ( is(vt, "signalSeries") || is(vt, "ts") )
	  dt <- deltat(vt)
	if ( ! is.finite(dt) )
	  dt <- 0.01

	start <- 0
	if ( is(vt, "signalSeries") ) {
	  in.units <- vt@units
	  units.position <- vt@units.position
  	if ( is.na(units.position) || is.null(units.position) )
  	  units.position <- "seconds"
	  start <- vt@positions@from
	}
	if ( is.na(in.units) || is.null(in.units) )
	  in.units <- "mks"
	if ( is(vt, "ts") )
	  start <- start(vt)[1]

	if ( unscaled ) {
	  coeff <- 1
	  if ( in.units == "cgs" || in.units == "cm/s" )
	    coeff <- coeff * 1e-4
	  EI.units="m^2/s"
	} else {
  	coeff <- dens * vs
  	if ( in.units == "cgs" || in.units == "cm/s" )
  	  coeff <- coeff * 1e3 * 1e-2 * 1e-4
  	EI.units="J/m^2"
	}

	multi.trace <- is.mts(vt) || is.matrix(vt) || length(dim(vt)) > 1 ||
	  ( is(vt, "signalSeries") && ! is.null(dim(vt)) )
	if ( multi.trace ) {
	  len <- dim(vt)[1]
	  ei.integral <- NULL
	  EI <- NULL
	  EI.dur <- NULL
	  cnames <- NULL
	  for ( ii in 1:dim(vt)[2] ) {
	    if ( is(vt, "signalSeries") ) {
	      xt <- vt@data[,ii]
	      cnames <- colnames(vt@data)
	    } else {
	      xt <- vt[,ii]
	      cnames <- colnames(vt)
	    }
	    ei.int <- coeff * cumsum(xt^2 * dt)
	    ei.val = ei.int[len]
	    ei.int <- ei.int / ei.val
	    ei.dur <- length(which(ei.int >= .05 & ei.int <= 0.95)) * dt
	    if ( is(vt, "signalSeries") ) {
	      if ( is.null(ei.integral) )
	        ei.integral <- signalSeries(ei.int, from=start, by=dt,
	                                    units.position=units.position)
	      else
	        ei.integral <- signalSeries(data.frame(ei.integral@data, ei.int),
	                                    from=start, by=dt, units.position=units.position)
	    } else if ( is(vt, "ts") ) {
	      if ( is.null(ei.integral) )
	        ei.integral <- ts(ei.int, start=start, deltat=dt)
	      else
	        ei.integral <- ts(data.frame(ei.integral, ei.int), start=start, deltat=dt)
	    } else {
	      if ( is.null(ei.integral) )
	        ei.integral <- data.frame(ei.int)
	      else
	        ei.integral <- data.frame(ei.integral, ei.int)
	    }
	    if ( is.null(EI) ) {
	      EI <- c(ei.val)
	      EI.dur <- c(ei.dur)
	    } else {
	      EI <- c(EI, ei.val)
	      EI.dur <- c(EI.dur, ei.dur)
	    }
	  }
	  names(EI) <- cnames
	  names(EI.dur) <- cnames
	  if ( is(vt, "signalSeries") )
	    colnames(ei.integral@data) <- cnames
	  else
	    colnames(ei.integral) <- cnames
	} else {
	  len <- length(vt)
	  if ( is(vt, "signalSeries") ) {
	    xt <- vt@data
	  } else {
	    xt <- vt
	  }
	  ei.int <- coeff * cumsum(xt^2 * dt)
	  EI = ei.int[len]
	  ei.int <- ei.int / EI
	  EI.dur <- length(which(ei.int >= .05 & ei.int <= 0.95)) * dt
	  if ( is(vt, "signalSeries") )
	    ei.integral <- signalSeries(ei.int, from=start, by=dt,
	                                units.position=units.position)
	  else if ( is(vt, "ts") )
	    ei.integral <- ts(ei.int, delta=dt, start=start)
	  else
	    ei.integral <- ei.int
	}

	ret <- list(EI=EI, EI.units=EI.units, dur=EI.dur, dur.units="seconds",
	            ei.integral=ei.integral)

	return(ret)
}