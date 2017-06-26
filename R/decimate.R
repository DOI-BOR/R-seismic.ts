#' Down-sample a numeric time series
#'
#' @description
#' \code{decimate} down-samples a numeric univariate or multivariate time series
#' after (optionally) low-pass filtering.
#'
#' @param xt Equally-sampled input time series, including a numeric \code{\link{vector}},
#' \code{\link{matrix}}, \code{\link{data.frame}}, \code{\link{ts}}, \code{\link{mts}}
#' or \code{\link{signalSeries}}.
#' @param factor Decimation factor (default is 1). If \code{factor >= 1}, then
#' \code{factor} will be rounded to the nearest integer, and the time-series
#' will be decimated. If \code{factor < 1}, then the time series will be
#' augmented by \code{(1/factor)} using \code{\link{augment}}.
#' @param lp.filter Set to \code{TRUE} to suitably low-pass filter the time
#' series before down-sampling, to avoid aliasing. Default is \code{FALSE},
#' which assumes the user has already done the necessary low-pass filtering.
#' @return The down-sampled time series is returned, with the same object type
#' as the input. The sampling interval \code{deltat} will be updated
#' for \code{\link{ts}} and \code{\link{signalSeries}} objects.
#' @details If \code{factor = N > 1}, then a new time series is constructed by
#' taking every \code{N}th point of the original time series. To avoid aliasing
#' of the down-sampled result, it is critical that the input be suitably
#' low-pass filtered. To have the greatest control, the user can do this
#' separately, before calling this function. Otherwise, set \code{lp.filter = TRUE},
#' in which case a 4-pole, zero-phase, low-pass Bessel filter at 2/3 of the new
#' Nyquist frequency will be applied to the original time series, using
#' \code{\link{filter.llnl}}.
#' @seealso \code{\link{filter.llnl}}, \code{\link{augment}}, \code{\link{ts}},
#' \code{\link{mts}}, and \code{\link{signalSeries}}.
#' @keywords ts

#' @describeIn decimate.default decimate a numeric \code{\link{vector}}, or the
#' columns of a numeric \code{\link{matrix}} or \code{\link{data.frame}}.
decimate.default <- function(xt, factor = NA, lp.filter=FALSE) {
  if ( is.na(factor) || (factor >= 1 && as.integer(round(factor)) == 1) )
    return(xt)

	multi.trace <- is.matrix(xt) || length(dim(xt)) > 1
	if ( multi.trace ) {
	  xt.len <- dim(xt)[1]
	  if ( xt.len < 2 )
	    stop("input time series must have at least 2 valid points")
	  xdt.mts = NULL
	  for ( ii in 1:dim(xt)[2] ) {
	    xdt <- dec.vec(xt[,ii], factor, lp.filter)
	    if ( is.null(xdt.mts) )
	      xdt.mts <- data.frame(xdt)
	    else
	      xdt.mts <- data.frame(xdt.mts, xdt)
	  }
	  colnames(xdt.mts) <- colnames(xt)
	  xt.dec <- xdt.mts
	} else {
	  len <- length(xt)
	  if ( len < 2 )
	    stop("input time series must have at least 2 valid points")
	  xt.dec <- dec.vec(xt, factor, lp.filter)
	}

	return(xt.dec)
}
setGeneric("decimate", def=decimate.default)

#' @describeIn decimate.default decimate a numeric \code{\link{ts}} or
#' \code{\link{mts}} time series.
decimate.ts <- function(xt, factor = NA, lp.filter=FALSE) {
  if ( is.na(factor) || (factor >= 1 && as.integer(round(factor)) == 1) )
    return(xt)

  dt <- deltat(xt)
  if ( ! is.finite(dt) || dt == 0 )
    dt <- 1
  if ( factor < 1 )
    dt <- dt * factor
  else
    dt <- dt * as.integer(round(factor))
  start <- start(xt)[1]

  multi.trace <- is.mts(xt)
  if ( multi.trace ) {
    xt.len <- dim(xt)[1]
    if ( xt.len < 2 )
      stop("input time series must have at least 2 valid points")
    xdt.mts = NULL
    for ( ii in 1:dim(xt)[2] ) {
      xdt <- dec.vec(xt[,ii], factor, lp.filter)
      if ( is.null(xdt.mts) )
        xdt.mts <- data.frame(xdt)
      else
        xdt.mts <- data.frame(xdt.mts, xdt)
    }
    colnames(xdt.mts) <- colnames(xt)
    xt.dec <- xdt.mts
  } else {
    len <- length(xt)
    if ( len < 2 )
      stop("input time series must have at least 2 valid points")
    xt.dec <- dec.vec(xt, factor, lp.filter)
  }
  xt.dec <- ts(xt.dec, start=start, deltat=dt)

  return(xt.dec)
}
setMethod("decimate", "ts", decimate.ts)

#' @describeIn decimate.default decimate a numeric \code{\link{signalSeries}}
#' time series.
decimate.signalSeries <- function(xt, factor = NA, lp.filter=FALSE) {
  if ( is.na(factor) || (factor >= 1 && as.integer(round(factor)) == 1) )
    return(xt)

  dt <- deltat(xt)
  if ( ! is.finite(dt) || dt == 0 )
    dt <- 1
  if ( factor < 1 )
    dt <- dt * factor
  else
    dt <- dt * as.integer(round(factor))
  start <- xt@positions@from
  units <- xt@units
  units.position <- xt@units.position

  multi.trace <- ! is.null(dim(xt))
  if ( multi.trace ) {
    xt.len <- dim(xt)[1]
    if ( xt.len < 2 )
      stop("input time series must have at least 2 valid points")
    xt.dec = NULL
    for ( ii in 1:dim(xt)[2] ) {
      xdt <- dec.vec(xt@data[,ii], factor, lp.filter)
      if ( is.null(xdt.mts) )
        xt.dec <- data.frame(xdt)
      else
        xt.dec <- data.frame(xt.dec, xdt)
    }
    colnames(xt.dec) <- colnames(xt@data)
  } else {
    len <- length(xt)
    if ( len < 2 )
      stop("input time series must have at least 2 valid points")
    xt.dec <- dec.vec(xt@data, factor, lp.filter)
  }
  xt.dec <- signalSeries(xt.dec, from=start, by=dt, units=units,
                          units.position=units.position)

  return(xt.dec)
}
setMethod("decimate", "signalSeries", decimate.signalSeries)

# decimate a numeric vector
dec.vec <- function(xt, factor=NA, lp.filter=FALSE) {
  if ( ! is.numeric(xt) )
    error("decimate: only numeric time series are supported")

  if ( factor < 1. ) {
    x.dec <- augment(xt, 1./factor)
  } else {
    xt <- as.vector(xt)
    skip <- as.integer(round(factor))
    if ( lp.filter ) {
      f.hi <- 0.67 * 0.5 / skip # Set high-cut frequency to 2/3 of the new Nyquist
      xt <- filter.llnl(xt, dt=1, order=4, pb.type="lp", filt.type="bessel",
                    f.hi=f.hi, dir="zp")
    }
    x.len <- length(xt)
    x.dec <- xt[seq(1,x.len,skip)]
  }
  return ( x.dec )
}
