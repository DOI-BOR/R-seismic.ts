#' Finite Difference Derivative of a Multivariate Time Series
#'
#' @description
#' \code{deriv_fd} is used to compute the derivative of a univariate or multivariate
#' time series using the finite difference method.
#'
#' @param xt Equally-sampled input series. Must convert to numeric vector.
#' @param dt Sample interval, in seconds. Default is 0.01 if input is not a
#' \code{\link{ts}} or \code{\link{signalSeries}}.
#' @param nd Integer order of the derivative. Currently only first and
#' second derivatives are supported. Default is 1.
#' @param order Integer order of the finite difference. Currently only
#' orders of 2, 4, 6, and 8 are supported. Default is 8.
#' @param pct Percentage of data window to apply a Hanning taper. Must be
#' between 0 and 50. Default is 0.
#' @return the derivative of the windowed data
#' @keywords ts
#'

#' @describeIn deriv_fd.default differentiates a numeric \code{vector} or \code{matrix}.
deriv_fd.default <- function(xt, dt=NA, nd=1, order=8, pct=NA) {
  order <- if ( order == 8 ) 3 else if ( order == 6 ) 2 else if ( order == 4 ) 1 else 0
  if ( is.na(dt) )
    dt <- 0.01

  multi.trace <- is.matrix(xt) || length(dim(xt)) > 1
  if ( multi.trace ) {
    xt.len <- dim(xt)[1]
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    dxdt.mts = NULL
    for ( ii in 1:dim(xt)[2] ) {
      ok <- ! is.na(xt[,ii])
      dxdt <- .Call("CALLfd_deriv",
                   as.double(xt[ok,ii]), as.double(dt), as.integer(nd),
                   as.integer(order), as.double(pct))
      if ( is.null(dxdt.mts) )
        dxdt.mts <- data.frame(dxdt)
      else
        dxdt.mts <- data.frame(dxdt.mts,dxdt)
    }
    dxdt <- dxdt.mts
  } else {
    xt <- as.double(xt[!is.na(xt)])
    xt.len <- length(xt)
  	if ( xt.len < 3 )
  		stop("input time series must have at least 3 valid points")
    dxdt <- .Call("CALLfd_deriv",
  							 xt, as.double(dt), as.integer(nd),
  							 as.integer(order), as.double(pct))
 }

	return(dxdt)
}
setGeneric("deriv_fd",def=deriv_fd.default)

#' @describeIn deriv_fd.default differentiates a \code{ts} or \code{mts} object.
deriv_fd.ts <- function(xt, dt=NA, nd=1, order=8, pct=NA) {
  order <- if ( order == 8 ) 3 else if ( order == 6 ) 2 else if ( order == 4 ) 1 else 0

  multi.trace <- is.mts(xt)

  dt <- deltat(x.data)
  if ( is.na(dt) )
    dt <- 0.01
  start <- start(xt)[1]

  if ( multi.trace == TRUE ) {
    xt.len <- dim(xt)[1]
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    dxdt.mts = NULL
    for ( ii in 1:dim(xt)[2] ) {
      ok <- ! is.na(xt[,ii])
      dxdt <- .Call("CALLfd_deriv",
                    as.double(xt[ok,ii]), as.double(dt), as.integer(nd),
                    as.integer(order), as.double(pct))
      if ( is.null(dxdt.mts) )
        dxdt.mts <- ts(dxdt, start = start, deltat = dt)
      else
        dxdt.mts <- ts(data.frame(dxdt.mts, dxdt), start = start, deltat = dt)
    }
    dxdt <- dxdt.mts
  } else {
    xt <- as.double(xt[!is.na(xt)])
    xt.len <- length(xt)
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    dxdt <- .Call("CALLfd_deriv",
                  xt, as.double(dt), as.integer(nd),
                  as.integer(order), as.double(pct))
    dxdt <- ts(dxdt, start = start, deltat = dt)
  }
  dimnames(dxdt) <- dimnames(xt)

  return(dxdt)
}
setMethod("deriv_fd","ts",deriv_fd.ts)

#' @describeIn deriv_fd.default differentiates a \code{signalSeries} object.
deriv_fd.signalSeries <- function(xt, dt=NA, nd=1, order=8, pct=NA) {
  order <- if ( order == 8 ) 3 else if ( order == 6 ) 2 else if ( order == 4 ) 1 else 0

  multi.trace <- ! is.null(dim(xt))

  dt <- deltat(x.data)
  if ( is.na(dt) )
    dt <- 0.01
  start <- xt@positions@from

  units <- xt@units
  if ( grep("cm/s/s",units) || grep("cm/s^2",units) )
    new.units <- "cm/s^3"
  else if ( grep("m/s/s",units) || grep("m/s^2",units) )
    new.units <- "m/s^3"
  else if ( grep("cm/s",units) )
    new.units <- "cm/s^2"
  else if ( grep("m/s",units) )
    new.units <- "m/s^2"
  else if ( grep("cm",units) )
    new.units <- "cm/s"
  else if ( grep("m",units) )
    new.units <- "m/s"
  else if ( grep("ft/s/s",units) || grep("ft/s^2",units) )
    new.units <- "ft/s^2"
  else if ( grep("ft/s",units) )
    new.units <- "ft/s^2"
  else if ( grep("ft",units) )
    new.units <- "ft/s"
  else
    new.units <- NULL

  if ( multi.trace == TRUE ) {
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    dxdt.mts = NULL
    for ( ii in 1:dim(xt)[2] ) {
      ok <- ! is.na(xt[,ii]@data)
      dxdt <- .Call("CALLfd_deriv",
                    as.double(xt[ok,ii]@data), as.double(dt), as.integer(nd),
                    as.integer(order), as.double(pct))
      if ( is.null(xt) )
        dxdt.mts <- signalSeries(dxdt, from = start, by = dt, units = new.units)
      else
        dxdt.mts <- signalSeries(data.frame(dxdt.mts@data, dxdt@data), from = start, by = dt, units = new.units)
    }
    dxdt <- dxdt.mts
  } else {
    ok <- ! is.na(xt@data)
    xt.len <- length(xt[ok])
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    dxdt <- .Call("CALLfd_deriv",
                  as.double(xt[ok]@data), as.double(dt), as.integer(nd),
                  as.integer(order), as.double(pct))
    dxdt <- signalSeries(dxdt, from = start, by = dt, units = new.units)
  }
  names(dxdt) <- names(xt)

  return(dxdt)
}
setMethod("deriv_fd","signalSeries",deriv_fd.signalSeries)


