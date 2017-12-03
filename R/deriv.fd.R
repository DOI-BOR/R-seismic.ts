#' Finite Difference Derivative of a Multivariate Time Series
#'
#' \code{deriv_fd} is used to compute the derivative of a univariate or multivariate
#' time series using the finite difference method.
#'
#' @param xt Equally-sampled input time series, including a numeric \code{\link{vector}},
#' \code{\link{matrix}}, \code{\link{data.frame}}, \code{\link{ts}}, \code{\link{mts}}
#' or \code{\link{signalSeries}}.
#' @param dt Sample interval, in seconds. Default is 0.01 if input is not a
#' \code{\link{ts}} or \code{\link{signalSeries}}.
#' @param nd Integer order of the derivative. Currently only first and
#' second derivatives are supported. Default is 1.
#' @param order Integer order of the finite difference. Currently only
#' orders of 2, 4, 6, and 8 are supported. Default is 8.
#' @param pct Percentage of data window to apply a Hanning taper. Must be
#' between 0 and 50. Default is 0.5 percent.
#' @return the derivative of the windowed data, with the same object type as
#' the input.
#' @seealso \code{\link{ts}}, \code{\link{mts}}, \code{\link{signalSeries}},
#' and \code{\link{integ}}.
#' @keywords ts
#'

#' @describeIn deriv_fd.default differentiates a numeric \code{vector} or \code{matrix}.
deriv_fd.default <- function(xt, dt=NA, nd=1, order=8, pct=0.5) {
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
                   as.integer(order))
      if ( is.null(dxdt.mts) )
        dxdt.mts <- data.frame(dxdt)
      else
        dxdt.mts <- data.frame(dxdt.mts,dxdt)
    }
    colnames(dxdt.mts) <- colnames(xt)
    dxdt <- dxdt.mts
  } else {
    xt <- as.double(xt[!is.na(xt)])
    xt.len <- length(xt)
  	if ( xt.len < 3 )
  		stop("input time series must have at least 3 valid points")
    dxdt <- .Call("CALLfd_deriv",
  							 xt, as.double(dt), as.integer(nd),
  							 as.integer(order))
  }
  if ( is.finite(pct) && pct > 0 )
    dxdt <- hanning(dxdt,pct)

  return(dxdt)
}
setGeneric("deriv_fd",def=deriv_fd.default)

#' @describeIn deriv_fd.default differentiates a \code{ts} or \code{mts} object.
deriv_fd.ts <- function(xt, dt=NA, nd=1, order=8, pct=0.5) {
  order <- if ( order == 8 ) 3 else if ( order == 6 ) 2 else if ( order == 4 ) 1 else 0

  multi.trace <- is.mts(xt)
  dt <- deltat(xt)
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
                    as.integer(order))
      if ( is.null(dxdt.mts) )
        dxdt.mts <- data.frame(dxdt)
      else
        dxdt.mts <- data.frame(dxdt.mts, dxdt)
    }
    colnames(dxdt.mts) <- colnames(xt)
    dxdt <- dxdt.mts
  } else {
    xt <- as.double(xt[!is.na(xt)])
    xt.len <- length(xt)
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    dxdt <- .Call("CALLfd_deriv",
                  xt, as.double(dt), as.integer(nd),
                  as.integer(order))
  }
  if ( is.finite(pct) && pct > 0 )
    dxdt <- hanning(dxdt,pct)
  dxdt <- ts(dxdt, start=start, deltat=dt)

  return(dxdt)
}
setMethod("deriv_fd","ts",deriv_fd.ts)

#' @describeIn deriv_fd.default differentiates a \code{signalSeries} object.
deriv_fd.signalSeries <- function(xt, dt=NA, nd=1, order=8, pct=0.5) {
  order <- if ( order == 8 ) 3 else if ( order == 6 ) 2 else if ( order == 4 ) 1 else 0

  multi.trace <- ! is.null(dim(xt))

  dt <- deltat(xt)
  if ( is.na(dt) )
    dt <- 0.01
  start <- xt@positions@from

  units <- xt@units
  units.position <- xt@units.position
  if ( is.na(units) || is.null(units) ) {
    units.new <- NA
  } else {
    units.parts <- unlist(strsplit(
      sub("^([[:alpha:]]+)(/+.+)","\\1 \\2", as.character(units)),
      " ", fixed=TRUE))
    if ( length(units.parts) > 1 ) {
      secs <- unlist(strsplit(
        sub("(.+)(\\^[0-9])*","\\1 \\2", units.parts[2]),
        " ", fixed=TRUE))
      if ( grep("^\\^",secs[2]) )
        ex <- sub("\\^([0-9])*","\\1", secs[2])
      else
        ex <- nchar(secs) / 2
    } else {
      ex <- 0
    }
    if ( ex + nd > 1 )
      units.new <- sprintf("%s/s^%d", units.parts[1], ex + nd)
    else
      units.new = sprintf("%s/s", units.parts[1])
  }

  if ( multi.trace == TRUE ) {
    xt.len <- dim(xt)[1]
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    dxdt.mts = NULL
    for ( ii in 1:dim(xt)[2] ) {
      ok <- ! is.na(xt@data[,ii])
      dxdt <- .Call("CALLfd_deriv",
                    as.double(xt@data[ok,ii]), as.double(dt), as.integer(nd),
                    as.integer(order))
      if ( is.null(dxdt.mts) )
        dxdt.mts <- data.frame(dxdt)
      else
        dxdt.mts <- data.frame(dxdt.mts, dxdt)
    }
    colnames(dxdt.mts) <- colnames(xt@data)
    dxdt <- dxdt.mts
  } else {
    ok <- ! is.na(xt@data)
    xt.len <- length(xt[ok])
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    dxdt <- .Call("CALLfd_deriv",
                  as.double(xt@data[ok]), as.double(dt), as.integer(nd),
                  as.integer(order))
  }
  if ( is.finite(pct) && pct > 0 )
    dxdt <- hanning(dxdt,pct)
  dxdt <- signalSeries(dxdt, from=start, by=dt, units=units.new,
                       units.position=units.position)

  return(dxdt)
}
setMethod("deriv_fd","signalSeries",deriv_fd.signalSeries)


