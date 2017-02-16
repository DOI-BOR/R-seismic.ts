#' Window and Taper a Multivariate Time Series.
#'
#' \code{windowTs} windows and tapers a time series.
#'
#' @param xt Equally-sampled input series, including a numeric \code{\link{vector}},
#' \code{\link{matrix}}, \code{\link{data.frame}}, \code{\link{ts}},
#' or \code{\link{signalSeries}}.
#' @param dt Sample interval, in seconds. Default is 0.01.
#' @param t0 Start time of window, in sec. Default is 0.
#' @param tw Window length, in sec. Default is len*dt - start.
#' @param demean Should windowed data be demeaned? Default is F.
#' @param pct Percentage of data window to apply a taper. Must be
#' between 0 and 50. Default is 0 (no taper).
#' @param type Taper type. Values are c("Hanning", "Bartlett", "Parzen",
#' "Blackmann-Harris", "Exact_Blackmann"). Default is "Hanning"
#' @param norm Normalize tapered data to preserve rms. Default is F
#'
#' @return the windowed, demeaned, and tapered data, , with the same object type as
#' the input.
#' @seealso \code{\link{hanning}}
#' @keywords ts

#' @describeIn windowTs.default windows a numeric \code{vector} or \code{matrix}.
windowTs.default <- function(xt, dt=NA, t0=NA, tw=NA, demean=NA,
                     pct=NA, type=NA, norm=NA) {
  if ( is.na(dt) )
    dt <- 0.01
	multi.trace <- is.matrix(xt) || length(dim(xt)) > 1
	if ( multi.trace ) {
	  xt.len <- dim(xt)[1]
	  if ( xt.len < 3 )
	    stop("input time series must have at least 3 valid points")
	  wt.mts = NULL
	  for ( ii in 1:dim(xt)[2] ) {
	    ok <- ! is.na(xt[,ii])
	    wt <- .Call("CALLwindow_ts",
	                as.double(xt[ok,ii]), as.double(dt), as.double(t0),
	                as.double(tw), as.logical(demean), as.double(pct),
	                as.character(type), as.logical(norm))
	    if ( is.null(wt.mts) )
	      wt.mts <- data.frame(wt)
	    else
	      wt.mts <- data.frame(wt.mts,wt)
	  }
	  wt <- wt.mts
	} else {
	  xt <- as.double(xt[!is.na(xt)])
	  len <- length(xt)
	  if ( len < 3 )
	    stop("input time series must have at least 3 valid points")
	  wt <- .Call("CALLwindow_ts",
	              xt, as.double(dt), as.double(t0),
	              as.double(tw), as.logical(demean), as.double(pct),
	              as.character(type), as.logical(norm))
	}

	return(wt)
}
setGeneric("windowTs",def=windowTs.default)

#' @describeIn windowTs.default windows a \code{ts} or \code{mts} object.
windowTs.ts <- function(xt, dt=NA, t0=NA, tw=NA, demean=NA,
                        pct=NA, type=NA, norm=NA) {
  multi.trace <- is.mts(xt)
  dt <- deltat(xt)
  if ( is.na(dt) )
    dt <- 0.01
  start <- start(xt)[1]

  if ( multi.trace == TRUE ) {
    xt.len <- dim(xt)[1]
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    wt.mts = NULL
    for ( ii in 1:dim(xt)[2] ) {
      ok <- ! is.na(xt[,ii])
      wt <- .Call("CALLwindow_ts",
                  as.double(xt[ok,ii]), as.double(dt), as.double(t0),
                  as.double(tw), as.logical(demean), as.double(pct),
                  as.character(type), as.logical(norm))
      if ( is.null(wt.mts) )
        wt.mts <- ts(wt, start = start, deltat = dt)
      else
        wt.mts <- ts(data.frame(wt.mts, wt), start = start, deltat = dt)
    }
    wt <- wt.mts
  } else {
    xt <- as.double(xt[!is.na(xt)])
    xt.len <- length(xt)
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    wt <- .Call("CALLwindow_ts",
                xt, as.double(dt), as.double(t0),
                as.double(tw), as.logical(demean), as.double(pct),
                as.character(type), as.logical(norm))
    wt <- ts(wt, start = start, deltat = dt)
  }
  dimnames(wt) <- dimnames(xt)

  return(wt)
}
setMethod("windowTs","ts",windowTs.ts)

#' @describeIn windowTs.default windows a \code{signalSeries} object.
windowTs.signalSeries <- function(xt, dt=NA, t0=NA, tw=NA, demean=NA,
                                  pct=NA, type=NA, norm=NA) {
  multi.trace <- ! is.null(dim(xt))

  dt <- deltat(xt)
  if ( is.na(dt) )
    dt <- 0.01
  start <- xt@positions@from
  units <- xt@units

  if ( multi.trace == TRUE ) {
    xt.len <- dim(xt)[1]
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    wt.mts = NULL
    for ( ii in 1:dim(xt)[2] ) {
      ok <- ! is.na(xt[,ii]@data)
      wt <- .Call("CALLwindow_ts",
                  as.double(xt[ok,ii]@data), as.double(dt), as.double(t0),
                  as.double(tw), as.logical(demean), as.double(pct),
                  as.character(type), as.logical(norm))
      if ( is.null(xt) )
        wt.mts <- signalSeries(wt, from = start, by = dt, units = new.units)
      else
        wt.mts <- signalSeries(data.frame(wt.mts@data, wt@data), from = start, by = dt, units = new.units)
    }
    wt <- wt.mts
  } else {
    ok <- ! is.na(xt@data)
    xt.len <- length(xt[ok])
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    wt <- .Call("CALLwindow_ts",
                as.double(xt[ok]@data), as.double(dt), as.double(t0),
                as.double(tw), as.logical(demean), as.double(pct),
                as.character(type), as.logical(norm))
    wt <- signalSeries(wt, from = start, by = dt, units = new.units)
  }
  names(wt) <- names(xt)

  return(wt)
}
setMethod("windowTs","signalSeries",windowTs.signalSeries)


