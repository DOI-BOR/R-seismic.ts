#' Window and Taper a Time Series.
#'
#' \code{windowTs} windows, tapers, and demeans a univariate or multivariate time series.
#'
#' @param xt Equally-sampled input time series, including a numeric \code{\link{vector}},
#' \code{\link{matrix}}, \code{\link{data.frame}}, \code{\link{ts}}, \code{\link{mts}}
#' or \code{\link{signalSeries}}.
#' @param dt Sample interval, in seconds. Normally only used for \code{\link{vector}},
#' \code{\link{matrix}}, or \code{\link{data.frame}} input objects. If input is a
#' \code{\link{ts}}, \code{\link{mts}}, or \code{\link{signalSeries}}, then the sample
#' interval will be taken from the object attributes. However, if \code{\link{deltat}}
#' is not set, then \code{dt} will be used. Default is 0.01 seconds.
#' @param t0 Start time of window, in seconds. Default is 0.
#' @param tw Window length, in seconds. Default is \code{length(xt) * dt - t0}.
#' @param demean Set to \code{TRUE} if the windowed data should be demeaned. Default is
#' \code{FALSE}.
#' @param pct Percentage of the data window to apply a taper to. Must be
#' between 0 and 50 percent. Default is 0 (no taper).
#' @param type Taper type. Supported types are \code{c("Hanning", "Bartlett",
#' "Parzen", "Blackmann-Harris", "Exact_Blackmann")}. Case-insensitive, and only
#' need to specify enough characters to be unique. Default is \code{"Hanning"}.
#' @param norm Set to \code{TRUE} to normalize tapered data so as to preserve
#' the rms average of the signal over the window. Default is \code{FALSE}.
#'
#' @return the optionally windowed, optionally demeaned, and optionally tapered data,
#' with the same object type as the input.
#' @seealso \code{\link{ts}}, \code{\link{mts}}, and \code{\link{signalSeries}},
#' and \code{\link{hanning}}.
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
	                as.character(type), as.logical(norm),
	                PACKAGE="seismic.ts")
	    if ( is.null(wt.mts) )
	      wt.mts <- data.frame(wt)
	    else
	      wt.mts <- data.frame(wt.mts, wt)
	  }
	  colnames(wt.mts) <- colnames(xt)
	  wt <- wt.mts
	} else {
	  xt <- as.double(xt[!is.na(xt)])
	  len <- length(xt)
	  if ( len < 3 )
	    stop("input time series must have at least 3 valid points")
	  wt <- .Call("CALLwindow_ts",
	              xt, as.double(dt), as.double(t0),
	              as.double(tw), as.logical(demean), as.double(pct),
	              as.character(type), as.logical(norm),
	              PACKAGE="seismic.ts")
	}

	return(wt)
}
setGeneric("windowTs", def=windowTs.default)

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
                  as.character(type), as.logical(norm),
                  PACKAGE="seismic.ts")
      if ( is.null(wt.mts) )
        wt.mts <- data.frame(wt)
      else
        wt.mts <- data.frame(wt.mts, wt)
    }
    colnames(wt.mts) <- colnames(xt)
    wt <- wt.mts
  } else {
    xt <- as.double(xt[!is.na(xt)])
    xt.len <- length(xt)
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    wt <- .Call("CALLwindow_ts",
                xt, as.double(dt), as.double(t0),
                as.double(tw), as.logical(demean), as.double(pct),
                as.character(type), as.logical(norm),
                PACKAGE="seismic.ts")
  }
  wt <- ts(wt, start=start, deltat=dt)

  return(wt)
}
setMethod("windowTs", "ts", windowTs.ts)

#' @describeIn windowTs.default windows a \code{signalSeries} object.
windowTs.signalSeries <- function(xt, dt=NA, t0=NA, tw=NA, demean=NA,
                                  pct=NA, type=NA, norm=NA) {
  multi.trace <- ! is.null(dim(xt))

  dt <- deltat(xt)
  if ( is.na(dt) )
    dt <- 0.01
  start <- xt@positions@from
  units <- xt@units
  units.position <- xt@units.position

  if ( multi.trace == TRUE ) {
    xt.len <- dim(xt)[1]
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    wt.mts = NULL
    for ( ii in 1:dim(xt)[2] ) {
      ok <- ! is.na(xt@data[,ii])
      wt <- .Call("CALLwindow_ts",
                  as.double(xt@data[ok,ii]), as.double(dt), as.double(t0),
                  as.double(tw), as.logical(demean), as.double(pct),
                  as.character(type), as.logical(norm),
                  PACKAGE="seismic.ts")
      if ( is.null(wt.mts) )
        wt.mts <- data.frame(wt)
      else
        wt.mts <- data.frame(wt.mts, wt)
    }
    colnames(wt.mts) <- colnames(xt@data)
    wt <- wt.mts
  } else {
    ok <- ! is.na(xt@data)
    xt.len <- length(xt[ok])
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    wt <- .Call("CALLwindow_ts",
                as.double(xt@data[ok]), as.double(dt), as.double(t0),
                as.double(tw), as.logical(demean), as.double(pct),
                as.character(type), as.logical(norm),
                PACKAGE="seismic.ts")
  }
  wt <- signalSeries(wt, from=start, by=dt, units=units, units.position=units.position)

  return(wt)
}
setMethod("windowTs", "signalSeries", windowTs.signalSeries)


