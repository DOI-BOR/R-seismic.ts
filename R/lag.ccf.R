#' Find the maximum and minimum of the cross-correlation function between two
#' univariate equally-spaced series. The time series can optionally be
#' augmented and/or windowed beforehand.
#'
#' @description
#' \code{lag.ccf} finds the lag of the maximum and minimum of the
#' cross-correlation function between two univariate series with uniform sample
#' interval. The sign of the returned lag gives the offset, in samples, of
#' the second series relative to the first. For example if x is a series, and
#' y an identical copy - but shifted to the right (left) by 10 (-10)
#' samples - then this function will return 10 (-10) for the lag of the
#' maximum. On the other hand, if y is a reversed polarity copy of x, i.e.,
#' \code{y = -x}, then 10 (-10) is returned for the lag of the minimum.
#'
#' @usage lag.ccf(
#'  xt, yt,
#'  xy.lag = NA,
#'  w.beg = NA,
#'  w.len = NA,
#'  pct = NA,
#'  lag.max = NA,
#'  aug.factor = NA )
#' @param xt required (first) equally-spaced series. Must convert to numeric vector.
#' @param yt required (second) equally-spaced series. Must convert to numeric vector.
#' @param xy.lag (optional) integer number of samples to delay (positive) or
#' advance (negative) yt relative to xt, used for optimizing searching. If the
#' lag is large, but is is known approximately, then this can be used as a hint
#' that will be faster then setting \code{lag.max} to a large value. Default is 0.
#' The returned lags will include the value of \code{xy.lag}.
#' @param w.beg (optional) beginning integer sample index of the comparison window.
#' The window is defined on the data after it has been lagged with \code{xy.lag},
#' and the non-overlapping parts of the vectors have been discarded. Default is
#' the first data point.
#' @param w.len (optional) length, in samples, of the comparison window. Default is
#' the full window length (after accounting for any lags).
#' @param pct (optional) percentage of data window to apply a Hanning taper. Must be
#' between 0 and 50. Default is 0.
#' @param lag.max (optional) maximum lag for which to compute the autocorrelation.
#' This effectively defines the maximum lag that can be found, but larger values
#' are slower. Default is the lesser of \code{w.len} or \code{30*log10(w.len)}.
#' @param aug.factor Augmentation factor (see \code{\link{augment}}). Default is NA.
#' @return Returns the list \code{(max.lag, max.ccf, min.lag, min.ccf, w.len,
#' lag.max, aug.factor)}
#' with the lags and values of the maximum an d minimum of the CCF, as well as
#' the lengths of the comparison window and maximum lag considered. The returned
#' lags incorporate whatever value \code{xy.lag} is set to, and therefore represent
#' the lags of the unshifted time series.
#' @keywords ts cross-correlation
#' @export lag.ccf
lag.ccf <- function(xt, yt, xy.lag = NA, w.beg = NA, w.len = NA,
                      pct = NA, lag.max = NA, aug.factor = NA) {
  # xy.lag is the delay, in samples, to apply to series y. A positive lag results in
  #		delaying series y relative to series x. A negative lag results in advancing
  #		series y relative to series x. Default lag is 0.
  if ( is.na(xy.lag) )
    xy.lag <- 0
  if ( xy.lag >= 0 ) {
    xt <- as.vector(xt[!is.na(xt)])
    yt <- as.vector(yt[!is.na(yt)])
    xy.swapped = FALSE
  } else {
    # to avoid conditionals later, we'll just redefine the series names
    #		so that lag is always non-negative (yt always delayed relative to xt)
    tmp <- xt[!is.na(xt)]
    xt <- as.vector(yt[!is.na(yt)])
    yt <- as.vector(tmp)
    xy.lag <- -xy.lag
    xy.swapped = TRUE
  }
  xt.len <- length(xt)
  yt.len <- length(yt)
  if ( yt.len - xy.lag <= 0 )
    stop(sprintf("lag (%d) >= length (%d)", xy.lag, yt.len))

  # set or check optional window parameters (default is max available window)
  # set or get w.beg, relative to current xt
  if ( is.na(w.beg) )
    w.beg <- 1
  else if ( w.beg <= 0 )
    stop("w.beg must be positive")
  else {
    # w.beg specified. Note that w.beg is assumed to be referenced w.r.t
    #	the original xt, so we need to adjust if swapped.
    if ( xy.swapped ) {
      w.beg <- w.beg - xy.lag
      if ( w.beg <= 0 )
        stop(sprintf("w.beg (%d) <= -xy.lag (%d)",w.beg,xy.lag))
    }
    if ( w.beg > xt.len )
      stop(sprintf("w.beg (%d) > length (%d)",
                   w.beg, xt.len))
    else if ( xy.lag + w.beg > yt.len )
      stop(sprintf("xy.lag (%d) + w.beg (%d) > length (%d)",
                   xy.lag, w.beg, yt.len))
  }
  if ( is.na(w.len) )
    w.len <- min(xt.len, yt.len - xy.lag) - w.beg + 1
  else if ( w.beg + w.len - 1 > xt.len )
    stop(sprintf("w.beg (%d) + w.len (%d) - 1 > length (%d)",
                 w.beg, w.len, xt.len))
  else if ( xy.lag + w.beg + w.len - 1 > yt.len )
    stop(sprintf("xy.lag (%d) + w.beg (%d) + w.len (%d) - 1 > length (%d)",
                 xy.lag, w.beg, w.len, yt.len))
  w.end <- w.beg + w.len - 1 # relative to current xt

  # make an mts with the windowed (possibly swapped) data
  xy.ts <- ts(cbind(xt[w.beg:w.end], yt[(xy.lag+w.beg):(xy.lag+w.end)]))

  # taper the windowed data
  if ( ! is.na(pct) && pct > 0 )
    xy.ts <- windowTs(xy.ts, dt=NA, type="Hanning", pct=pct)

  # augment the windowed and tapered data if requested
  if ( is.na(aug.factor) )
    aug.factor <- 1
  if ( aug.factor > 1 ) {
    xy.ts <- augment(xy.ts, aug.factor)
    w.len <- w.len * aug.factor
    if ( ! is.na(xy.lag) )
      xy.lag <- xy.lag * aug.factor
    if ( ! is.na(lag.max) )
      lag.max <- lag.max * aug.factor
  }

  # if not set, set lag.max = max lag at which to compute the acf. If this is
  # too long, the acf will be slow
  if ( is.na(lag.max) )
    lag.max <- max(1, 30*log10(w.len))
  # lag.max must not exceed the window length
  lag.max <- round(min(w.len, lag.max))

  # get the auto and cross correlations
  # xz.acvf <- acf(xy.ts, type="correlation", lag.max=lag.max, plot=TRUE) # debug
  xy.acvf <- acf(xy.ts, type="correlation", lag.max=lag.max, plot=FALSE)

  # get the lag of the max value from the 2 parts of the cross correlation
  # returned by acf
  acf.21.max <- max(xy.acvf$acf[,2,1])
  acf.12.max <- max(xy.acvf$acf[,1,2])
  if ( acf.12.max > acf.21.max ) {
    max.lag <- -which(xy.acvf$acf[,1,2] == acf.12.max) + 1
    max.ccf <- acf.12.max
    #print(sprintf("lagMaxCCF: (-) w.len=%d lag.max=%d xy.lag=%d max.lag=%d+%d",
    #              w.len,lag.max,xy.lag,max.lag,xy.lag))
  } else {
    max.lag <- which(xy.acvf$acf[,2,1] == acf.21.max) - 1
    max.ccf <- acf.21.max
    #print(sprintf("lagMaxCCF: (+) w.len=%d lag.max=%d xy.lag=%d max.lag=%d+%d",
    #              w.len,lag.max,xy.lag,max.lag,xy.lag))
  }
  if ( ! is.na(xy.lag) )
    max.lag <- max.lag + xy.lag

  # get the lag of the min value of the cross correlation (max negative)
  acf.21.min <- min(xy.acvf$acf[,2,1])
  acf.12.min <- min(xy.acvf$acf[,1,2])
  if ( acf.12.min < acf.21.min ) {
    min.lag <- -which(xy.acvf$acf[,1,2] == acf.12.min) + 1
    min.ccf <- acf.12.min
  } else {
    min.lag <- which(xy.acvf$acf[,2,1] == acf.21.min) - 1
    min.ccf <- acf.21.min
  }
  if ( ! is.na(xy.lag) )
    min.lag <- min.lag + xy.lag

  if ( xy.swapped ) {
    max.lag <- -max.lag
    min.lag <- -min.lag
  }

  results <- list(max.lag=max.lag, max.ccf=max.ccf, min.lag=min.lag,
                  min.ccf=min.ccf, w.len=w.len, lag.max=lag.max,
                  aug.factor=aug.factor)
  return(results)
}
