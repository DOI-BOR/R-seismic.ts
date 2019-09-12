#' Combine and align univariate time series into a multivariate time series
#'
#' @description
#' \code{combine.traces} combines a set of real, univariate or multivariate
#' time series into a multivariate time series. The traces are aligned using
#' cross-correlations methods. Reversed-polarity traces are inverted.
#'
#' @param ... An input set of equally-sampled univariate or multivariate time
#' series, each of which can be represented as a \code{\link{matrix}},
#' \code{\link{data.frame}}, \code{\link{signalSeries}}, \code{\link{vector}},
#' \code{\link{ts}} or \code{\link{mts}} object. Different types can be mixed.
#' The first time series will be considered the reference time series.
#' @param aug.factor Augmentation factor. Allows for sub-sample alignment
#' of input traces if greater than 1. For the cross-correlations, the traces
#' are interpolated to a sample rate of N / dt by zero-padding the FFT. The
#' frequency content is unchanged. 10 is often a reasonable number. Default
#' is 1 (no augmentation). The returned traces are decimated by the same factor.
#' See \code{\link{augment}}.
#' @param ref.lag An integer number of samples to delay (positive) or advance
#' (negative) traces relative to reference trace. Used for optimizing searching.
#' See \code{\link{lag.ccf} xy.lag}.
#' @param w.beg Beginning integer sample index of the comparison window. Default
#' is the first sample. See \code{\link{lag.ccf}}.
#' @param w.len Length, in samples, of the comparison window. Default is
#' the full window length (after accounting for any lags). See \code{\link{lag.ccf}}.
#' @param pct Percentage of data window to apply a Hanning taper.
#' Default is 0. See \code{\link{lag.ccf}}.
#' @param lag.max Maximum lag of the maximum lag for which to compute the
#' autocorrelation. See \code{\link{lag.ccf}}.
#' @param ccf.sig Least significant value of CCF to be used for aligning
#' traces. Must be be between 0 and 1. Traces will not be aligned unless
#' the max (or -min) of the CCF is greater than this. Default is 0.5.
#' @param ccf.rev For detecting reversed-ploarity traces, the ratio of the
#' minimum CCF to the maximum CCF must be greater than \code{ccf.rev}.
#' Default is 1.5
#' @param same.lag Use the lag determined for the first trace after the
#' reference trace for all other traces
#' @param no.dec Don't decimate the returned traces if the traces have been
#' augmented. Default is the returned traces will be decimated if augmented
#' @param id An optional id (numeric or string) prepended to messages. Default is "".
#' @param col.names A vector of names for the combined traces. If not specified,
#' then the name is determined from the name specified on input, or if none,
#' a system-generated name.
#' @return Returns \code{list(xt.all, xt.lags, xt.revs, ccfs, aug.factor, no.dec)}
#' where \code{xt.all} is a \code{\link{matrix}} of the combined time series (each
#' column is a univariate time series), \code{xt.lags} is a vector of the
#' lags, in augmented samples, relative to the reference trace, \code{xt.revs}
#' is a vector of boolean values indicating if the trace was detected to have
#' reversed polarity, \code{ccfs} is the raw CCF output from
#' \code{\link{lag.ccf}}, and the other values are the input parameters.
#' @details Input time series must have the same sample rate, but can be
#' of varing length. Both univariate and multivariate time series can be
#' input. Multivariate time series will be split into their
#' univariate components before aligning and combining.
#' @seealso \code{\link{lag.ccf}}, \code{\link{augment}}, \code{\link{decimate}}
#' @keywords ts

combine.traces <- function (..., demean=FALSE, aug.factor=NA, ref.lag=NA,
                            w.beg=NA, w.len=NA, pct=NA, lag.max=NA,
                            ccf.sig=NA, ccf.rev=NA, same.lag=FALSE, no.dec=FALSE,
                            id=NA, col.names=NA ) {
  n.opt.args <- ...length()
  if ( n.opt.args < 1 )
    stop("No time series specified")

  if ( is.na(aug.factor) )
    aug.factor <- 1
  if ( is.na(ccf.sig) )
    ccf.sig <- 0.5
  if ( is.na(ccf.rev) )
    ccf.rev <- 1.5
  if ( is.na(id) )
    id <- ""
  if ( is.na(same.lag) )
    same.lag <- FALSE
  if ( is.na(no.dec) )
    no.dec <- FALSE

  # get the optional args as language
  opt.args <- eval(substitute(alist(...)))
  opt.args.names <- names(opt.args)

  # loop over the univariate or multivariate time series specified in ...
  # and get the lags relative to the first (reference) trace
  xt.ref <- NULL
  lags <- list()
  zt.names <- NULL
  zt.lengths <- NULL
  for ( ii in 1:n.opt.args ) {
    xt.name <- if ( is.null(opt.args.names) || opt.args.names[ii] == "" )
        paste0("xt.",ii) else opt.args.names[ii]
    xt <- ...elt(ii) # eval(opt.args[[ii]]) won't work here
    multi.trace <- is.mts(xt) ||
      ( ! is.null(dim(xt)) && length(dim(xt)) > 1 && dim(xt)[2] > 1 )
    if ( multi.trace ) {
      xt.len <- dim(xt)[1]
      for ( jj in 1:dim(xt)[2] ) {
        if ( is(xt, "signalSeries") ) {
          zt <- xt@data[,jj]
          zt.name <- colnames(xt@data)
        } else {
          zt <- xt[,jj]
          zt.name <- colnames(xt)
        }
        if ( is.null(zt.name) )
          zt.name <- paste0(xt.name,".",jj)
        if ( demean )
          zt <- windowTs(zt,demean=TRUE)
        if ( is.null(xt.ref) ) {
          xt.ref <- zt
        } else {
          # get cross-correlations with reference trace
          lag.xz <- lag.ccf(xt.ref, zt, xy.lag=ref.lag, w.beg=w.beg, w.len=w.len,
                  pct=pct, lag.max=lag.max, aug.factor=aug.factor)
          lags[[length(lags)+1]] <- lag.xz
        }
        zt.names <- c(zt.names, zt.name)
        zt.lengths <- c(zt.lengths, length(zt) * aug.factor)
      }
    } else {
      if ( is(xt, "signalSeries") ) {
        zt <- xt@data
        zt.name <- colnames(xt@data)
      } else {
        zt <- xt
        zt.name <- colnames(xt)
      }
      if ( is.null(zt.name) )
        zt.name <- xt.name
      if ( demean )
        zt <- windowTs(zt,demean=TRUE)
      if ( is.null(xt.ref) ) {
        xt.ref <- zt
      } else {
        # get cross-correlations with reference trace
        lag.xz <- lag.ccf(xt.ref, zt, xy.lag=ref.lag, w.beg=w.beg, w.len=w.len,
                          pct=pct, lag.max=lag.max, aug.factor=aug.factor)
        lags[[length(lags)+1]] <- lag.xz
      }
      zt.names <- c(zt.names, zt.name)
      zt.lengths <- c(zt.lengths, length(zt) * aug.factor)
    }
  }

  # find the minimum lag
  maxCCF.lag.el <- 0
  min.lag <- NA
  zt.lags <- NULL
  zt.revs <- NULL
  n.lags <- length(lags)
  for ( ii in 1:n.lags ) {
    lag.xz <- lags[[ii]]
    zt.name <- zt.names[ii+1]
    zt.rev <- FALSE
    if ( -lag.xz$min.ccf < ccf.sig && lag.xz$max.ccf < ccf.sig ) {
      zt.lag <- 0
      warning(sprintf("%s %s: CCF below minimum significance level", id, zt.name),
              call.=FALSE)
      print(lag.xz)
    } else {
      if ( -lag.xz$min.ccf > ccf.rev * lag.xz$max.ccf ) {
        zt.lag <- lag.xz$min.lag
        warning(sprintf("%s %s: polarity is reversed: inverting", id, zt.name),
                call.=FALSE)
        zt.rev <- TRUE
      } else {
        zt.lag <- lag.xz$max.lag
      }
    }
    if ( same.lag && ii > 1 )
      zt.lag <- zt.lags[1]
    zt.lags <- c(zt.lags, zt.lag)
    zt.revs <- c(zt.revs, zt.rev)
    if ( is.na(min.lag) || zt.lag < min.lag ) {
      maxCCF.lag.el <- ii
      min.lag <- zt.lag
    }
  }

  # find the min length
  min.lag <- min(min.lag, 0)
  n.tr <- length(zt.names)
  # print(sprintf("min.lag=%d n.lag=%d n.ztlags=%d n.tr=%d ",min.lag,n.lags,length(zt.lags),n.tr))
  min.len <- NA
  for ( ii in 1:n.tr ) {
    lag.ref <- if ( ii == 1 ) 0 else zt.lags[ii-1]
    tr.len <- zt.lengths[ii] + min.lag - lag.ref
    # print(sprintf("ii=%d lag.ref=%d zt.length=%d tr.len=%d",
    #               ii,lag.ref,zt.lengths[ii],tr.len))
    if ( is.na(min.len) || tr.len < min.len )
      min.len <- tr.len
  }


  # loop over the univariate or multivariate time series specified in ...
  # and get the traces lagged relative to the minimum lag trace and
  # reference trace
  xt.all <- NULL
  kk <- 1
  for ( ii in 1:n.opt.args ) {
    xt <- ...elt(ii) # eval(opt.args[[ii]]) won't work here
    multi.trace <- is.mts(xt) ||
      ( ! is.null(dim(xt)) && length(dim(xt)) > 1 && dim(xt)[2] > 1 )
    if ( multi.trace ) {
      xt.len <- dim(xt)[1]
      for ( jj in 1:dim(xt)[2] ) {
        lag.ref <- if ( kk == 1 ) 0 else zt.lags[kk+1]
        t.s <- 1 + lag.ref
        t.e <- min.len + lag.ref
        zt <- if ( is(xt, "signalSeries") ) xt@data[,jj] else xt[,jj]
        if ( kk > 1 && zt.revs[kk-1] )
          zt <- -zt
        if ( aug.factor > 1 )
          zt <- augment(zt, aug.factor)
        zt <- zt[t.s:t.e]
        if ( demean )
          zt <- windowTs(zt,demean=TRUE)
        xt.all <- cbind(xt.all,zt)
      }
    } else {
      lag.ref <- if ( kk == 1 ) 0 else zt.lags[kk-1]
      t.s <- 1 + lag.ref - min.lag
      t.e <- min.len + lag.ref - min.lag
      # print(sprintf("kk=%d lag.ref=%d min.len=%d t.s=%d t.e=%d zt.len=%d",
      #               kk,lag.ref,min.len,t.s,t.e,zt.lengths[kk]))
      zt <- if ( is(xt, "signalSeries") ) xt@data else xt
      if ( kk > 1 && zt.revs[kk-1] )
        zt <- -zt
      if ( aug.factor > 1 )
        zt <- augment(zt, aug.factor)
      zt <- zt[t.s:t.e]
      if ( demean )
        zt <- windowTs(zt,demean=TRUE)
      xt.all <- cbind(xt.all,zt)
    }
    kk <- kk + 1
  }

  # decimate if augmented. set start so as to get original ref points
  if ( aug.factor > 1 && ! no.dec ) {
    t.s <- 1
    if ( min.lag < 0 )
      t.s <- aug.factor - ( -min.lag ) %% aug.factor + 1
    xt.all <- decimate(xt.all, factor=aug.factor, first.t=t.s, lp.filter=FALSE)
  }

  colnames(xt.all) <- if ( anyNA(col.names) || n.tr != length(col.names) )
                          zt.names else col.names

  ret <- list(xt.all=xt.all, xt.lags=zt.lags, xt.revs=zt.revs,
              aug.factor=aug.factor, no.dec=no.dec, ccfs=lags)
  return (ret)
}
