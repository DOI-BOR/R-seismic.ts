#' Filter a Multivariate Time Series
#'
#' @description
#' \code{filter.llnl} filters a univariate or multivariate time series using
#' the digital  equivalents of common analog filter types. It is derived from the
#' Seismic Analysis Code (SAC) package from Lawrence Livermore National Labs.
#'
#' @param xt Equally-sampled input time series, including a numeric \code{\link{vector}},
#' \code{\link{matrix}}, \code{\link{data.frame}}, \code{\link{ts}}, \code{\link{mts}}
#' or \code{\link{signalSeries}}.
#' @param dt Sample interval, in seconds. Default is 0.01 if input is not a
#' \code{\link{ts}} or \code{\link{signalSeries}}.
#' @param order Order of filter, between 1 and 8. Default is 4.
#' @param pb.type \describe{
#' \item{band pass}{\code{c("bp","bandpass","band-pass")}}
#' \item{band reject (notch)}{\code{c("br","bandreject","band-reject","notch")}}
#' \item{low pass}{\code{c("lp","lowpass","low-pass")}}
#' \item{high pass}{\code{c("hp","highpass","high-pass")}}
#' }
#' Characters are case-insensitive, and only need to uniquely
#' define the option (e.g., "L" is sufficient to select low pass). Default is
#' bandpass.
#' @param filt.type \describe{
#' \item{Butterworth}{\code{c("bu","butterworth")}}
#' \item{Bessel}{\code{c("be","bessel")}}
#' \item{Chebyshev Type I}{\code{c("c1","Chebyshev1","Chebyshev-Type-I")}}
#' \item{Chebyshev Type II}{\code{c("c2","Chebyshev2","Chebyshev-Type-II")}}
#' }
#' Characters are case-insensitive, and only need to uniquely
#' define the option (e.g., "BE" is sufficient to select Bessel). Default is
#' Butterworth.
#' @param f.lo Low-cut filter corner frequency, in Hz. Default is 2 / (len(xt) * dt).
#' @param f.hi Hi-cut filter corner frequency, in Hz. Default is 1 / (3 * dt), which
#' is 1/3 of the Nyquist frequency.
#' @param dir Filter direction: \code{c("forward", "reverse", "zp", "both", "zerophase")}.
#' Characters are case-insensitive, and only need to uniquely define
#' the option (e.g., "F" is sufficient to select Forward). Default is zero-phase.
#' @param cheb.sb.atten Chebyshev stop band attenuation (ignored for
#' others), such that the maximum stop band amplitude is 1/atten. Default is 30.
#' @param cheb.tr.bw Chebyshev transition bandwidth between stop and pass
#' bands(ignored for others), as a fraction of the passband width. Default is 0.3.
#' @return The filtered time series, with the same object type as the input.
#' @details Wraps a modified version of the SAC filter code, which in turn was
#' based on a conversion of the original Fortran code to C at the University
#' of Washington.
#' @examples
#' dt <- 0.01
#' dirac <- ts(c(rep(0,1000),1,rep(0,1000)), start=0, deltat=dt)
#' f.lo <- 20 / (length(dirac) * dt)
#' f.hi <-  (2 / 5) * 1 / (2 * dt)
#' response <- filter.llnl(dirac, order=8, pb.type="bp", filt.type="c2",
#'                        f.lo=f.lo, f.hi=f.hi, dir="zp", cheb.sb.atten=100,
#'                        cheb.tr.bw=0.4)
#' tsplot(response) # response is a ts object because input was
#' filt.spec <- spec.pgram(response,plot=FALSE)
#' ymax <- 1.2 * 10 ^ round(log10(max(filt.spec$spec)),digits=0)
#' ymin <- 0.8 * 10 ^ (round(log10(max(filt.spec$spec)),digits=0) - 10)
#' plot(filt.spec,ylab="Spectral Power",xlab="Frequency, Hz",ylim=c(ymin,ymax))
#' @seealso \href{https://ds.iris.edu/files/sac-manual/commands/bandpass.html}{SAC Manual}
#' @keywords ts
#'

#' @describeIn filter.llnl.default filters a numeric \code{vector} or \code{matrix}.
filter.llnl.default <- function(xt, dt=NA, order=NA, pb.type=NA, filt.type=NA,
												f.lo=NA, f.hi=NA, dir=NA, cheb.sb.atten=NA, cheb.tr.bw=NA) {
	multi.trace <- is.matrix(xt) || length(dim(xt)) > 1
	if ( is.na(dt) )
	  dt <- 0.01

	# do basic sanity checking, and silently fix obviously bad values
	if ( ! is.na(order) )
		order <- max(min(order,8),1)
	if ( ! is.na(f.lo) && f.lo <= 0 )
		f.lo <- NA
	if ( ! is.na(f.hi) && f.hi >= 2 / dt )
		f.hi <- NA

	if ( multi.trace ) {
	  xt.len <- dim(xt)[1]
	  if ( xt.len < 3 )
	    stop("input time series must have at least 3 valid points")
	  ft.mts = NULL
	  for ( ii in 1:dim(xt)[2] ) {
	    ok <- ! is.na(xt[,ii])
	    ft <- .Call("CALLfilter_ts",
	                 as.double(xt[ok,ii]), as.double(dt), as.integer(order), as.character(pb.type),
	                 as.character(filt.type), as.double(f.lo), as.double(f.hi),
	                 as.character(dir), as.double(cheb.sb.atten), as.double(cheb.tr.bw))
	    if ( is.null(dxdt.mts) )
	      ft.mts <- data.frame(ft)
	    else
	      ft.mts <- data.frame(ft.mts, ft)
	  }
	  colnames(ft.mts) <- colnames(xt)
	  ft <- ft.mts
	} else {
	  xt <- as.double(xt[!is.na(xt)])
	  xt.len <- length(xt)
	  if ( xt.len < 3 )
	    stop("input time series must have at least 3 valid points")
	  ft <- .Call("CALLfilter_ts",
	              xt, as.double(dt), as.integer(order), as.character(pb.type),
	              as.character(filt.type), as.double(f.lo), as.double(f.hi),
	              as.character(dir), as.double(cheb.sb.atten), as.double(cheb.tr.bw))
	}

	return(ft)
}
setGeneric("filter.llnl",def=filter.llnl.default)

#' @describeIn filter.llnl.default filters a \code{ts} or \code{mts} object.
filter.llnl.ts <- function(xt, dt=NA, order=NA, pb.type=NA, filt.type=NA,
                           f.lo=NA, f.hi=NA, dir=NA, cheb.sb.atten=NA, cheb.tr.bw=NA) {
  multi.trace <- is.mts(xt)
  dt <- deltat(xt)
  if ( is.na(dt) )
    dt <- 0.01
  start <- start(xt)[1]

  # do basic sanity checking, and silently fix obviously bad values
  if ( ! is.na(order) )
    order <- max(min(order,8),1)
  if ( ! is.na(f.lo) && f.lo <= 0 )
    f.lo <- NA
  if ( ! is.na(f.hi) && f.hi >= 2 / dt )
    f.hi <- NA

  if ( multi.trace == TRUE ) {
    xt.len <- dim(xt)[1]
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    ft.mts = NULL
    for ( ii in 1:dim(xt)[2] ) {
      ok <- ! is.na(xt[,ii])
      ft <- .Call("CALLfilter_ts",
                  as.double(xt[ok,ii]), as.double(dt), as.integer(order), as.character(pb.type),
                  as.character(filt.type), as.double(f.lo), as.double(f.hi),
                  as.character(dir), as.double(cheb.sb.atten), as.double(cheb.tr.bw))
      if ( is.null(ft.mts) )
        ft.mts <- data.frame(ft)
      else
        ft.mts <- data.frame(ft.mts, ft)
    }
    colnames(ft.mts) <- colnames(xt)
    ft <- ft.mts
  } else {
    xt <- as.double(xt[!is.na(xt)])
    xt.len <- length(xt)
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    ft <- .Call("CALLfilter_ts",
                xt, as.double(dt), as.integer(order), as.character(pb.type),
                as.character(filt.type), as.double(f.lo), as.double(f.hi),
                as.character(dir), as.double(cheb.sb.atten), as.double(cheb.tr.bw))
  }
  ft <- ts(ft, start=start, deltat=dt)

  return(ft)
}
setMethod("filter.llnl","ts",filter.llnl.ts)

#' @describeIn filter.llnl.default filters a \code{signalSeries} object.
filter.llnl.signalSeries <- function(xt, dt=NA, order=NA, pb.type=NA, filt.type=NA,
                                     f.lo=NA, f.hi=NA, dir=NA, cheb.sb.atten=NA, cheb.tr.bw=NA) {
  multi.trace <- ! is.null(dim(xt))
  if ( is.na(dt) )
    dt <- deltat(xt)
  start <- xt@positions@from
  units <- xt@units
  units.position <- xt@units.position

  # do basic sanity checking, and silently fix obviously bad values
  if ( ! is.na(order) )
    order <- max(min(order,8),1)
  if ( ! is.na(f.lo) && f.lo <= 0 )
    f.lo <- NA
  if ( ! is.na(f.hi) && f.hi >= 2 / dt )
    f.hi <- NA

  if ( multi.trace == TRUE ) {
    xt.len <- dim(xt)[1]
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    ft.mts = NULL
    for ( ii in 1:dim(xt)[2] ) {
      ok <- ! is.na(xt@data[,ii])
      ft <- .Call("CALLfilter_ts",
                  as.double(xt[ok,ii]@data), as.double(dt), as.integer(order), as.character(pb.type),
                  as.character(filt.type), as.double(f.lo), as.double(f.hi),
                  as.character(dir), as.double(cheb.sb.atten), as.double(cheb.tr.bw))
      if ( is.null(ft.mts) )
        ft.mts <- data.frame(ft)
      else
        ft.mts <- data.frame(ft.mts, ft)
    }
    colnames(ft.mts) <- colnames(xt@data)
    ft <- ft.mts
  } else {
    ok <- ! is.na(xt@data)
    xt.len <- length(xt[ok])
    if ( xt.len < 3 )
      stop("input time series must have at least 3 valid points")
    ft <- .Call("CALLfilter_ts",
                as.double(xt[ok]@data), as.double(dt), as.integer(order), as.character(pb.type),
                as.character(filt.type), as.double(f.lo), as.double(f.hi),
                as.character(dir), as.double(cheb.sb.atten), as.double(cheb.tr.bw))
  }
  names(ft) <- names(xt)
  ft <- signalSeries(ft, from=start, by=dt, units=units, units.position=units.position)

  return(ft)
}
setMethod("filter.llnl","signalSeries",filter.llnl.signalSeries)


