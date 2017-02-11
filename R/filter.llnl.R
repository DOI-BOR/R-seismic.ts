#' Filter a Univariate Time Series
#'
#' @description
#' \code{filter.llnl} filters a univariate time series using the digital
#' equivalents of common analog filter types. It is derived from the
#' Seismic Analysis Code (SAC) package from Lawrence Livermore National Labs.
#'
#' @param xt Equally-sampled input series. Must convert to numeric vector.
#' @param dt Sample interval, in seconds. Default is 0.01.
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
#' @param f.lo Low-pass filter corner frequency, in Hz. Default is 2 / (len(xt) * dt).
#' @param f.hi Hi-pass filter corner frequency, in Hz. Default is 1 / (3 * dt), which
#' is 1/3 of the Nyquist frequency.
#' @param dir Filter direction: \code{c("forward", "reverse", "zp", "both", "zerophase")}.
#' Characters are case-insensitive, and only need to uniquely define
#' the option (e.g., "F" is sufficient to select Forward). Default is zero-phase.
#' @param cheb.sb.atten Chebyshev stop band attenuation (ignored for
#' others), such that the maximum stop band amplitude is 1/atten. Default is 30.
#' @param cheb.tr.bw Chebyshev transition bandwidth between stop and pass
#' bands(ignored for others), as a fraction of the passband width. Default is 0.3.
#' @return The filtered time series.
#' @details Wraps a modified version of the SAC filter code, which in turn was
#' based on a conversion of the original Fortran code to C at the University
#' of Washington.
#' @examples
#' dirac <- c(rep(0,1000),1,rep(0,1000))
#' dt <- 0.01
#' f.lo <- 20 / (length(dirac) * dt)
#' f.hi <-  (2 / 5) * 1 / (2 * dt)
#' response <- filter.llnl(dirac, dt, order=8, pb.type="bp", filt.type="c2",
#'                        f.lo=f.lo, f.hi=f.hi, dir="zp", cheb.sb.atten=100,
#'                        cheb.tr.bw=0.4)
#' response <- ts(response, start=0, deltat=dt)
#' tsplot(response)
#' filt.spec <- spec.pgram(response,plot=FALSE)
#' ymax <- 1.2 * 10 ^ round(log10(max(filt.spec$spec)),digits=0)
#' ymin <- 0.8 * 10 ^ (round(log10(max(filt.spec$spec)),digits=0) - 10)
#' plot(filt.spec,ylab="Spectral Power",xlab="Frequency, Hz",ylim=c(ymin,ymax))
#' @seealso \href{https://ds.iris.edu/files/sac-manual/commands/bandpass.html}{SAC Manual}
#' @keywords ts
#'
filter.llnl <- function(xt, dt=0.01, order=NA, pb.type=NA, filt.type=NA,
												f.lo=NA, f.hi=NA, dir=NA, cheb.sb.atten=NA, cheb.tr.bw=NA) {

	xt <- as.double(xt[!is.na(xt)])
	len <- length(xt)
	if ( len < 3 )
		stop("input time series must have at least 3 valid points")

	# do basic sanity checking, and silently fix obviously bad values
	if ( ! is.na(order) )
		order <- max(min(order,8),1)
	if ( ! is.na(f.lo) && f.lo <= 0 )
		f.lo <- NA
	if ( ! is.na(f.hi) && f.hi >= 2 / dt )
		f.hi <- NA

	# call C function
	out <- .Call("CALLfilter_ts",
							 as.double(xt), as.double(dt), as.integer(order), as.character(pb.type),
							 as.character(filt.type), as.double(f.lo), as.double(f.hi),
							 as.character(dir), as.double(cheb.sb.atten), as.double(cheb.tr.bw))
	return(out)
}
