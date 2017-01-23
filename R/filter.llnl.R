#' Filter a time series using the LLNL-based functions.
#'
#' \code{filter.llnl} filters a time series.
#'
#' @param xt required (actual) equally-sampled input series. Must convert to numeric vector.
#' @param dt (optional) sample interval, in seconds. Default is 0.01.
#' @param order (optional) integer order of filter, between 1 and 10. Default is 4.
#' @param pb.type (optional) passband type:
#' band pass - c("bp","bandpass","band-pass"),
#' band reject (notch) - c("br","bandreject",band-reject","notch"),
#' low pass - c("lp","lowpass","low-pass"),
#' high pass - c("hp","highpass","high-pass").
#' CHaracters are case-insensitive, and only enough need be supplied to uniquely
#' define the string (e.g., "L" is sufficient to select low pass). Default is bandpass.
#' @param filt.type (optional) filter type:
#' Butterworth - c("bu","butterworth"),
#' Bessel - c("be","bessel"),
#' Chebyshev Type I - c("c1","Chebyshev1","Chebyshev-Type-I"),
#' Chebyshev Type II - c("c2","Chebyshev2","Chebyshev-Type-II"),
#' Characters are case-insensitive, and only enough need be supplied to uniquely
#' define the string (e.g., "BE" is sufficient to select Bessel). Default is Butterworth
#' @param f.lo (optional) low-pass filter corner frequency, in Hz. Default is 2 / (len(xt) * dt).
#' @param f.hi (optional) hi-pass filter corner frequency, in Hz. Default is 1 / (3 * dt), which
#' is 33% of the Nyquist frequency.
#' @param dir (optional) filter direction: "forward", "reverse", or
#' both (zero-phase) - c("zp","both","zerophase"). Characters are case-insensitive,
#' and only enough need be supplied to uniquely define the string (e.g., "F" is
#' sufficient to select Forward). Default is zero-phase.
#' @param cheb.sb.atten (optional) Chebyshev stop band attenuation (ignored for
#' others), such that the maximum stop band amplitude is 1/atten. Default is 30.
#' @param cheb.tr.bw (optional) Chebyshev transition bandwidth between stop and pass
#' bands(ignored for others), as a fraction of the passband width. Default is 0.3.
#
#' @return the filtered time series
filter.llnl <- function(xt, dt=0.01, order=NA, pb.type=NA, filt.type=NA,
												f.lo=NA, f.hi=NA, dir=NA, cheb.sb.atten=NA, cheb.tr.bw=NA) {

	xt <- as.double(xt[!is.na(xt)])
	len <- length(xt)
	if ( len < 3 )
		stop("input time series must have at least 3 valid points")

	# do basic sanity checking, and silently fix obviously bad values
	if ( ! is.na(order) )
		order <- max(min(order,10),1)
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
