#' Window and taper a time series.
#'
#' \code{windowTs} windows and tapers a time series.
#'
#' @param xt Equally-sampled input series. Must convert to
#' numeric vector.
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
#' @return the windowed, demeaned, and tapered data
#' @seealso \code{\link{hanning}}
#' @keywords ts

windowTs <- function(xt, dt=0.01, t0=NA, tw=NA, demean=NA, pct=NA, type=NA, norm=NA) {

	xt <- as.double(xt[!is.na(xt)])
	len <- length(xt)
	if ( len < 3 )
		stop("input time series must have at least 3 valid points")

	# call C function
	out <- .Call("CALLwindow_ts", as.double(xt), as.double(dt), as.double(t0),
							 as.double(tw), as.logical(demean), as.double(pct),
							 as.character(type), as.logical(norm))
	return(out)
}
