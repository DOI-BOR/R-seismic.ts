#' Elastic response spectrum
#'
#' \code{rspec} computes the elastic response spectrum a univariate
#' time series.
#'
#' @param ts Equally-sampled input series. Must convert to numeric vector.
#' @param dt Sample interval, in seconds. Default is 0.01.
#' @param ts.type String indicating input type: \code{c("acceleration", "velocity",
#' "displacement")}. Can be abbreviated to just one character, and is case independent.
#' Default is \code{"velocity"}.
#' @param rs.type Response spectrum type: \code{c("acceleration", "velocity",
#' "displacement")} for pseudo-absolute acceleration, psuedo-relative velocity, or
#' relative displacement response spectra, respectively. Default is
#' \code{"acceleration"}.
#' @param ptap Percentage of data window to apply a Hanning taper. Must be
#' between 0 and 50. Default is 0.
#' @param damp Damping value(s) to use. Default is 0.05
#' @param tau.range Vector of length 2 defining the range of response
#' spectral periods, in seconds, to consider. Default is \code{c(0,10)}.
#' @param periods Vector of response spectral periods, in seconds, to
#' use. If not defined, a set of periods in the range tau.range is generated. The
#' periods generated are a superset of the 2008 NGA periods.
#' @param rs.meth Method to use for computing response spectra. Currently
#' implemented methods are: \enumerate{
#'  \item recursive IIR (default)
#'  \item bilinear Z transform
#'  \item original method of Niagam and Jennings (1969)
#'  \item circular convolution (slow)
#' }
#' Note that these methods may give very different results for periods < ~5 * dt.
#' @param tau.si.range Period range to compute spectrum intensity. Default
#' is \code{c(0.1,0.5)} seconds for acceleration, \code{c(0.1,2.5)} seconds
#' for velocity, and \code{c(0.1,max_period)} seconds for displacement.
#'
#' @return A \code{list} with the complete response spectra information, including:
#' damping values, spectral periods, response spectra and spectrum intensity
#' for each damping value, method, integration limits for SI, and string labels for
#' response spectra type, units, etc. The specific list returned is:
#' \code{list(damping, periods, rspect, rs.type, rs.units, rs.method,
#' 			 SI, SI.per.range, SI.units)}
#'
#' @keywords ts
rspec <- function(ts, dt=0.01, ts.type="vel", rs.type="acc", ptap=NA,
									damp=NA, tau.range=NA, periods=NA,
									rs.meth=NA, tau.si.range=NA) {
	ts <- as.double(ts[!is.na(ts)])
	if ( length(ts) < 3 )
		stop("input time series must have at least 3 points")
	out <- .Call("CALLresp",
			as.character(ts.type), as.character(rs.type),
			as.double(ts), as.double(dt), as.double(ptap), as.double(damp),
			as.double(tau.range), as.double(periods),
			as.integer(rs.meth), as.double(tau.si.range))
	return(out)
}
