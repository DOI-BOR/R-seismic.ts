#' Elastic response spectrum
#'
#' \code{rspec} computes the elastic response spectrum a univariate
#' time series.
#'
#' @param ts Equally-sampled input series. Must convert to numeric vector.
#' @param units.ts Units of time series. One of \code{c("cgs", "mks", "SI")}.
#' Default is \code{"cgs"}.
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
#' @param rs.meth Numeric value indicating the method to use for computing response
#' spectra. Currently implemented methods are: \tabular{rl}{
#'  0 \tab recursive IIR (default) \cr
#'  1 \tab bilinear Z transform \cr
#'  2 \tab original method of Nigam and Jennings (1969) \cr
#'  3 \tab circular convolution (slow) \cr
#'  4 \tab FFT method \cr
#' }
#' These methods generally give similar results, except for periods < ~5 * dt,
#' when the results may substantially differ.
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
#' @seealso
#' \itemize{
#' \item  \code{\link{IMRot}}
#' \item  \href{https://pubs.geoscienceworld.org/ssa/bssa/article/59/2/909/116747/calculation-of-response-spectra-from-strong-motion}
#' {Nigam and Jennings (1969)} Nigam, N.C., and Jennings, P.C., 1969, Calculation of
#' response spectra from strong-motion earthquake records: Bulletin of the Seismological
#' Society of America, v. 59, no. 2, p. 909-922.
#' }
#' @keywords ts
rspec <- function(ts, dt=0.01, units.ts=NA, ts.type="vel", rs.type="acc", ptap=NA,
									damp=NA, tau.range=NA, periods=NA, rs.meth=NA, tau.si.range=NA) {
	ts <- as.double(ts[!is.na(ts)])
	if ( length(ts) < 3 )
		stop("input time series must have at least 3 points")

	g <- 9.80665 # standard gravity, in m/s^2
	m2cm <- 100
	if ( ! is.na(units.ts) ) {
	  if ( units.ts == "mks" || units.ts == "SI" || units.ts == "si" )
	    ts <- ts * m2cm
	  else if ( units.ts == "g" )
	    ts <- ts * g * m2cm
	}

	out <- .Call("CALLresp",
			as.character(ts.type), as.character(rs.type),
			as.double(ts), as.double(dt), as.double(ptap), as.double(damp),
			as.double(tau.range), as.double(periods),
			as.integer(rs.meth), as.double(tau.si.range),
			PACKAGE="seismic.ts")

	return(out)
}
