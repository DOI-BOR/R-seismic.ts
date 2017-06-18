#' Polarization attributes from instantaneous polarization analysis
#'
#' \code{ipa} determines polarization attributes from an input
#' 3-component time series, using the instantaneous polarization method of
#' Morozov and Smithson (1996).
#'
#' @param xt,yt,zt Equally-sampled univariate input series for the X, Y,
#' and Z directions. Alternatively, \code{xt} can be a multivariate time series
#' with X, Y, and Z data in the first 3 positions (and in that order). Input must
#' convert to numeric \code{\link{vector}}, \code{\link{matrix}},
#' \code{\link{data.frame}}, \code{\link{ts}}, or \code{\link{signalSeries}}.
#' @param dt Sample interval, in seconds. Default is 0.01, or \code{deltat}
#' if \code{xt} is a \code{\link{ts}} or \code{\link{signalSeries}} object.
#' @param demean TRUE if the means should first be removed input data. Default is TRUE.
#' @param pct Percentage of data window to apply a \code{\link{hanning}} taper to.
#' Must be between 0 and 50. Default is 0 (no taper).
#' @param w.len Smooth the input time series using a rectangular window of length
#' \code{w.len} points. Default is no smoothing.
#'
#' @details Computes the analytic signal for each component, and then finds a
#' common instantaneous phase factor across all components by maximizing a selected
#' functional. Once the common phase factor has been determined, it is easy to directly
#' compute polarization attributes such as the semi-major and semi-minor axes. This
#' method is in contrast to principal component analysis, which determines directions
#' through an eigenvalue-eigenvector decomposition of the covariance matrix of the
#' analytic signal vector. See \code{\link{pca}} for a method that yields similar
#' results using principal component analysis of the analytic signal. This method
#' provides polarization attributes for each point in time, but is independent of
#' frequency. See \code{\link{super.ellips}}, \code{\link{rayleigh.filter}} and
#' \code{\link{nip.filter}} for polarization analyses based on frequency and time.
#' See \code{\link{rayleigh.dir}} and \code{\link{pca}} for polarization estimates
#' that are independent of time and frequency.
#'
#' @return List of multicomponent time-series for the instantaneous polarization
#' attributes, including: \code{a} and \code{b} - semi-major and semi-minor ellipse
#' axes, \code{n} - unit vector normal to \code{a} and \code{b}, \code{pa} is an
#' \code{\link{mts}} of polarization attributes, including: \code{pe} - ellipticity,
#' \code{az} - azimuth of the semi-major axis, \code{plunge} - plunge of the semi-major
#' axis, \code{dip} - dip of the plane containing the \code{a} and \code{b} axes.
#' @seealso \code{\link{pca}}, \code{\link{super.ellips}},
#' \code{\link{rayleigh.filter}}, \code{\link{nip.filter}},
#' \code{\link{hanning}}, \code{\link{analytic.ts}}
#' \itemize{
#' \item \href{https://en.wikipedia.org/wiki/Analytic_signal}{Analytic Signal (Wikipedia)}
#' \item \href{http://geophysics.geoscienceworld.org/content/61/3/872.abstract}{Morozov and Smithson (1996)}
#' Instantaneous polarization attributes and directional filtering.
#' \item \href{http://gji.oxfordjournals.org/content/155/2/653.abstract}{Schimmel and Gallart (2003)}
#' The use of instantaneous polarization attributes for seismic signal detection and image enhancement.
#' }
#' @keywords ts

ipa <- function(xt, yt=NA, zt=NA, dt=NA, demean=TRUE, pct=NA, w.len=NA)
{
  if ( is.complex(xt) )
    stop("ipa: input is already complex-valued")

  # get mts representing input data
  xyz <- get.xyz(xt, yt, zt, dt, demean, pct, w.len)
  len <- dim(xyz)[1]
  dt <- deltat(xyz)

  # get the analytic signal vector
  xyz.a <- analytic.ts(xyz)

  # get phase using damped quadratic regularization
  A <- 0.5 * (xyz.a[,1]^2 + xyz.a[,2]^2 + xyz.a[,3]^2)
  B <- 0.5 * (xyz.a[,1] + xyz.a[,2] + xyz.a[,3])^2
  eps <- 0.01
  psi0 <- 0.5 * Arg(A + eps * B)

  # get the major and minor ellipse axes
  a <- Re(exp(-1i*psi0) * xyz.a) # major axis
  b <- Re(exp(-1i*(psi0 + 0.5*pi)) * xyz.a) # minor axis

  # get ts of the unit normal a x b
  n <- unorm(a, b, dt)

  # get the multicomponent signal amplitude
  At <- sqrt(a[,1]^2 + a[,2]^2 + a[,3]^2)

  # Instaneous Reciprocal Ellipticity (IRE)
  ire = sqrt(b[,1]^2 + b[,2]^2 + b[,3]^2) / sqrt(a[,1]^2 + a[,2]^2 + a[,3]^2)

  # azimuth of a-axis, measured from X axis
  az <- atan2(a[,2], a[,1])
  az <- az + 0.5 * pi * (1 - sign(cos(az))) # direction is ambiguous
  az <- normalize.angle(az)

  # get plunge of a-axis, measured from horizontal
  # and defined by tan(plunge) = (z/sqrt(x^2+y^2))
  # Vidale (1986)
  plunge <- atan2(a[,3], sqrt(a[,1]^2 + a[,2]^2))
  # normalize so that 0 < plunge <= pi/2
  plunge <- plunge + 0.5 * pi * (1 - sign(cos(plunge))) # direction is ambiguous
  plunge <- abs(normalize.angle(plunge))

  # get the strike of plane containing a and b, defined by tan(dip) = (sqrt(x^2+y^2)/-z)
  strike <- atan2(-n[,1], n[,2])

  # get the dip of plane containing a and b, defined by tan(dip) = (sqrt(x^2+y^2)/-z)
  dip <- atan2(sqrt(n[,1]^2 + n[,2]^2), -n[,3])

  # bind polarization attributes into an mts
  pa <- ts(cbind(ire=ire, az=rad2deg*az, plunge=plunge*rad2deg, strike=rad2deg*strike,
                 dip=rad2deg*dip, At=At, psi0=rad2deg*psi0), deltat=dt)

  # get weighted averages
  avg.ire <- sum(pa[,"ire"] * At) / sum(At)
  avg.az <- sum(pa[,"az"] * At) / sum(At)
  avg.plunge <- sum(pa[,"plunge"] * At) / sum(At)
  avg.strike <- sum(pa[,"strike"] * At) / sum(At)
  avg.dip <- sum(pa[,"dip"] * At) / sum(At)

  #----- return list of attributes
  ret = list(a=a, b=b, pa=pa, pa.avg=list(ire=avg.ire,az=avg.az,plunge=avg.plunge,
                                               strike=avg.strike,dip=avg.dip))

  return(ret)
}

# get the unit normal a x b
unorm <- function(a, b, dt) {
  nx <-  a[,2]*b[,3] - a[,3]*b[,2]
  ny <- -a[,1]*b[,3] + a[,3]*b[,1]
  nz <-  a[,1]*b[,2] - a[,2]*b[,1]
  nn <- sqrt(nx^2 + ny^2 + nz^2)
  nx <- nx / nn
  ny <- ny / nn
  nz <- nz / nn
  n <- ts(cbind(nx=nx,ny=ny,nz=nz), deltat=dt)

  return(n)
}
