#' Polarization attributes from superposition of ellipses
#'
#' \code{sea} determines instantaneous polarization attributes from an input 3-component
#' time series, using the superposition of ellipses method of Pinnegar (2006).
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
#' @details This function computes the analytic signal for each component, and then
#' uses the method of Pinnegar (2006) to compute the instantaneous semi-major and
#' semi-minor axes, as well as other elliptical elements. The elliptical elements are
#' computed by calling \code{\link{super.ellips}}. Instantaneous polarization
#' attributes are then determined from these elements. This method can be contrasted
#' to principal component analysis of Vidale (1996), which determines polarization
#' attributes through an eigenvalue-eigenvector decomposition of the covariance
#' matrix of the analytic signal vector (see \code{\link{pca}}). It can also be
#' contrasted with the instantaneous polarization attribute method (see \code{\link{ipa}}),
#' which determines elliptical elements using the method of Morozov and Smithson (1996).
#' This function provides polarization attributes for each point in time, however
#' the method also can be applied to the S-transform (see \code{\link{rayleigh.filter}}).
#'
#' @return List of multicomponent time-series for the instantaneous polarization
#' attributes, including: \code{a} and \code{b} - semi-major and semi-minor ellipse
#' axes, \code{n} - unit vector normal to \code{a} and \code{b}, \code{pa} is an
#' \code{\link{mts}} of polarization attributes, including: \code{pe} - ellipticity,
#' \code{az} - azimuth of the semi-major axis, \code{plunge} - plunge of the semi-major
#' axis, \code{dip} - dip of the plane containing the \code{a} and \code{b} axes.
#' @seealso \code{\link{analytic.ts}}, \code{\link{super.ellips}}, \code{\link{pca}},
#' \code{\link{ipa}}, \code{\link{rayleigh.filter}}, \code{\link{nip.filter}}
#' \itemize{
#' \item \href{https://en.wikipedia.org/wiki/Analytic_signal}{Analytic Signal (Wikipedia)}
#' \item \href{http://gji.oxfordjournals.org/content/165/2/596.abstract}{Pinnegar (2006)}
#' Polarization analysis and polarization filtering of three-component signals with
#' the time-frequency S transform.
#' \item \href{http://www.bssaonline.org/content/76/5/1393.abstract}{Vidale (1986)}
#' Complex polarization analysis of particle motion.
#' \item \href{http://geophysics.geoscienceworld.org/content/61/3/872.abstract}{Morozov and Smithson (1996)}
#' Instantaneous polarization attributes and directional filtering.
#' }
#' @keywords ts

sea <- function(xt, yt=NA, zt=NA, dt=NA, demean=TRUE, pct=NA, w.len=NA)
{
  if ( is.complex(xt) )
    stop("sea: input is already complex-valued")

  # get mts representing input data
  xyz <- get.xyz(xt, yt, zt, dt, demean, pct, w.len)
  len <- dim(xyz)[1]
  dt <- deltat(xyz)
  start <- start(xyz)

  # get the analytic signal vector
  xyz.a <- analytic.ts(xyz)

  # get the elliptical elements
  se <- super.ellips(xyz.a[,1], xyz.a[,2], xyz.a[,3])
  a <- se$a
  b <- se$b
  little.omega <- se$little.omega
  big.omega <- se$big.omega

  # Instaneous Reciprocal Ellipticity (IRE)
  ire = b / a

  # azimuth of a-axis, measured from X axis
  az <- se$zeta
  az <- az + 0.5 * pi * (1 - sign(cos(az))) # direction is ambiguous
  az <- normalize.angle(az)

  # get plunge of a-axis, measured from horizontal
  plunge <- se$theta
  # normalize so that 0 < plunge <= pi/2
  plunge <- plunge + 0.5 * pi * (1 - sign(cos(plunge))) # direction is ambiguous
  plunge <- abs(normalize.angle(plunge))

  # get angle of incidence
  dip <- se$I
  #dip <- dip + 0.5 * pi * (1 - sign(cos(dip))) # direction is ambiguous
  #dip <- abs(normalize.angle(dip))

  # bind polarization attributes into an mts
  pa <- ts(cbind(ire=ire, az=az*rad2deg, plunge=plunge*rad2deg,
                 strike=big.omega*rad2deg, dip=dip*rad2deg,
                 rake=little.omega*rad2deg), start=start, deltat=dt)

  # get weighted averages
  avg.ire <- sum(pa[,"ire"] * a) / sum(a)
  avg.az <- sum(pa[,"az"] * a) / sum(a)
  avg.plunge <- sum(pa[,"plunge"] * a) / sum(a)
  avg.strike <- sum(pa[,"strike"] * a) / sum(a)
  avg.dip <- sum(pa[,"dip"] * a) / sum(a)
  avg.rake <- sum(pa[,"rake"] * a) / sum(a)

  #----- return list of attributes
  ret = list(a=a, b=b, pa=pa, pa.avg=list(ire=avg.ire,az=avg.az,plunge=avg.plunge,
                                          strike=avg.strike,dip=avg.dip,rake=avg.rake))

  return(ret)
}
