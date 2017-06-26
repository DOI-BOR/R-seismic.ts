#' Rayleigh wave filtering
#'
#' \code{rayleigh.filter} determines a filter for identifying Rayleigh waves using
#' a previously computed superposition-of-ellipses decomposition.
#'
#' @param a major axis of ellipse, from \code{\link{super.ellips}}
#' @param b minor axis of ellipse, from \code{\link{super.ellips}}
#' @param dip of the ellipse plane (\code{I} from \code{\link{super.ellips}}).
#' @param strike strike of ellipse plane (\code{big-omega} from
#' \code{\link{super.ellips}}). If provided, the filter will select strike
#' angles ~ 0, which only makes sense if the original data was rotated so
#' that X=radial, and Y=transverse directions.
#' @param pitch pitch of the major axis of the ellipse (\code{little-omega}
#' from \code{\link{super.ellips}}). If provided, the filter
#' will select pitch angles ~ pi/2.
#' @param reject set to \code{TRUE} to reject Rayleigh waves, or \code{FALSE} to
#' reject everything but Rayleigh waves. Default is \code{TRUE}
#'
#' @details Computes a filter that selects or rejects Rayleigh waves based on
#' criteria defined in Pinnegar (2006). This will be a \code{vector} or
#' \code{matrix} of values between 0 and 1 that can be used to multiply the
#' original transformed 3-component time series. The filtered time series
#' typically are obtained by then inverse-transforming these products.
#' @return Filter for rejecting or selecting Rayleigh waves.  If the \code{dim}
#' attribute is set on the inputs, it will be set on the returned filter as well.
#' @seealso \code{\link{super.ellips}}
#' \itemize{
#' \item \href{http://gji.oxfordjournals.org/content/165/2/596.abstract}{Pinnegar (2006)}
#' Polarization analysis and polarization filtering of three-component signals with the
#' time-frequency S transform.
#' \item \href{http://www.bssaonline.org/content/105/1/210.abstract}{Meza-Fajardo et al (2015)}
#' Identification and Extraction of Surface Waves from Three-Component Seismograms
#' Based on the Normalized Inner Product.
#' }
#' @keywords ts

rayleigh.filter <- function(a, b, dip, strike, pitch, reject=TRUE)
{
  if ( missing(a) || missing(b) || missing(dip) )
    stop("Must provide input a, b, and dip")

  la <- length(a)
  lb <- length(b)
  ldip <- length(dip)
  if ( la != lb || la != ldip  )
    stop("Length of a, b, and dip must be the same")

  ire.selector <- b / a # aka Instaneous Reciprocal Ellipticity (IRE)
  # Rayleigh waves should have IRE > 0.5
  v1 <- 0.5 # Pinnegar uses 0.6, 0.5
  v2 <- 0.4
  F <- cos.filter(ire.selector, v1, v2)

  # Rayleigh waves should have dip - pi/2 ~ 0
  #dip.selector <- abs(dip - pi/2)
  #F <- F * cos.filter(dip.selector, pi/8, pi/4) # Pinnegar uses pi/10, pi/5

  # try to distinguish prograde vs. retrograde. may need to use strike as well
  # Rayleigh waves should have dip - pi/2 ~ 0
  motion <- "all" # can't seem to distinguish prograde vs retrograde using dip alone
  v0 <- 0
  v1 <- pi/8 # Pinnegar uses pi/10, pi/5
  v2 <- pi/4
  dip.selector <- dip - pi/2
  if ( startsWith(motion,"a") ) {
    dip.selector <- abs(dip.selector)
    v0 <- NA
  } else if ( startsWith(motion,"r") ) {
    v1 <- -v1
    v2 <- -v2
  }
  F <- F * cos.filter(dip.selector, v1, v2, v0)

  if ( ! missing(strike) ) {
    if ( length(strike) != la )
      stop("Length of a, b, dip, and strike must be the same")
    plus.x <- TRUE
    sign.x <- 1
    if ( ! plus.x )
      sign.x <- -1
    strike <- strike + 0.5 * pi * (1 - sign.x * sign(cos(strike)))
    strike <- normalize.angle(strike)

    # Rayleigh waves should have strike ~ 0
    # strike.selector <- abs(strike)
    # F <- F * cos.filter(strike.selector, pi/6, pi/3)
  }

  if ( ! missing(pitch) ) {
    if ( length(pitch) != la )
      stop("Length of a, b, dip, and pitch must be the same")
    # Rayleigh waves should have pitch - pi/2 ~ 0
    pitch.selector <- abs(pitch - pi/2)
    F <- F * cos.filter(pitch.selector, pi/8, pi/4)
  }

  if ( reject )
    F <- 1 - F

  # apply same matrix structure to F as a (if any)
  if ( ! is.null(dim(a)) )
    dim(F) <- dim(a)

  return(F)
}

cos.filter <- function(selector, v1, v2, v0=NA) {
  len <- length(selector)
  conv.filter <- as.double(rep(1, len))
  if ( v1 < v2 ) {
    if ( is.na(v0) )
      conv.filter[which(selector < v1)] <- 0
    else
      conv.filter[which(v0 <= selector & selector < v1)] <- 0
    sel <- which(v1 <= selector & selector <= v2)
    conv.filter[sel] <- 0.5 * (1 + cos(pi*(v2 - selector[sel]) / (v2 - v1)))
  } else {
    if ( is.na(v0) )
      conv.filter[which(v1 < selector)] <- 0
    else
      conv.filter[which(v1 < selector & selector <= v0)] <- 0
    sel <- which(v2 <= selector & selector <= v1)
    conv.filter[sel] <- 0.5 * (1 + cos(pi*(selector[sel] - v2) / (v1 - v2)))
  }
  return(1 - conv.filter)
}