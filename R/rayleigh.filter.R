#' Rayleigh wave filtering
#'
#' \code{rayleigh.filter} determines a filter for identifying Rayleigh waves using
#' a previously computed superposition-of-ellipses decomposition.
#'
#' @param a major axis of ellipse, from \code{\link{superlips}}
#' @param b minor axis of ellipse, from \code{\link{superlips}}
#' @param I inclination (dip) of ellipse plane, from \code{\link{superlips}}
#' @param strike strike of ellipse plane (\code{big-omega} from
#' \code{\link{superlips}}). If provided, the filter will select strike
#' angles ~ 0, which only makes sense if the original data was rotated so
#' that X=radial, and Y=transverse directions.
#' @param pitch pitch of the major axis of the ellipse (\code{little-omega}
#' from \code{\link{superlips}}). If provided, the filter
#' will select pitch angles ~ pi/2.
#' @param reject.rayleigh set to \code{TRUE} to reject Rayleigh waves, or \code{FALSE} to
#' reject everything but Rayleigh waves. Default is \code{TRUE}
#'
#' @details Computes a filter that selects or rejects Rayleigh waves based on
#' criteria defined in Pinnegar (2006). This will be a \code{vector} or
#' \code{matrix} of values between 0 and 1 that can be used to multiply the
#' original transformed 3-component time series. The filtered time series
#' typically are obtained by then inverse-transforming these products.
#' @return Filter for rejecting or selecting Rayleigh waves.  If the \code{dim}
#' attribute is set on the inputs, it will be set on the returned filter as well.
#' @seealso \code{\link{superlips}}
#' \itemize{
#' \item \href{http://gji.oxfordjournals.org/content/165/2/596.abstract}{Pinnegar (2006)}
#' Polarization analysis and polarization filtering of three-component signals with the
#' time-frequency S transform.
#' \item \href{http://www.bssaonline.org/content/105/1/210.abstract}{Meza-Fajardo et al (2015)}
#' Identification and Extraction of Surface Waves from Three-Component Seismograms
#' Based on the Normalized Inner Product.
#' }
#' @keywords ts

rayleigh.filter <- function(a, b, I, strike, pitch, reject.rayleigh=TRUE)
{
  if ( missing(a) || missing(b) || missing(I) )
    stop("Must provide input a, b, and I")

  la <- length(a)
  lb <- length(b)
  lI <- length(I)
  if ( la != lb || la != lI  )
    stop("Length of a, b, and I must be the same")

  ire.selector <- b / a # aka Instaneous Reciprocal Ellipticity (IRE)
  # Rayleigh waves should have IRE < 0.5
  F <- cos.filter(ire.selector, 0.5, 0.6)

  # Rayleigh waves should have dip - pi/2 ~ 0
  dip.selector <- abs(I - pi/2)
  F <- F * cos.filter(dip.selector, pi/10, pi/5)

  if ( ! missing(strike) ) {
    if ( length(strike) != la )
      stop("Length of a, b, I, and strike must be the same")
    # Rayleigh waves should have strike ~ 0
    strike.selector <- abs(strike)
    F <- F * cos.filter(strike.selector, pi/6, pi/3)
  }

  if ( ! missing(pitch) ) {
    if ( length(pitch) != la )
      stop("Length of a, b, I, and pitch must be the same")
    # Rayleigh waves should have pitch - pi/2 ~ 0
    pitch.selector <- abs(pitch - pi/2)
    F <- F * cos.filter(pitch.selector, pi/8, pi/4)
  }

  if ( reject.rayleigh )
    F <- 1 - F

  # apply same matrix structure to F as a (if any)
  if ( ! is.null(dim(a)) )
    dim(F) <- dim(a)

  return(F)
}

cos.filter <- function(selector, v1, v2) {
  len <- length(selector)
  conv.filter <- as.double(rep(1, len))
  conv.filter[which(selector < v1)] <- 0
  sel <- which(v1 <= selector & selector <= v2)
  conv.filter[sel] <- 0.5 * (1 + cos(pi*(v2 - selector[sel]) / (v2 - v1)))
  return(1 - conv.filter)
}