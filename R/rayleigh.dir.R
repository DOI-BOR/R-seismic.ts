#' Get propagation azimuth from Rayleigh-wave polarization
#'
#' \code{rayleigh.dir} estimates the propagation azimuth from source to receiver
#' using the Rayleigh wave polarization direction from a 3-component seismic
#' waveform, assuming retrograde (default) or prograde motion.
#'
#' @param xt,yt,zt Equally-sampled univariate input series for the X, Y,
#' and Z directions. Alternatively, \code{xt} can be a multivariate time series
#' with X, Y, and Z data in the first 3 positions (and in that order). Input must
#' convert to numeric \code{\link{vector}}, \code{\link{matrix}},
#' @param plus.x \code{TRUE} if Rayleigh wave propagation angle \code{phi} is in the \code{+X}
#' half-plane (\code{-pi/2 < phi < pi/2}), and \code{FALSE} if it is in the \code{-X}
#' half-plane. Needed because this method can't distinguish between retrograde motion
#' in the \code{+X} direction and prograde motion in the \code{-X} direction. Default
#' is \code{TRUE}.
#'
#' @details Computes rayleigh wave polarization direction based on the method
#' from Meza-Fajardo et al (2015). The directions for the 3-component seismogram
#' must form a right-handed coordinate system.
#' @return The direction from the source to the receiver, relative to the positive
#' X direction, in radians.

#' @seealso \code{\link{super.ellips}}
#' \itemize{
#' \item \href{http://www.bssaonline.org/content/105/1/210.abstract}{Meza-Fajardo et al (2015)}
#' Identification and Extraction of Surface Waves from Three-Component Seismograms
#' Based on the Normalized Inner Product.
#' }
#' @keywords ts

rayleigh.dir <- function(xt, yt=NA, zt=NA, plus.x=TRUE) {
  # get mts representing the input data
  xyz <- get.xyz(xt, yt, zt, NA, FALSE, 0, 0)

  # get the Hilbert transform of z
  zh <- hilbert(xyz[,3], op = "hilbert", zero.pad = FALSE)

  # get raw source to receiver angle
  phi <- atan2(sum(xyz[,2]*zh), sum(xyz[,1]*zh))

  # Retrograde motion in the +X direction is indistinguisable from prograde
  # motion in the -X direction, and vice versa. Therefore, we need to
  # restrict the raw azimuth from its range of -pi < phi < pi to:
  # (a) -pi/2 < phi < pi/2 (propagation in +X direction), or (b) pi/2 < phi < 3*pi/2
  # (propagation in the -X direction). The sign of NIP will then indicate
  # prograde (+) or retrograde (-) motion. The function sign(cos(phi)) indicates
  # which half-plane we're in, and thus whether we need to add a factor of pi
  # to the raw phi to get to the correct quadrant
  sign.x = 1
  if ( ! plus.x )
    sign.x = -1
  phi <- phi + 0.5 * pi * (1 - sign.x * sign(cos(phi)))
  phi <- normalize.angle(phi)


  return(phi)
}