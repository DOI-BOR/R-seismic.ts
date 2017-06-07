#' Get propagation azimuth from Rayleigh-wave polarization
#'
#' \code{rayleigh.dir} estimates the propagation azimuth from source to receiver
#' using the Rayleigh wave polarization direction from a 3-component seismic
#' waveform, assuming retrograde (default) or prograde motion.
#'
#' @param x,y,z equally-spaced seismograms in the X, Y, and Z directions.
#' @param retrograde Set to \code{TRUE} to assume retrograde motion (default).
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

rayleigh.dir <- function(x, y, z, retrograde = TRUE) {
  if ( missing(x) || missing(y) || missing(z) )
    stop("Must provide input x, y, and z")

  lx <- length(x)
  ly <- length(y)
  lz <- length(z)
  if ( lx != ly && lx != lz )
    stop("Length of x, y, and z must be the same")
  # get the Hilbert transform of z
  zh <- hilbert(z, op = "hilbert", zero.pad = FALSE)

  # get source to receiver angle (azimuth relative to positive X-direction)
  phi <- atan2(sum(y*zh), sum(x*zh))
  if ( retrograde )
    phi <- normalize.angle(phi - pi)

  return(phi)
}