#' Get the Rayleigh wave polarization direction
#'
#' \code{rayleigh.dir} estimates the Rayleigh wave polarization direction
#' from a 3-component seismic waveform.
#'
#' @param x,y,z equally-spaced seismograms in the X, Y, and Z directions.
#'
#' @details Computes rayleigh wave polarization direction using the method
#' from Meza-Fajardo et al (2015). The directions for the 3-component seismogram
#' must form a right-handed coordinate system.
#' @return The direction of polarization, in radians. The polarization direction
#' is measured (positive direction) from the positive X axis towards the positive
#' Y axis, about the positive Z axis.
#' @seealso \code{\link{super.ellips}}
#' \itemize{
#' \item \href{http://www.bssaonline.org/content/105/1/210.abstract}{Meza-Fajardo et al (2015)}
#' Identification and Extraction of Surface Waves from Three-Component Seismograms
#' Based on the Normalized Inner Product.
#' }
#' @keywords ts

rayleigh.dir <- function(x, y, z) {
  if ( missing(x) || missing(y) || missing(z) )
    stop("Must provide input x, y, and z")

  lx <- length(x)
  ly <- length(y)
  lz <- length(z)
  if ( lx != ly && lx != lz )
    stop("Length of x, y, and z must be the same")
  # get the Hilbert transform of z
  zh <- hilbert(z, op = "hilbert", zero.pad = FALSE)
  theta0 <- atan2(sum(x*zh), sum(y*zh))
  theta <- theta0 + pi * (1 - sign(sin(theta0))) +
    0.5 * pi * (1 - sign(cos(theta0))) * sign(sin(theta0))
  phi <- (0.5 * pi - theta) # positive angle measured from x towards y
  if ( phi < -pi )
    phi <- phi + 2 * pi
  if ( phi > pi )
    phi <- phi - 2 * pi
  return(list(theta=theta, phi=phi))
}