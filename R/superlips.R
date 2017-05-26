#' 3D complex waveform decomposition using superposition of ellipses
#'
#' \code{superlips} decomposes the complex transform of a 3-component waveform
#' into an equivalent superposition of ellipses, from which polarization
#' attributes can be determined. Typically, inputs are the complex Fourier
#' or S Transforms of a 3-component seismogram.
#'
#' @param X,Y,Z Equal-length \code{complex} vectors or matrices for the X, Y,
#' and Z directions representing the complex Fourier or S Transforms
#' of the original 3-component data.
#' @param nrow,ncol If specified, and \code{X, Y} and \code{Z} are matrices,
#' then convert the ouput elliptical parameters to matrix form, with
#' \code{\link{dim} = c(nrow,ncol)}, and filled by row.
#'
#' @details The code implements the methods of Pinnegar (2006) to compute the
#' following elliptical elements for each input point:
#' \describe{
#' \item{\code{a}}{Semi-major axis of the ellipse. \code{a >= 0}}
#' \item{\code{b}}{Semi-minor axis of the ellipse. \code{a >= b >= 0}}
#' \item{\code{I}}{Inclination (dip) of the ellipse to the horizontal (X-Y) plane.
#' \code{0 < I < pi}}
#' \item{\code{big-omega}}{Azimuth (strike) of ascending node. \code{-pi < big-omega < pi}}
#' \item{\code{little-omega}}{Pitch angle between the ascending node and the
#' semi-major axis. \code{0 < little-omega < pi}}
#' \item{\code{phi}}{Phase, with respect to the maximum direction of motion.
#' \code{-pi < phi < pi}}
#' }
#' See Pinnegar (2006) for details. The elliptical elements can be used to
#' compute polarization attributes. For example, the Instaneous Reciprocal
#' Ellipticity (IRE) is determined from \code{MRE = b / a}. Filters also can
#' be determined from the elemts to select or reject Rayleigh waves
#' (e.g., see \code{\link{rayleigh.filter}})
#' @return List with the elliptical elements described above.
#' @seealso \code{\link{pca}}, \code{\pkg{ngft}}, \code{\link{rayleigh.filter}}
#' \itemize{
#' \item \href{https://en.wikipedia.org/wiki/S_transform}{S-Transform (Wikipedia)}
#' \item \href{http://gji.oxfordjournals.org/content/165/2/596.abstract}{Pinnegar (2006)}
#' Polarization analysis and polarization filtering of three-component signals with
#' the time-frequency S transform.
#' \item \href{http://www.bssaonline.org/content/105/1/210.abstract}{Meza-Fajardo et al (2015)}
#' Identification and Extraction of Surface Waves from Three-Component Seismograms
#' Based on the Normalized Inner Product.
#' }
#' @keywords ts

superlips <- function(X, Y, Z, nr=NA, nc=NA)
{
  if ( missing(X) || missing(Y) || missing(Z) )
    stop("Must provide input X, Y, and Z")
  if ( ! is.complex(X) || ! is.complex(Y) || ! is.complex(Z) )
    stop("X, Y, and Z must be complex")
  dimX <- dim(X)
  lx <- length(X)
  ly <- length(Y)
  lz <- length(Z)
  if ( lx != ly || lx != lz )
    stop("Length of X, Y, and Z must be the same")

  A <- Mod(X)^2 + Mod(Y)^2 + Mod(Z)^2
  B <- Re(X)^2 - Im(X)^2 + Re(Y)^2 - Im(Y)^2 + Re(Z)^2 - Im(Z)^2
  C <- -2 * (Re(X) * Im(X) + Re(Y) * Im(Y) + Re(Z) * Im(Z))

  a <- sqrt(0.5 * (A + sqrt(B^2 + C^2))) # semi-major axis
  b <- A - sqrt(B^2 + C^2)
  b[which(b < 0)] <- 0
  b <- sqrt(0.5 * b) # semi-minor axis

  # inclination (dip)
  I <- atan2(sqrt((Re(Z) * Im(Y) - Im(Z) * Re(Y))^2 + (Re(Z) * Im(X) - Im(Z) * Re(X))^2),
             Re(Y) * Im(X) - Im(Y) * Re(X))

  # azimuth of ascending node (strike)
  big.omega <- atan2(Re(Z) * Im(Y) - Im(Z) * Re(Y), Re(Z) * Im(X) - Im(Z) * Re(X))

  phi0 <- 0.5 * atan2(C, B)

  # angle between the ascending node and the long axis
  little.omega0 <- atan2(b * (Re(Z) * cos(phi0) - Im(Z) * sin(phi0)),
                         -a * (Re(Z) * sin(phi0) - Im(Z) * cos(phi0)))
  little.omega <- little.omega0 - 0.5 * pi * (sign(little.omega0) - 1)

  # phase, measured w.r.t. the time of maximum displacement
  phi <- phi0 + 0.5 * pi * (sign(little.omega0) - 1) * sign(phi0)

  # set output dimensions, filling by row
  if ( is.null(dimX) && ! is.na(nr) && ! is.na(nc)) {
    dimOut <- c(nc,nr)
    dim(a) <- dimOut
    a <- tr(a)
    dim(b) <- dimOut
    b <- tr(b)
    dim(I) <- dimOut
    I <- tr(I)
    dim(big.omega) <- dimOut
    big.omega <- tr(big.omega)
    dim(little.omega) <- dimOut
    little.omega <- tr(little.omega)
    dim(phi) <- dimOut
    phi <- tr(phi)
  }

  #----- return list of values
  ret = list(a=a, b=b, I=I, big.omega=big.omega, little.omega=little.omega, phi=phi)

  return(ret)
}