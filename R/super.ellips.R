#' 3D waveform decomposition using superposition of ellipses
#'
#' \code{super.ellips} decomposes a complex transform of a 3-component
#' seismic waveform into an equivalent superposition of ellipses, from which
#' polarization attributes can be determined. It implements the method of
#' Pinnegar (2006).
#'
#' @param X,Y,Z \code{complex} transformed input data in the X, Y, and Z directions
#' representing one of the following types: (1) analytic signal; (2) Fourier transform;
#' or, (3) S transform. The inputs must be of the same type and length. Input class
#' can be \code{\link{ts}}, \code{\link{vector}}, or \code{\link{matrix}}.
#' @param nrow,ncol If specified, and \code{X, Y} and \code{Z} are of class
#' \code{\link{matrix}}, then convert the ouput elliptical attributes to matrix
#' form, with \code{\link{dim} = c(nrow,ncol)}, and filled by row. Otherwise, the
#' outputs will have the same class as the inputs.
#'
#' @details The code implements the methods of Pinnegar (2006). The elliptical
#' elements computed by this function can subsequently be used to compute
#' polarization attributes. For example, the instaneous reciprocal ellipticity
#' (IRE) is determined from \code{IRE = b/a}. Filters can be implemented using
#' the elliptical elements or derived polarization attributes to select or reject
#' certain wave types (e.g., see \code{\link{rayleigh.filter}}).
#' @return A \code{\link{list}} is returned with the following elliptical elements
#' for each input point:
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
#' Each element has the same class as the input, unless the \code{nrow} and
#' \code{ncol} arguments are provided, in which case the returned elements will
#' be of class \code{link{matrix}}.
#' @seealso \code{\link{rayleigh.filter}}, \code{\link{analytic.ts}}, \pkg{\link{ngft}},
#' \code{\link{pca}}, \code{\link{ipa}}
#'
#' \itemize{
#' \item \href{https://en.wikipedia.org/wiki/S_transform}{S-Transform (Wikipedia)}
#' \item \href{http://gji.oxfordjournals.org/content/165/2/596.abstract}{Pinnegar (2006)}
#' Polarization analysis and polarization filtering of three-component signals with
#' the time-frequency S transform.
#' }
#' @keywords ts

super.ellips <- function(X, Y, Z, nr=NA, nc=NA)
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

  # angle between the ascending node and the long axis (rake)
  little.omega0 <- atan2(b * (Re(Z) * cos(phi0) - Im(Z) * sin(phi0)),
                         -a * (Re(Z) * sin(phi0) + Im(Z) * cos(phi0)))
  little.omega <- little.omega0 - 0.5 * pi * (sign(little.omega0) - 1)

  # phase, measured w.r.t. the time of maximum displacement
  phi <- phi0 + 0.5 * pi * (sign(little.omega0) - 1) * sign(phi0)

  # plunge of the semi-major axis
  theta <- asin(sin(I) * sin(little.omega))

  # Get the azimuth of the semi-major axis.
  # Pinnegar (2006), which seems to be incorrect:
  #zeta <- atan2(tan(big.omega) + tan(little.omega) * tan(I),
  #              1 - tan(big.omega) * cos(I) * tan(big.omega))
  # CKW (2017), from inspection of the geometry:
  zeta <- big.omega + atan(tan(little.omega) * cos(I))

  # set output dimensions, filling by row
  if ( is.null(dimX) && ! is.na(nr) && ! is.na(nc)) {
    dimOut <- c(nc,nr)
    dim(a) <- dimOut
    a <- t(a)
    dim(b) <- dimOut
    b <- t(b)
    dim(I) <- dimOut
    I <- t(I)
    dim(big.omega) <- dimOut
    big.omega <- t(big.omega)
    dim(little.omega) <- dimOut
    little.omega <- t(little.omega)
    dim(theta) <- dimOut
    theta <- t(theta)
    dim(zeta) <- dimOut
    zeta <- t(zeta)
    dim(phi) <- dimOut
    phi <- t(phi)
  }

  #----- return list of values
  ret = list(a=a, b=b, I=I, big.omega=big.omega, little.omega=little.omega,
             phi=phi, theta=theta, zeta=zeta)

  return(ret)
}