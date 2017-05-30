#' Rayleigh and Love wave filtering
#'
#' \code{nip.filter} make a filter from the Fourier or S-transform of a 3-component
#' seismogram using the normalized inner product method.
#'
#' @param X,Y,Z Equal-length \code{complex} vectors or matrices for the X, Y,
#' and Z directions representing the complex Fourier or S Transforms
#' of the original 3-component data. Each row is for a specific frequency
#' @param nf Number of frequencies (rows)
#' @param reject.rayleigh set to \code{TRUE} to reject Rayleigh waves, or \code{FALSE} to
#' reject everything but Rayleigh waves. Default is \code{TRUE}
#' @param motion One of \code{c("prograde", "retrograde", "all")}.
#' Only the first character is used. Default is \code{"all"}.
#'
#' @details The code implements the methods of Meza-Fajardo et al (2015) to compute
#' the filter. See Meza-Fajardo et al (2015) for details. Either the \code{dim}
#' attribute must be set on the inputs, or the \code{nf} argument must be set.
#' @return List that includes a filter \code{F} for rejecting or selecting Rayleigh
#' waves, as a matrix with \code{dim = c(nf,nt)}.
#' @seealso \code{\pkg{ngft}}, \code{\link{rayleigh.filter}}
#' \itemize{
#' \item \href{https://en.wikipedia.org/wiki/S_transform}{S-Transform (Wikipedia)}
#' \item \href{http://www.bssaonline.org/content/105/1/210.abstract}{Meza-Fajardo et al (2015)}
#' Identification and Extraction of Surface Waves from Three-Component Seismograms
#' Based on the Normalized Inner Product.
#' }
#' @keywords ts

nip.filter <- function(X, Y, Z, nf=NA, reject.rayleigh=TRUE, motion="all") {
  if ( missing(X) || missing(Y) || missing(Z) )
    stop("Must provide input X, Y, and Z")
  if ( ! is.complex(X) || ! is.complex(Y) || ! is.complex(Z) )
    stop("X, Y, and Z must be complex")
  lx <- length(X)
  ly <- length(Y)
  lz <- length(Z)
  if ( lx != ly || lx != lz )
    stop("Length of X, Y, and Z must be the same")
  nt <- NA
  dimx <- dim(X)
  if ( is.null(dimx) ) {
    # convert input vectors to matrices with nf rows and nt columns
    if ( is.na(nf) )
      stop("Must specify nf")
    if ( lx %% nf != 0 )
      stop("Input vector length is not a multiplier of nf")
    nt <- lx / nf
    dimxt <- c(nt, nf)
    dim(X) <- dimxt
    X <- t(X)
    dim(Y) <- dimxt
    Y <- t(Y)
    dim(Z) <- dimxt
    Z <- t(Z)
  } else {
    # assume already have input matrices with nf rows and nt columns
    nf <- dimx[1]
    nt <- dimx[2]
  }

  # Need to shift positive frequencies by -pi/2 and negative frequencies
  # by +pi/2. Since each row of the S-transform is the time history at that
  # row's frequency, to apply the phase shift, we just need to multiply each
  # element in the row by exp(-/+ i pi/2) = -/+ i. Positive frequenciy terms
  # are therefore multiplied by -i, and negative frequency terms are multiplied
  # by +i.
  # get FT of Hilbert Transform operator
  Hf <- complex(nf)
  Hf[1] <- 1 # Set 0-frequency value to 1 (rather than 0) for now to avoid NaN
  if ( nf %% 2 == 0 ) {
    # original number of frequencies is even
    Hf[2:(nf/2+1)] = -1i # positive frequencies
    Hf[(nf/2+2):nf] = 1i # negative frequencies
  } else {
    # original number of frequencies is odd
    Hf[2:((nf+1)/2)] = -1i # positive frequencies
    Hf[((nf+1)/2+1):nf] = 1i # negative frequencies
  }

  # Get Fourier Transform of Hilbert Transform of Z. Each column of Z
  # is multiplied by Hf
  ZH <- Z * Hf # FT(H(z))

  # theta is the angle measured from Y towards X
  theta.r <- atan2(Re(X)*Re(ZH) + Im(X)*Im(ZH), Re(Y)*Re(ZH) + Im(Y)*Im(ZH))
  theta.r[1,] = atan2(Re(X[1,]), Re(Y[1,])) # handle 0-frequency case
  theta.l <- theta.r + pi * (1 - sign(sin(theta.r))) +
    0.5 * pi * (1 - sign(cos(theta.r))) * sign(sin(theta.r))
  x.pr <- 1 # assume -X to +X propagation (-1 for +X to -X)
  theta <- theta.l + 0.5 * pi * (sign(sin(theta.l)) - sign(x.pr))

  # get radial and transverse components
  R <- Y * cos(theta) + X * sin(theta)
  T <- -Y * sin(theta) + X * cos(theta)

  # get NIPR and filter
  NIP.R <- (Re(R) * Re(ZH) + Im(R) * Im(ZH)) / (Mod(R) * Mod(ZH))
  # Rayleigh waves should have |NIP| > 0.8
  v1 = 0.9
  v2 = 0.8
  nip.selector <- NIP.R
  if ( startsWith(motion,"a") ) {
    nip.selector <- abs(NIP.R)
  } else if ( startsWith(motion,"r") ) {
    v1 <- -v1
    v2 <- -v2
  }
  F <- cos.filter(nip.selector, v1, v2)
  if ( reject.rayleigh )
    F <- 1 - F
  dim(F) <- dim(NIP.R)
  if ( ! reject.rayleigh )
    ZH <- ZH * F

  # get positive angle measured from x towards y
  phi <- 0.5 * pi - theta
  ind <- which(phi < -pi)
  phi[ind] <- phi[ind] + 2 * pi
  ind <- which(phi > pi)
  phi[ind] <- phi[ind] - 2 * pi

  phi.r <- 0.5 * pi - theta.r
  ind <- which(phi.r < -pi)
  phi.r[ind] <- phi.r[ind] + 2 * pi
  ind <- which(phi.r > pi)
  phi.r[ind] <- phi.r[ind] - 2 * pi

    # get reference direction transverse to theta
  rho.r <- atan2(Re(X), Re(Y))

  #----- return list of values
  ret = list(F=F, nip=NIP.R, R=R, T=T, ZH=ZH, theta=theta, theta.r=theta.r,
             phi=phi, phi.r=phi.r, rho.r=rho.r)

  return(ret)
}