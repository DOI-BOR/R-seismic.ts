#' Rayleigh and Love wave filtering
#'
#' \code{nip.filter} generates a multiplicative filter from the Fourier or
#' S-transform of a 3-component seismogram, using the normalized inner product (NIP)
#' method.
#'
#' @param X,Y,Z Equal-length \code{complex} input vectors or matrices for the X, Y,
#' and Z directions, and representing the complex Fourier or S Transforms of the
#' original 3-component data. The layout of the input can vary with the type of
#' the transform, but this function only needs to be able to ascertain which elements
#' correspond to positive, negative, and zero frequencies. If the input is a
#' matrix, then it is assumed to be in "standard" form with \code{nf} frequency rows
#' and \code{nt} time columns, filled by row (e.g., \code{dim(X) = c(nf,nt)}), and
#' no further parameters are needed to specify the layout. If the input is a vector,
#' then additional parameters are required to specify the layout: (1) for input vectors
#' representing a standard S-transform matrix in vector form,
#' \code{X[1:N] = C(X[f0,1:nt],X[f1,1:nt],...X[n/2;1:nt],X[-n/2,1:nt]...X[-1,1:nt])},
#' the \code{nf} parameter must be specified, and \code{nt} will be computed;
#' (2) for input vectors representing a variable time and frequency partitioning
#' (e.g., a fast S-transform), the layout must have zero-frequency values first,
#' followed by the positive and then negative frequency values, and the
#' \code{posf.start} and \code{negf.start} parameters must be used to specify the
#' starting indices of the positive and negative frequency values.
#' @param nf Number of frequencies (rows). Must be specified for the case of Fourier
#' transforms, or for S-Transforms which are input as N-length column vectors (i.e,
#' for each time index \code{1:nt}, there are \code{1:nf} frequency values, and
#' \code{N = nf * nt}). Frequencies are assumed to be in standard order (zero,
#' positive, negative). Not used otherwise.
#' @param reject set to \code{TRUE} to reject Rayleigh waves, or \code{FALSE} to
#' reject everything but Rayleigh waves. Default is \code{TRUE}
#' @param motion Type of motion to filter. One of \code{c("prograde", "retrograde",
#' "all")}. Only the first character is used. Default is \code{"all"}.
#' @param posf.start,negf.start Index into \code{X}, \code{Y}, and \code{Z} of
#' the starting index for positive and negative frequencies. Only needed for
#' S-transforms that are NOT in standard (rectangular) form, such as fast S-transforms.
#'
#' @details The code implements methods based on Meza-Fajardo et al (2015) to
#' compute the filter. See Meza-Fajardo et al (2015) for details. Either the \code{dim}
#' attribute must be set on the \code{ X, Y, Z} inputs, or the layout must be
#' specified using the \code{nf} or \code{posf.start, negf.start} paramaters.
#' @return List that includes the filter and other values. The layout of these values will
#' be a row-filled matrix for matrix inputs, or a vector for variable frequency-time
#' partitions. The elements of the returned list are as follows:
#' \describe{
#' \item{F}{Multiplicative filter in the transform domain. This filter can by used to multiply select
#' or reject data in the transform domain. Inverse transform the data after multiplication
#' by the filter to get the time-domain result.}
#' \item{nip}{Normalized inner-product. Will range from -1 to 1.}
#' \item{R}{Radial values. The radial direction depends on frequency and time}
#' \item{T}{Transverse values. The transverse direction depends on frequency and time}
#' \item{ZH}{Hilbert transform of vertical values. \code{ZH} is correlated with
#' \code{R}} to get the radial direction.
#' \item{theta}{The complement of \code{phi}: \code{theta = pi/2 - phi}}
#' \item{phi}{The angle to the radial direction, measured from the positive \code{X}
#' axis towards the positive \code{Y} asis, about the positive \code{Z} axis.}
#' \item{rho}{another estimate of theta. The radial direction depends on frequency and time}
#' }
#' @seealso \code{\pkg{ngft}}, \code{\link{rayleigh.filter}}
#' \itemize{
#' \item \href{https://en.wikipedia.org/wiki/S_transform}{S-Transform (Wikipedia)}
#' \item \href{http://www.bssaonline.org/content/105/1/210.abstract}{Meza-Fajardo et al (2015)}
#' Identification and Extraction of Surface Waves from Three-Component Seismograms
#' Based on the Normalized Inner Product.
#' }
#' @keywords ts

nip.filter <- function(X, Y, Z, nf=NA, reject=TRUE, motion="all",
                       posf.start=NA, negf.start=NA) {
  if ( missing(X) || missing(Y) || missing(Z) )
    stop("Must provide input X, Y, and Z")
  if ( ! is.complex(X) || ! is.complex(Y) || ! is.complex(Z) )
    stop("X, Y, and Z must be complex")
  lx <- length(X)
  ly <- length(Y)
  lz <- length(Z)
  if ( lx != ly || lx != lz )
    stop("Length of X, Y, and Z must be the same")

  # Approach: get indices defining where positive and negative frequencies start
  # in the input Fourier or S-transforms.
  # This code needs to handle S-transforms in several forms. As long as the input
  # transforms are arranged as consecutive indices of: (1) zero-frequency elements;
  # (2) positive-frequency elements; and, (3) negative-frequency elements, then
  # this code won't need to know any other details.
  matrix.in <- FALSE
  if ( is.na(posf.start) || is.na(negf.start) ) {
    # handle standard rectangular S-transform input stored as a column-filled
    # vector, or as a row-filled matrix
    nt <- NA
    dimx <- dim(X)
    if ( is.null(dimx) ) {
      # assume inputs are vectors in the form C(X[f0,1:nt],X[f1,1:nt],...
      # X[n/2;1:nt],X[-n/2,1:nt]...X[-1,1:nt]). Compute nt = number
      # of time-steps from the input length and number of frequencies
      if ( is.na(nf) )
        stop("Must specify nf")
      if ( lx %% nf != 0 )
        stop("Input vector length is not a multiplier of nf")
      nt <- lx / nf
    } else {
      # assume inputs are row-filled matrices with nf rows and nt columns
      nf <- dimx[1]
      nt <- dimx[2]
      # convert matrix to vectors in the form C(X[fo,1:nt],X[f1,1:nt],...
      # X[n/2;1:nt],X[-n/2,1:nt]...X[-1,1:nt]).
      X <- as.vector(t(X))
      Y <- as.vector(t(Y))
      Z <- as.vector(t(Z))
    }
    matrix.in <- TRUE
    posf.start <- nt + 1
    if ( nf %% 2 == 0 ) {
      # number of frequencies is even
      negf.start <- (nf/2 + 1) * nt + 1
    } else {
      # number of frequencies is odd
      negf.start <- ((nf + 1)/2) * nt + 1
    }
  }
  f0.end <- posf.start - 1
  posf.end <- negf.start - 1
  negf.end <- lx

  #warning("fo.end=",f0.end," posf.start=",posf.start," posf.end=",posf.end,
  #  " negf.start=",negf.start," negf.end=",negf.end)

  # This function uses the Hilbert transform of the Z component. In frequency space,
  # this requires shifting positive frequencies by -pi/2 and negative frequencies
  # by +pi/2. To apply the phase shift, we need to multiply each element of
  # the transform by exp(-/+ i pi/2) = -/+ i, for +/- frequencies
  ZH <- complex(lx)
  ZH[1:f0.end] <- Z[1:f0.end] # for now, just copy 0-frequency elements
  ZH[posf.start:posf.end] <- Z[posf.start:posf.end] * -1i
  ZH[negf.start:negf.end] <- Z[negf.start:negf.end] * 1i

  # get the Rayleigh wave azimuth at each frequency and time
  phi <- atan2(Re(Y)*Re(ZH) + Im(Y)*Im(ZH), Re(X)*Re(ZH) + Im(X)*Im(ZH))
  phi[1:f0.end] = atan2(Re(Y[1:f0.end]), Re(X[1:f0.end])) # handle 0-frequency case
  if ( ! startsWith(motion,"p") ) {
    phi <- phi - pi
    ind <- which(phi < -pi)
    phi[ind] <- phi[ind] + 2 * pi
    ind <- which(phi > pi)
    phi[ind] <- phi[ind] - 2 * pi
  }

  # get radial and transverse components at each frequency and time
  R <- X * cos(phi) + Y * sin(phi)
  T <- -X * sin(phi) + Y * cos(phi)

  # get NIP between R and Hilbert transform of Z
  NIP.RZH <- (Re(R) * Re(ZH) + Im(R) * Im(ZH)) / (Mod(R) * Mod(ZH))

  # Construct filter from NIP. Rayleigh waves should have |NIP| > 0.8
  v1 = 0.9
  v2 = 0.8
  nip.selector <- NIP.RZH
  if ( startsWith(motion,"a") ) {
    nip.selector <- abs(NIP.RZH)
  } else if ( startsWith(motion,"r") ) {
    v1 <- -v1
    v2 <- -v2
  }
  F <- cos.filter(nip.selector, v1, v2)
  if ( reject )
    F <- 1 - F

  theta <- 0.5 * pi - phi
  ind <- which(theta < -pi)
  theta[ind] <- theta[ind] + 2 * pi
  ind <- which(theta > pi)
  theta[ind] <- theta[ind] - 2 * pi

  # get reference direction transverse to theta
  rho.r <- atan2(Re(Y), Re(X))

  if ( matrix.in ) {
    dimtf <- c(nt, nf)
    dim(F) <- dimtf
    F <- t(F)
    dim(NIP.RZH) <- dimtf
    NIP.RZH <- t(NIP.RZH)
    dim(R) <- dimtf
    R <- t(R)
    dim(T) <- dimtf
    T <- t(T)
    dim(ZH) <- dimtf
    ZH <- t(ZH)
    dim(theta) <- dimtf
    theta <- t(theta)
    dim(phi) <- dimtf
    phi <- t(phi)
    dim(rho.r) <- dimtf
    rho.r <- t(rho.r)
  }

  #----- return list of values
  ret = list(F=F, nip=NIP.RZH, R=R, T=T, ZH=ZH, theta=theta,
             phi=phi, rho=rho.r)

  return(ret)
}