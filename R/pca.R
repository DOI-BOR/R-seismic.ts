#' Polarization attributes from principal component analysis
#'
#' \code{pca} determines polarization attributes of an input 3-component
#' time series using principal component analysis (PCA).
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
#' @param w.len If using the \code{analytic} option, average the complex covariance
#' matrix over \code{w.len} points. Default is 5 points.
#' @param analytic Use the analytic signal of each component instead of the
#' component itself. Default is \code{FALSE}.
#'
#' @details Computes polarization attributes from principal component analysis (PCA),
#' following the references listed below. Notes: (1) if \code{analytic = FALSE},
#' then average polarization attributes are computed for the entire window; (2) if
#' \code{analytic = TRUE}, then polarization attributes are computed for each point
#' in the window by forming the analytic signal vector of the input data (i.e., the
#' estimates are instantaneous); (3) ellipticity can only be determined if
#' \code{analytic = TRUE}; (4) if using the analytic signal method, then
#' \code{w.len >= 3} to determine the three non-zero eigenvalues needed for
#' estimating planarity and rectilinearity; (5) If \code{w.len <= 1}, then only a
#' single non-zero eigenvalue can be determined; (6) planarity and rectilinearity
#' estimates (\code{rect}, \code{plan}, and \code{dop}) are strongly dependent on
#' \code{w.len}, with values close to unity for small values of \code{w.len}; (7)
#' Using \code{analytic = TRUE} and \code{w.len = 5} will provide nearly identical
#' results to \code{\link{ipa}}.
#'
#' @return List of the polarization attributes and other data:
#' \code{pa} is a list with the following attributes:
#' \code{ire} - instantaneous reciprical ellipticity,
#' \code{az} - azimuth for principal eigenvector,
#' \code{plunge} - plunge angle of the principal eigenvector,
#' \code{rect} - rectilinearity,
#' \code{plan} - planarity,
#' \code{dop} - degree of planarity.
#' \code{ev} - return value(s) from \code{eigen}.
#' \code{alpha} - phase rotations to maximize the real part of the principal
#' eigenvector (see Vidale, 1986; Only provided if \code{analytic = TRUE}).
#' If \code{analytic = TRUE}, then \code{pa} and \code{alpha} will be returned
#' as \code{ts} objects.
#'
#' @seealso \code{\link{hanning}}, \code{\link{analytic.ts}}, \code{\link{eigen}},
#' \code{\link{ipa}}, \code{\link{super.ellips}}, \code{\link{rayleigh.filter}},
#' \code{\link{nip.filter}}, \code{\link{rayleigh.dir}}
#' \itemize{
#' \item \href{https://en.wikipedia.org/wiki/Analytic_signal}{Analytic Signal (Wikipedia)}
#' \item \href{http://www.bssaonline.org/content/76/5/1393.abstract}{Samson and Olson (1980)}
#' Some comments on the descriptions of the polarization states of waves.
#' \item \href{http://www.bssaonline.org/content/76/5/1393.abstract}{Vidale (1986)}
#' Complex polarization analysis of particle motion.
#' \item \href{http://bssa.geoscienceworld.org/content/78/5/1725.abstract}{Jurkevics (1988)}
#' Polarization analysis of three-component array data.
#' \item \href{http://www.bssaonline.org/content/81/2/622.abstract}{Bataille and Chiu (1991)}
#' Polarization analysis of high-frequency, three-component seismic data.
#' \item \href{http://www.bssaonline.org/content/94/3/1016.abstract}{Schimmel and Gallart (2004)}
#' Degree of Polarization Filter for Frequency-Dependent Signal Enhancement Through Noise Suppression.
#' \item \href{http://geophysics.geoscienceworld.org/content/71/5/V99.abstract}{Diallo et al (2006)}
#' Instantaneous polarization attributes based on an adaptive approximate covariance method.
#' \item \href{http://www.bssaonline.org/content/102/2/854.abstract}{Amoroso et al (2012)}
#' S-Wave Identification by Polarization Filtering and Waveform Coherence Analyses.
#' \item \href{http://gji.oxfordjournals.org/content/204/2/968.abstract}{Jones (2016)}
#' Quantifying the similarity of seismic polarizations.
#' }
#' @keywords ts
#' @export pca
#' @export
pca <- function(xt, yt=NA, zt=NA, dt=NA, demean=TRUE, pct=NA,
                w.len=5, analytic=FALSE)
{
  # get mts representing input data
	xyz <- get.xyz(xt, yt, zt, dt, demean, pct)
	len <- dim(xyz)[1]
	dt <- deltat(xyz)
	start <- start(xyz)

	if ( analytic ) {
	  # use analytic signal method for each point in time, if requested
	  if ( is.complex(xyz) )
	    stop("pca: input is already complex-valued")
	  xyz <- analytic.ts(xyz)

	 	# get the auto and cross covariances. Note: covariance matrix is Hermitian
  	xx <- xyz[,1] * Conj(xyz[,1])
  	xy <- xyz[,1] * Conj(xyz[,2])
  	xz <- xyz[,1] * Conj(xyz[,3])
  	yx <- xyz[,2] * Conj(xyz[,1])
  	yy <- xyz[,2] * Conj(xyz[,2])
  	yz <- xyz[,2] * Conj(xyz[,3])
  	zx <- xyz[,3] * Conj(xyz[,1])
  	zy <- xyz[,3] * Conj(xyz[,2])
  	zz <- xyz[,3] * Conj(xyz[,3])

  	# average the covariance matrix over multiple points, if requested.
  	# Note: if w.len=0, then only a single non-zero eigenvalue will be
  	# returned. Need w.len >= 3 to get non-zero eigenvalues.
  	if ( ! is.na(w.len) && w.len > 1 ) {
  	  xx <- mavg(xx,w.len)
  	  xy <- mavg(xy,w.len)
  	  xz <- mavg(xz,w.len)
  	  yx <- mavg(yx,w.len)
  	  yy <- mavg(yy,w.len)
  	  yz <- mavg(yz,w.len)
  	  zx <- mavg(zx,w.len)
  	  zy <- mavg(zy,w.len)
  	  zz <- mavg(zz,w.len)
  	}

  	# get cross-correlations at lag 0: rho.xy, rho.xz, and rho.yz
  	#rho <- sqrt(c((xy*xy)/(xx*yy),
  	#              (xz*xz)/(xx*zz),
  	#              (yz*yz)/(yy*zz)))

  	# get eigenvalue decomposition of the covariance matrix at each point in time
  	# note: eigen takes real or complex inputs
  	ev <- lapply(1:len, function(t) eigen(matrix(c(xx[t],xy[t],xz[t],
  	                                               yx[t],yy[t],yz[t],
  	                                               zx[t],zy[t],zz[t]),
  	                                             nrow=3,ncol=3,byrow=TRUE),
  	                                      symmetric=TRUE))
  	# ev[[t]]$values # 1x3 vector of eigenvalues, from largest to smallest
  	# ev[[t]]$vectors # 3x3 matrix whose columns are the normalized eigenvectors
  	# (in the same order as the eigenvalues)

  	# get angle needed to rotate eigenvector with largest eigenvalue so as to
  	# maximize its real part, following Vidale (1986)
  	alpha <- sapply(1:len, function(t) max.re(ev[[t]]$vectors[1,1],
  	                                          ev[[t]]$vectors[2,1],
  	                                          ev[[t]]$vectors[3,1]))

  	# rotate principal eignevector by alpha
  	for ( t in 1:len ) {
  	  if ( ! is.finite(alpha[t]) )
  	    warning("bad alpha for t=",t,": ",alpha[t])
      ev[[t]]$vectors[,1] <- ev[[t]]$vectors[,1] * exp(1i*alpha[t])
  	}

  	# get polarization attributes at each point in time, and put into an mts
    pa.values <- sapply(1:len, function(t) unlist(pca.attributes(ev[[t]]$vectors, ev[[t]]$values)))
    dim(pa.values) <- c(6,len)
    pa.values <- t(pa.values)
    colnames(pa.values) <- c("ire", "az", "plunge", "rect", "plan", "dop")

    ret <- list(pa=ts(pa.values, start=start, deltat=dt), ev=ev, alpha=alpha)
	} else {
	  # use the auto and cross covariances at lag 0 of the entire data series. Note:
	  # acf accepts only real inputs
	  acf.xyz <- acf(xyz, lag.max = 10, type="covariance", demean=demean, plot=FALSE)

  	# get cross-correlations at lag 0: rho[1,2], rho[1,3], and rho[2,3]
  	rho <- c((acf.xyz$acf[1,1,2]*acf.xyz$acf[1,1,2])/(acf.xyz$acf[1,1,1]*acf.xyz$acf[1,2,2]),
  					 (acf.xyz$acf[1,1,3]*acf.xyz$acf[1,1,3])/(acf.xyz$acf[1,1,1]*acf.xyz$acf[1,3,3]),
  					 (acf.xyz$acf[1,2,3]*acf.xyz$acf[1,2,3])/(acf.xyz$acf[1,2,2]*acf.xyz$acf[1,3,3]))
  	rho <- sqrt(rho)

  	# get eigenvalue decomposition of the covariance at lag 0
  	# note: eigen takes real or complex inputs
  	ev <- eigen(acf.xyz$acf[1,,], symmetric=TRUE)
  	# ev$values # 1x3 vector of eigenvalues, from largest to smallest
  	# ev$vectors # 3x3 matrix whose columns are the normalized eigenvectors (same order)

  	ret <- list(pa=pca.attributes(ev$vectors, ev$values), ev=ev)
	}

	return(ret)
}

# extract common polarization attributes from PCA eigenvectors and eigenvalues
#' @export pca.attributes
#' @export
pca.attributes <- function(vectors, values) {
  if ( nrow(vectors) != 3 && ncol(vectors) != 3 )
    error("pca.attributes: vectors must be a 3 x 3 matrix")
  if ( length(values) != 3 )
    error("pca.attributes: values must be a vector of length 3")
  if ( is.complex(values) )
    error("pca.attributes: eigenvalues should be real")

  # just need the real parts of the eigenvectors
  if ( is.complex(vectors) )
    vectors <- Re(vectors)

	#----- get azimuth for the principal eigenvector, measured in the positive
	# direction from X to Y, and defined by tan(az) = y/x
	# s3 = sign(vectors[3,1]) # Jurkevics (1988)
	# az <- atan2(vectors[2,1] * s3, vectors[1,1] * s3)
	az <- atan2(vectors[2,1], vectors[1,1])
	az <- az + 0.5 * pi * (1 - sign(cos(az))) # direction is ambiguous
	az <- normalize.angle(az) * rad2deg

	#----- get plunge angle for the principal eigenvector, measured from the horizontal,
	# and defined by tan(plunge) = (z/sqrt(x^2+y^2))
	# Vidale (1986)
	plunge <- atan2(vectors[3,1], sqrt(vectors[1,1]^2 + vectors[2,1]^2))
	# Jurkevics (1988)
	# plunge <- acos(vectors[3,1])

	# normalize so that 0 < plunge <= pi/2
	plunge <- plunge + 0.5 * pi * (1 - sign(cos(plunge))) # direction is ambiguous
	plunge <- abs(normalize.angle(plunge)) * rad2deg

	#------ get rectilinearity measure from eigenvalues
	# Jones et al (2016), pg 969, citing Jurkevics (1988); Also Vidale (1986), eqn. 8
	rect <- 1 - (values[2] + values[3]) / values[1]
	# Jurkevics (1988), pg. 1728
	# rect <- 1 - 0.5 * (values[2] + values[3]) / values[1]

	#------- degreee of polarization
	# Samson & Olson (1980) eqn. 18, Bataille & Chiu (1991) eqn. 8,
	# Schimmel & Gallart (2004) eqn, 5, Amoroso et al (2012) eqn. 3.
	# Note to self: the 0.5 belongs inside the sqrt, e.g., as written out
	# in B&C (1991), and follows from the double-summation (6 non-zero terms)
	dop <- sqrt(0.5 * ((values[1] - values[2])^2 +
	                      (values[1] - values[3])^2 +
	                      (values[2] - values[3])^2)) /
            abs(values[1] + values[2] + values[3])

	#----- get planarity measure from eigenvalues

	# Jurkevics (1988), pg. 1728
	plan <- 1 - 2 * values[3] / (values[1] + values[2])

	# Vidale (1986)
	# plan <- 1 - values[2] / values[3]

	#----- get elliptical component of polarization (analgous to b/a = IRE)
	# Vidale (1986)
  X <- sqrt(vectors[1,1]^2 + vectors[2,1]^2 + vectors[3,1]^2)
  ire <- sqrt(1 - X^2) / X

  #----- return list of values
	ret = list(ire=ire, az=az, plunge=plunge, rect=rect, plan=plan, dop=dop)

	return(ret)
}


# Find the phase that maximizes the real part of a
# complex vector, to implement Vidale (1986) method
#' @export max.re
#' @export
max.re <- function(x,y,z) {
  if ( ! is.complex(x) || ! is.complex(y) || ! is.complex(z) )
    error("max.re: inputes must be complex")

  # minimize -(X^2), where X is the real part of the vector
  x2 <- function(alpha) -(Re(x*exp(1i*alpha))^2 +
                            Re(y*exp(1i*alpha))^2 + Re(z*exp(1i*alpha))^2)

  #opt.val <- optim(0.5 * pi, x2, method="L-BFGS-B", lower=0, upper=pi)
  opt.val <- optim(0.5 * pi, x2, method="Brent", lower=0, upper=pi)

  if ( opt.val$convergence != 0 )
    warning("max.re: did not converge")

  alpha.max <- unlist(opt.val$par)

  return(alpha.max)
}


# circular moving average (for when FFT convolution with a boxcar is just too fast)
#' @export mavg
#' @export
mavg <- function(xt, w.len=1) {
  if ( ! is.numeric(xt) && ! is.complex(xt) )
    stop("mavg: xt must be numeric or complex")
  xt <- as.vector(xt)

  lx <- length(xt)
  if ( lx < 2 || w.len < 2 )
    return(xt)

  # create a w.len-by-lx matrix whose columns are the indices into xt for each window
  w.len <- min(w.len,lx)
  xt.index0 <- 1:w.len
  xt.index <- sapply(0:(lx-1), function(ix) (((xt.index0 - 1) + ix) %% lx) + 1)

  # for each window, sum xt over the window indices, and divide by w.len to get average
  xt.avg <- apply(xt.index, 2, function(ix) sum(xt[ix])) * (1. / as.numeric(w.len))

  # center the average by right-shifting by a half-window
  if ( w.len %% 2 == 0 )
    rshift <- w.len / 2 - 1
  else
    rshift <- (w.len + 1) / 2 - 1
  xt.ind <- c((lx-rshift+1):lx, 1:(lx-rshift))
  xt.avg <- xt.avg[xt.ind]

  return(xt.avg)
}


# make an mts from inputs
#' @export get.xyz
#' @export
get.xyz <- function(xt, yt=NA, zt=NA, dt=NA, demean=TRUE, pct=NA, w.len=NA)
{
  if ( ! is.na(pct) && (pct < 0 || pct > 50) )
    stop(sprintf("pct (%.2f) must be between 0 and 50",pct))

  # split out the three components
  multi.trace <- is.mts(xt) || is.matrix(xt) || length(dim(xt)) > 1 ||
    ( is(xt, "signalSeries") && ! is.null(dim(xt)) )
  if ( multi.trace ) {
    if ( dim(xt)[2] < 3 )
      stop("multivariate time series must have at least 3 components (x, Y, and z)")
    len <- dim(xt)[1]
    if ( is(xt, "signalSeries") ) {
      x <- xt@data[,1]
      y <- xt@data[,2]
      z <- xt@data[,3]
    } else {
      x <- xt[,1]
      y <- xt[,2]
      z <- xt[,3]
    }
    if ( is(xt, "signalSeries") || is(xt, "ts") )
      dt <- deltat(xt)
  } else {
    if ( is.na(yt) || is.na(zt) )
      stop("Missing yt and/or zt inputs")
    x <- as.vector(xt)
    y <- as.vector(yt)
    z <- as.vector(zt)
    # trim vectors to minimum length
    len <- min(length(x), length(y), length(z))
    x <- x[1:len]
    y <- y[1:len]
    z <- z[1:len]
  }

  # de-mean, if requested
  if ( demean ) {
    x <- x - mean(x)
    y <- y - mean(y)
    z <- z - mean(z)
  }

  # smooth input data, if requested
  if ( ! is.na(w.len) && w.len > 1 ) {
    x <- mavg(x,w.len)
    y <- mavg(y,w.len)
    z <- mavg(z,w.len)
  }

  start <- 0
  if ( is.ts(xt) || is.mts(xt) )
    start <- start(xt)[1]
  else if ( is(xt, "signalSeries") )
    start <- xt@positions@from

  if ( is.na(dt) )
    dt <- 0.01

  # bind components into an mts
  xyz <- ts(cbind(x,y,z), start=start, deltat=dt)

  # taper, if requested
  if ( ! is.na(pct) )
    xyz <- hanning(xyz, pct=pct, demean=demean)

  return(xyz)
}
