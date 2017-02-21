#' Polarization analysis using principal components
#'
#' \code{pca} determines several polarization measures of an input 3-component
#' time series using principal component analysis.
#'
#' @param xt,yt,zt Equally-sampled univariate input series for the X, Y,
#' and Z directions. Alternatively, \code{xt} can be a multivariate time series
#' with X, Y, and Z data in the first 3 positions (and in that order). Input must
#' convert to numeric \code{\link{vector}}, \code{\link{matrix}},
#' \code{\link{data.frame}}, \code{\link{ts}}, or \code{\link{signalSeries}}.
#' @param dt Sample interval, in seconds. Default is 0.01, or \code{deltat}
#' if \code{xt} is a \code{\link{ts}} or \code{\link{signalSeries}} object.
#' @param demean TRUE if the windowed data should be demeaned? Default is FALSE.
#' @param pct Percentage of data window to apply a \code{\link{hanning}} taper.
#' Must be between 0 and 50. Default is 0 (no taper).
#'
#' @details Computes polarization measures
#' @return List with the polarization directions and measures.
#' @seealso \code{\link{hanning}}, \code{\link{analytic.ts}},
#' \code{\link{acf.xyz}}, \code{\link{eigen}}
#'
#' @keywords ts

pca <- function(xt, yt=NA, zt=NA, dt=NA, demean=TRUE, pct=NA)
{
	# pca(ltoe.data$ltoe.l.disp[1850:2365],ltoe.data$ltoe.t.disp[1850:2365],ltoe.data$ltoe.z.disp[1850:2365],pct=20)
	if ( ! is.na(pct) && (pct < 0 || pct > 50) )
		stop(sprintf("pct (%.2f) must be between 0 and 50",pct))

  multi.trace <- is.mts(xt) || is.matrix(xt) || length(dim(xt)) > 1 ||
    ( is(xt, "signalSeries") && ! is.null(dim(xt)) )

  if ( multi.trace ) {
    if ( dim(xt)[2] < 3 )
      stop("multivariate time series must have at least 3 components (x, Y, and z)")
    len <- dim(xt)[1]
    if ( is(xt, "signalSeries") ) {
      x <- xt[,1]@data
      y <- xt[,2]@data
      z <- xt[,3]@data
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

  if ( is.na(dt) )
    dt <- 0.01

	# bind data into an mts
	xyz <- ts(cbind(x,y,z), deltat=dt)

	# taper, if requested
	if ( ! is.na(pct) )
		xyz <- hanning(xyz, pct=pct, demean=demean)

	# convert to analytic time series, if requested
	# TODO: acf only works with real inputs; need replacement to work with complex
	#if ( analytic )
	#  xyz <- analytic.ts(xyz)

	# get the auto and cross covariances. Note: acf accepts only real inputs
	acf.xyz <- acf(xyz, lag.max = 10, type="covariance", demean=demean, plot=FALSE)

	# get cross-correlations at lag 0: rho[1,2], rho[1,3], and rho[2,3]
	rho <- c((acf.xyz$acf[1,1,2]*acf.xyz$acf[1,1,2])/(acf.xyz$acf[1,1,1]*acf.xyz$acf[1,2,2]),
					 (acf.xyz$acf[1,1,3]*acf.xyz$acf[1,1,3])/(acf.xyz$acf[1,1,1]*acf.xyz$acf[1,3,3]),
					 (acf.xyz$acf[1,2,3]*acf.xyz$acf[1,2,3])/(acf.xyz$acf[1,2,2]*acf.xyz$acf[1,3,3]))
	rho <- sqrt(rho)

	# get SVD
	#s <- La.svd(acf.xyz$acf[10,,])
	#U <- s$u
	#D <- diag(s$d)
	#Vt <- s$vt

	# get eigenvalue decomposition of the covariance at lag 0
	# note: eigen takes real or complex inputs
	ev <- eigen(acf.xyz$acf[1,,], symmetric=TRUE)
	#ev$values # 1x3 vector of eigenvalues
	#ev$vectors # 3x3 matrix whose columns are the eigenvectors

	#----- get azimuth for principal eigenvector, defined by tan(az) = y/x,
	# recalling that 1 = X, 2 = Y, 3 = Z
	s3 = sign(ev$vectors[3,1]) # Jurkevics (1988)
	az <- atan2(ev$vectors[2,1] * s3, ev$vectors[1,1] * s3)
	az <- az * 180 / pi # convert to degrees

	# assume that X = E, Y = N, Z = up; az is ccl from X (E), but we
	# want azimuth cl from N, so need to adjust
	az <- 90 - az

	#----- get angle of incidence for the principal eigenvector,
	# defined by tan(ain) = (sqrt(x^2+y^2)/z)
	# Vidale (1986)
	ain <- atan2(sqrt(ev$vectors[2,1] * ev$vectors[2,1] + ev$vectors[1,1] * ev$vectors[1,1]), ev$vectors[3,1])
	# Jurkevics (1988)
	# ain <- acos(ev$vectors[3,1])
	ain <- ain * 180 / pi # convert to degrees

	# normalize so that 0 < ain <= 90, and 0 <= az < 360
	if ( ain > 90 ) {
	  az <- az + 180
	  ain <- 180 - ain
	}
	if ( ain <= 0 ) {
	  az <- az + 180
	  ain <- 180 + ain
	}
	while ( az > 360 )
	  az <- az - 360
	while ( az < 0 )
	  az <- az + 360

	#------ get rectilinearity measure from eigenvalues
	# Jones et al (2016), pg 969, citing Jurkevics (1988)
	rect <- 1 - (ev$values[2] + ev$values[3]) / ev$values[1]

	# Jurkevics (1988), pg. 1728
	# rect <- 1 - 0.5 * (ev$values[2] + ev$values[3]) / ev$values[1]

	#------- degreee of polarization
	# Samson & Olson (1980) eqn. 18, Bataille & Chiu (1991) eqn. 8,
	# Schimmel & Gallart (2004) eqn, 5, Amoroso et al (2012) eqn. 3.
	# Note to myself: the 0.5 belongs inside the sqrt, e.g., as written out
	# in B&C (1991), and follows from the double-summation (6 non-zero terms)
	dop <- sqrt(0.5 * ((ev$values[1] - ev$values[2])^2 +
	                      (ev$values[1] - ev$values[3])^2 +
	                      (ev$values[2] - ev$values[3])^2)) /
            abs(ev$values[1] + ev$values[2] + ev$values[3])

	#----- get planarity measure from eigenvalues

	# Jurkevics (1988), pg. 1728
	plan <- 1 - 2 * ev$values[3] / (ev$values[1] + ev$values[2])

	# Vidale (1986)
	# plan <- 1 - ev$values[2] / ev$values[3]

	#----- get ellipticity
	# Vidale (1986)


	#----- return list of values
	ret = list(rho=rho, az=az, ain=ain, rect=rect, plan=plan, dop=dop, acf=acf.xyz$acf[1,,], ev=ev)

	return(ret)
}