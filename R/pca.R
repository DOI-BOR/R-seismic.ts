#' Polarization analysis using principal components
#'
#' \code{pca} determines several polarization measures of an input 3-component
#' time series using principal component analysis
#'
#' @param xt,yt,zt Equally-sampled univariate input series for the the X, Y,
#' and Z directions. Alternatively, \code{xt} can be a multivariate time series
#' with X, Y, and Z data in the first 3 positions (and in that order). Must
#' convert to numeric \code{\link{vector}}, \code{\link{matrix}},
#' \code{\link{data.frame}}, \code{\link{ts}}, or \code{\link{signalSeries}}.
#' @param dt Sample interval, in seconds. Default is 0.01, or \code{deltat}
#' if \code{xt} is a \code{\link{ts}} or \code{\link{signalSeries}} object.
#' @param demean TRUE if the windowed data should be demeaned? Default is FALSE.
#' @param pct Percentage of data window to apply a \code{\link{hanning}} taper.
#' Must be between 0 and 50. Default is 0 (no taper).
#'
#' @details Here's how it works...
#' @return List with the polarization directions and measures.
#' @seealso \code{\link{hanning}}, \code{\link{analytic.ts}},
#' \code{\link{acf.xyz}}, \code{\link{eigen}}
#'
#' @keywords ts

pca <- function(xt, yt=NULL, zt=NULL, dt=NULL, demean=TRUE, pct=NA)
{
	# pca(ltoe.data$ltoe.l.disp[1850:2365],ltoe.data$ltoe.t.disp[1850:2365],ltoe.data$ltoe.z.disp[1850:2365],pct=20)
	if ( ! is.na(pct) && (pct < 0 || pct > 50) )
		stop(sprintf("pct (%.2f) must be between 0 and 50",pct))

  multi.trace <- is.mts(xt) || is.matrix(xt) || length(dim(xt)) > 1 ||
    ( is(xt, "signalSeries") && ! is.null(dim(xt)) )

  if ( multi.trace == TRUE ) {
    if ( dim(xt)[2] < 3 )
      stop("multivariate time series must have at least 3 components (x, Y, and z)")
    len <- dim(xt)[1]
    if ( is(xt, "signalSeries") ) {
      x <- xt[,1]@data
      y <- yt[,2]@data
      z <- zt[,3]@data
    } else {
      if ( is.na(yt) || is.na(zt) )
        stop("Missing yt and/or zt inputs")
      x <- xt[,1]
      y <- yt[,2]
      z <- zt[,3]
    }
    if ( is(xt, "signalSeries") || is(xt, "ts") )
      dt <- deltat(xt)
  } else {
  	x <- as.vector(x.data)
  	y <- as.vector(y.data)
  	z <- as.vector(z.data)
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
		xyz <- hanning(xyz, pct = pct, demean = TRUE)

	# convert to analytic time series, if requested
	#if ( analytic )
	#	xyz <- analytic.ts(xyz)

	# get the auto and cross covariances. Note: acf only accepts real inputs
	acf.xyz <- acf(xyz, lag.max = 10, type="covariance", demean=TRUE, plot=FALSE)
	rho <- c((acf.xyz$acf[1,1,2]*acf.xyz$acf[1,1,2])/(acf.xyz$acf[1,1,1]*acf.xyz$acf[1,2,2]),
					 (acf.xyz$acf[1,1,3]*acf.xyz$acf[1,1,3])/(acf.xyz$acf[1,1,1]*acf.xyz$acf[1,3,3]),
					 (acf.xyz$acf[1,2,3]*acf.xyz$acf[1,2,3])/(acf.xyz$acf[1,2,2]*acf.xyz$acf[1,2,2]))
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

	# get azimuth for principal eigenvector, defined by tan(az) = y/x
	s3 = sign(ev$vectors[3,1]) # Jurkevics (1988)
	az <- atan2(ev$vectors[2,1] * s3, ev$vectors[1,1] * s3)
	az <- az * 180 / pi

	# get angle of incidence for principal eigenvector, defined by tan(ain) = (sqrt(x^2+y^2)/z)
	ain <- atan2(sqrt(ev$vectors[2,1] * ev$vectors[2,1] + ev$vectors[1,1] * ev$vectors[1,1]), ev$vectors[3,1])
	# ain <- acos(ev$vectors[3,1]) # Jurkevics (1988)
	ain <- ain * 180 / pi

	# get rectilinearity measure from eigenvalues
	rect <- 1 - (ev$values[2] + ev$values[3]) / ev$values[1]
	# rect <- 1 - 0.5 * (ev$values[2] + ev$values[3]) / ev$values[1] # Jurkevics (1988)
	# rect <- sqrt(0.5 * ((ev$values[1] - ev$values[2])^2 + (ev$values[1] - ev$values[3])^2 + (ev$values[2] - ev$values[3])^2)) / (ev$values[1] + ev$values[2] + ev$values[3]) # Samson & Olson (1980), Bataille & Chiu (1991)

	# get planarity measure from eigenvalues
	plan <- 1 - 2 * ev$values[3] / (ev$values[1] + ev$values[2]) # Jurkevics (1988)
	# plan <- 1 - ev$values[2] / ev$values[3] # Vidale (1986)

	ret = list(rho=rho, az=az, ain=ain, rect=rect, plan=plan, acf=acf.xyz$acf[1,,], ev=ev)

	return(ret)
}