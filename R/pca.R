pca <-
function(x.data, y.data, z.data, dt=0.01, demean = TRUE, pct = NA)
{
	# pca(ltoe.data$ltoe.l.disp[1850:2365],ltoe.data$ltoe.t.disp[1850:2365],ltoe.data$ltoe.z.disp[1850:2365],pct=20)
	if ( ! is.na(pct) && (pct < 0 || pct > 50) )
		stop(sprintf("pct (%.2f) must be between 0 and 50",pct))

	# treat input as vector, and trim to minimum length
	x <- as.vector(x.data)
	y <- as.vector(y.data)
	z <- as.vector(z.data)
	len <- min(length(x), length(y), length(z))
	x <- x[1:len]
	y <- y[1:len]
	z <- z[1:len]

	# bind the input series into a multivariate time series
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