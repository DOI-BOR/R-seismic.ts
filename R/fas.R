#' Compute a multitaper amplitude spectrum.
#'
#' \code{fas} is used to compute a multitaper Fourier amplitude spectrum
#' of an equally-spaced numerical vector.
#'
#' @param xt required (first) equally-spaced series. Must convert to numeric vector.
#' @param dt (optional) sample interval, in seconds. Default is 0.01
#' @param freq.range (optional) to limit output spectrum to. Default is 0-Nyquist
#' @return list with the amplitude spectrum and frequency increment
fas <- function(xt, dt = 0.01, freq.range = NA) {

	kind.hiRes <- 1
	kind.avgVar <- 2 # mtapspec default is 2
	nwin <- 5 # number of taper windows. mtapspec default is 5
	npi <- 3 # order of the slepian functions. mtapspec default is 3
	inorm.none <- 0 # no normalization for kind=1, or by sqrt(1/npoints) for kind=2
	inorm.len <- 1 # normalize by 1/npoints (mtaspec default)
	inorm.dt <- 2 # normalize by dt
	inorm.rtlen <- 3 # normalize by sqrt(1/npoints)
	MTP <- list(kind=kind.hiRes, nwin=nwin, npi=npi, inorm=inorm.len)
	xt.mta <- RSEIS::mtapspec(xt, dt, MTP=MTP)

	len <- length(xt.mta$spec) # same as mta$klen
	df <- 1 / (len * dt) # same as mta$df
	nf <- 1 + len / 2 # same as mta$numfreqs
	xt.spec <- xt.mta$spec[1:nf] # just need spectrum at non-negative frequencies

	f <- seq(0, nf - 1, 1) * df
	if ( ! anyNA(freq.range) ) {
		f.ok <- which( (f >= freq.range[1]) & (f <= freq.range[2]), arr.ind=TRUE )
		f <- f[f.ok]
		xt.spec <- xt.spec[f.ok]
	}

	ret <- list(aspec=xt.spec, df=df, f=f)

	return(ret)
}
