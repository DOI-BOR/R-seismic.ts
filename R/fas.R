#' Compute a Fourier amplitude spectrum.
#'
#' \code{fas} is used to compute a Fourier amplitude spectrum
#' of a univariate or multivariate time series.
#'
#' @param xt required equally-spaced time series,
#' represented as a numeric \code{\link{vector}}, \code{link{matrix}},
#' \code{\link{data.frame}}, \code{\link{ts}}, \code{\link{mts}}, or
#' \code{\link{signalSeries}} object.
#' @param dt sample interval, in seconds. Default is 0.01
#' @param freq.range to limit output spectrum to. Default is 0-Nyquist
#' @param multi.taper is set to \code{TRUE} if a multi-taper spectrum is
#' to be computed; otherwise a standard FFT is computed.
#' @param demean is set to \code{TRUE} to de-mean the input series. Default
#' is \code{FALSE}.
#' @param power Set to \code{TRUE} for power spectum, \code{FALSE} for amplitude
#' spectrum. Default is \code{FALSE}.
#' @param phase Set to \code{TRUE} for phase spectrum, \code{FALSE} for amplitude
#' or power spectrum. Default is \code{FALSE}.
#' @return A list with \code{spec, df, f}, where \code{spec} is a
#' \code{\link{vector}} (univariate input time series), or a
#' \code{\link{data.frame}} (multivariate input time series) of the
#' spectra for the non-negative frequencies. The
#' frequency increment is \code{df}, and the discrete frequencies
#' are in the \code{\link{vector}} \code{f}.
#' @keywords ts

fas <- function(xt, dt = NA, freq.range = NA, multi.taper=TRUE,
                demean=FALSE, power=FALSE, phase=FALSE) {
  if ( is(xt, "signalSeries") || is(xt, "ts") )
    dt <- deltat(xt)
  if ( ! is.finite(dt) )
    dt <- 0.01

  start <- 0
  in.units <- NA
  units.position <- NA
  if ( is(xt, "signalSeries") ) {
    in.units <- xt@units
    units.position <- xt@units.position
    start <- xt@positions@from
  }
  if ( is.na(in.units) || is.null(in.units) )
    in.units <- "mks"
  if ( is.na(units.position) || is.null(units.position) )
    units.position <- "seconds"
  if ( is(xt, "ts") )
    start <- start(xt)[1]

  MTP <- NULL
  if ( multi.taper ) {
    kind.hiRes <- 1
  	kind.avgVar <- 2 # mtapspec default is 2
  	nwin <- 5 # number of taper windows. mtapspec default is 5
  	npi <- 3 # order of the slepian functions. mtapspec default is 3
  	inorm.none <- 0 # no normalization for kind=1, or by sqrt(1/npoints) for kind=2
  	inorm.len <- 1 # normalize by 1/npoints (mtaspec default)
  	inorm.dt <- 2 # normalize by dt
  	inorm.rtlen <- 3 # normalize by sqrt(1/npoints)
  	MTP <- list(kind=kind.hiRes, nwin=nwin, npi=npi, inorm=inorm.none)
  }
  PLAN <- NULL

	multi.trace <- is.mts(xt) ||
	    ( ! is.null(dim(xt)) && length(dim(xt)) > 1 && dim(xt)[2] > 1 )
	if ( multi.trace ) {
	  if ( ! multi.taper ) {
  	  xt.len <- dim(xt)[1]
  	  PLAN <- fftw::planFFT(xt.len)
	  }
	  spec <- NULL
	  cnames <- NULL
	  for ( ii in 1:dim(xt)[2] ) {
	    if ( is(xt, "signalSeries") ) {
	      zt <- xt@data[,ii]
	      cnames <- colnames(xt@data)
	    } else {
	      zt <- xt[,ii]
	      cnames <- colnames(xt)
	    }
	    if ( demean )
	      zt <- zt - mean(zt)
	    if ( multi.taper ) {
  	    xt.mta <- RSEIS::mtapspec(zt, dt, MTP=MTP)
  	    len <- length(xt.mta$spec) # same as mta$klen
  	    df <- 1 / (len * dt) # same as mta$df
  	    nf <- 1 + len / 2 # same as mta$numfreqs
  	    if ( phase ) {
  	      xt.spec <- -atan2(xt.mta$Ispec,-xt.mta$Rspec)
  	      xt.spec <- unwrap.phase(xt.spec)
  	    } else {
    	    xt.spec <- xt.mta$spec[1:nf]
    	    if ( ! power )
    	      xt.spec <- sqrt(xt.spec)
  	    }
	    } else {
	      fft <- FFT(zt, plan=PLAN)
	      len <- length(fft)
	      df <- 1 / (len * dt)
	      nf <- 1 + len / 2
	      if ( phase ) {
	        xt.spec <- -Arg(fft[1:nf])
	        xt.spec <- unwrap.phase(xt.spec)
	      } else {
	        xt.spec <- Mod(fft[1:nf])
	        if ( power )
	          xt.spec <- xt.spec^2
	      }
	    }
	    f <- seq(0, nf - 1, 1) * df
	    if ( ! anyNA(freq.range) ) {
	      f.ok <- which( (f >= freq.range[1]) & (f <= freq.range[2]), arr.ind=TRUE )
	      f <- f[f.ok]
	      xt.spec <- xt.spec[f.ok]
	    }
	    if ( is.null(spec) )
	      spec <- xt.spec
	    else
        spec <- data.frame(spec, xt.spec)
	  }
	  colnames(spec) <- cnames
	} else {
	  if ( is(xt, "signalSeries") ) {
	    if ( ! is.null(dim(xt)) && length(dim(xt)) > 1 && dim(xt)[2] == 1 )
	      zt <- xt@data[,1]
	    else
	      zt <- xt@data
	  } else
	    zt <- xt
	  if ( demean )
	    zt <- zt - mean(zt)
	  if ( multi.taper ) {
	    xt.mta <- RSEIS::mtapspec(zt, dt, MTP=MTP)
	    len <- length(xt.mta$spec) # same as mta$klen
	    df <- 1 / (len * dt) # same as mta$df
	    nf <- 1 + len / 2 # same as mta$numfreqs
	    if ( phase ) {
	      xt.spec <- -atan2(xt.mta$Ispec,-xt.mta$Rspec)
	      xt.spec <- unwrap.phase(xt.spec)
	    } else {
	      xt.spec <- xt.mta$spec[1:nf]
	      if ( power )
	        xt.spec <- xt.spec^2
	      #xt.spec <- xt.mta$Rspec[1:nf]^2 + xt.mta$Ispec[1:nf]^2
	      #if ( ! power )
	      #  xt.spec <- sqrt(xt.spec)
	    }
	  } else {
	    zt.len <- length(zt)
	    PLAN <- fftw::planFFT(zt.len)
	    fft <- FFT(zt, plan=PLAN)
	    len <- length(fft)
	    df <- 1 / (len * dt)
	    nf <- 1 + len / 2
	    if ( phase ) {
	      xt.spec <- -Arg(fft[1:nf])
	      xt.spec <- unwrap.phase(xt.spec)
	    } else {
	      xt.spec <- Mod(fft[1:nf])
	      if ( power )
	        xt.spec <- xt.spec^2
	    }
	  }
	  f <- seq(0, nf - 1, 1) * df
	  if ( ! anyNA(freq.range) ) {
	    f.ok <- which( (f >= freq.range[1]) & (f <= freq.range[2]), arr.ind=TRUE )
	    f <- f[f.ok]
	    xt.spec <- xt.spec[f.ok]
	  }
	  spec <- xt.spec
	}

	ret <- list(spec=spec, df=df, f=f)

  return(ret)
}
