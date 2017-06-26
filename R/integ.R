#' Integrate a time series
#'
#' @description
#' \code{integ} performs discrete-time integration of a univariate or
#' multivariate time series.
#'
#' @param x.data Equally-sampled input time series, including a numeric \code{\link{vector}},
#' \code{\link{matrix}}, \code{\link{data.frame}}, \code{\link{ts}}, \code{\link{mts}}
#' or \code{\link{signalSeries}}.
#' @param dt Sample interval, in seconds. Normally only used for \code{\link{vector}},
#' \code{\link{matrix}}, or \code{\link{data.frame}} input objects. If input is a
#' \code{\link{ts}}, \code{\link{mts}}, or \code{\link{signalSeries}}, then the sample
#' interval will be taken from the object attributes. However, if \code{\link{deltat}}
#' is not set, then \code{dt} will be used. Default is 0.01 seconds.
#' @param order Order of the discrete-time integrator to use: 0 = Backwards
#' rectangular, 1 = Trapezoidal, 2 = Simpson. Default is 2.
#'
#' @details These are generic functions. Implementation is in R,
#' using appropriate calls to \code{\link{cumsum}}.
#' @return List containing the cumulative sum, and the integrated time series
#' (which has the same object type as the input).
#' @seealso \code{\link{ts}}, \code{\link{mts}}, \code{\link{signalSeries}},
#' \code{\link{deltat}}, \code{\link{cumsum}}, and \code{\link{deriv.fd}}.
#' @keywords ts
#'

#' @describeIn integ.default integ.default integrates a \code{\link{vector}}
#' univariate time series, and \code{\link{matrix}} or \code{\link{data.frame}}
#' multivariate time series.
integ.default <- function(x.data, dt=NA, order=NA) {
  if ( is.na(dt) )
    dt <- 0.01

  multi.trace <- is.matrix(x.data) || length(dim(x.data)) > 1

	ss <- NULL
	if ( multi.trace ) {
		ix.mts = NULL
		for ( cn in 1:dim(x.data)[2] ) {
			ok <- ! is.na(x.data[,cn])
			ix <- dt_integ(x.data[ok,cn], order) * dt
			ss <- c(ss, ix[length(ix)])
			if ( is.null(ix.mts) )
			  ix.mts <- data.frame(ix)
			else
			  ix.mts <- data.frame(ix.mts, ix)
		}
	  colnames(ix.mts) <- colnames(x.data)
	  Ix <- ix.mts
	} else {
		ok <- ! is.na(x.data)
		Ix <- dt_integ(x.data[ok], order) * dt
		ss <- c(ss, Ix[length(Ix)])
	}
	ret <- list(sum=ss, Ix.data=Ix)

	return(ret)
}
setGeneric("integ",def=integ.default)

#' @describeIn integ.default integrates a \code{\link{ts}} or \code{\link{mts}}
#' univariate or multivariate time series
integ.ts <- function(x.data, dt=NA, order=NA) {
	multi.trace <- is.mts(x.data)
	dt <- deltat(x.data)
	if ( is.na(dt) )
	  dt <- 0.01
	start <- start(x.data)[1]

		ss <- NULL
	if ( multi.trace == TRUE ) {
		ix.mts = NULL
		for ( cn in 1:dim(x.data)[2] ) {
			ok <- ! is.na(x.data[,cn])
			ix <- dt_integ(x.data[ok,cn], order) * dt
			ss <- c(ss, ix[length(ix)])
			if ( is.null(ix.mts) )
			  ix.mts <- data.frame(ix)
			else
			  ix.mts <- data.frame(ix.mts, ix)
		}
	  colnames(ix.mts) <- colnames(x.data)
	  Ix <- ix.mts
	} else {
		ok <- ! is.na(x.data)
		Ix <- dt_integ(x.data[ok], order) * dt
		ss <- c(ss, Ix[length(Ix)])
	}
	Ix <- ts(Ix, start=start, deltat=dt)

	ret <- list(sum=ss, Ix.data=Ix)

	return(ret)
}
setMethod("integ","ts",integ.ts)

#' @describeIn integ.default integrates a \code{\link{signalSeries}}
#' univariate or multivariate time series
integ.signalSeries <- function(x.data, dt=NA, order=NA) {

	multi.trace <- ! is.null(dim(x.data))

	dt <- deltat(x.data)
	if ( is.na(dt) )
	  dt <- 0.01
	start <- x.data@positions@from
	units.position <- x.data@units.position

	units <- x.data@units
	if ( is.na(units) || is.null(units) ) {
	  new.units <- NULL
	} else if ( units == "cm/s/s" || units == "cm/s^2" ) {
		new.units <- "cm/s"
	} else if ( units == "cm/s" ) {
		new.units <- "cm"
	} else if ( units == "m/s/s" || units == "m/s^2" ) {
	  new.units <- "m/s"
	} else if ( units == "m/s" ) {
	  new.units <- "m"
	} else if ( units == "ft/s/s" || units == "ft/s^2" ) {
		new.units <- "ft/s"
	} else if ( units == "ft/s" ) {
		new.units <- "ft"
	} else {
		new.units <- NULL
	}

	ss <- NULL
	if ( multi.trace == TRUE ) {
		ix.mts = NULL
		for ( cn in 1:dim(x.data)[2] ) {
			ok <- ! is.na(x.data@data[,cn])
			ix <- dt_integ(x.data@data[ok,cn], order) * dt
			ss <- c(ss, ix[length(ix)])
			if ( is.null(ix.mts) )
			  ix.mts <- data.frame(ix)
			else
			  ix.mts <- data.frame(ix.mts, ix)
		}
		colnames(ix.mts) <- colnames(x.data@data)
		Ix <- ix.mts
	} else {
		ok <- ! is.na(x.data@data)
		Ix <- dt_integ(x.data@data[ok], order) * dt
		ss <- c(ss, Ix[length(Ix)])
	}
	Ix <- signalSeries(Ix, from=start, by=dt, units=new.units,
	                        units.position=units.position)

	ret <- list(sum=ss,Ix.data=Ix)

	return(ret)
}
setMethod("integ","signalSeries",integ.signalSeries)

# discrete-time integration filters (digital integrators) for order 0-2
dt_integ <- function(xt, order=NA) {
  xt <- as.numeric(xt)

  # use default if order not set; otherwise, silently force 0 <= order <= 2
  if ( is.na(order) )
    order <- 2
  else
    order <- max(0,min(as.integer(order),2))

  # get series length, and handle 0-length case
  lx <- length(xt)
  if ( lx == 0 )
    return( c(0.) )

  # for length(xt) > 0, silently force order <= length(xt) - 1
  order <- min(order, lx - 1)

  # use cumsum to implement the discrete integration filter
  Ix <- NA
  if ( order == 0 ) {
    # backward rectangular integrator
    Ix <- cumsum(xt)
  } else if ( order == 1 ) {
    # trapezoid integrator
    ax <- 0.5 * (c(xt[1],xt) + c(xt,xt[lx]))
    Ix = cumsum(ax[1:lx])
    # note: dropped last term to preseve length
  } else if ( order == 2 ) {
    # simpson's integrator
    ax <- (c(xt,xt[lx],2*xt[lx]) + 4. * c(0,xt,0) + c(2*xt[1],xt[1],xt)) / 3.
    # Need to sum separately over the even and odd indices. Note: this code
    # handles both even and odd-length input series
    lx2_odd <- as.integer(floor((lx + 1) / 2))
    lx2_even <- as.integer(floor(lx / 2))
    iodd <- seq(1, 2 * lx2_odd - 1, 2)
    ieven <- seq(2, 2 * lx2_even, 2)
    Ix[iodd] <- cumsum(ax[iodd])
    Ix[ieven] <- cumsum(ax[ieven])
    Ix <- Ix[1:lx]
    # note: dropped last 2 terms to preseve length
  }
  Ix
}
