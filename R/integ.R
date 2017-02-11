#' Integrate a time series
#'
#' \code{integ} integrates a univariate or multivariate time series.
#'
#' @param x.data Equally-sampled input series. Must convert to a numeric
#' \code{\link{vector}}, \code{\link{signalSeries}}, or \code{\link{ts}}.
#' @param dt Sample interval. Default is 0.01 seconds.
#'
#' @details These are generic functions. Integration is performed by \code{\link{cumsum}}.
#' @return List containing the cumulative sum, and the integrated time series.
#' @seealso \code{\link{cumsum}}
#' @keywords ts
#'

integ.default <- function(x.data, dt=0.01) {
	multi.trace <- is.matrix(x.data) || length(dim(x.data)) > 1

	ss <- NULL
	if ( multi.trace ) {
		x.len <- dim(x.data)[1]
		xi.data = NULL
		for ( cn in 1:dim(x.data)[2] ) {
			ok <- ! is.na(x.data[,cn])
			Ix <- cumsum(x.data[ok,cn]) * dt
			ss <- c(ss, Ix[length(Ix)])
			if ( is.null(xi.data) )
				xi.data <- data.frame(Ix)
			else
				xi.data <- data.frame(xi.data,Ix)
		}
	} else {
		x.len <- length(x.data)
		ok <- ! is.na(x.data)
		Ix <- cumsum(x.data[ok]) * dt
		ss <- c(ss, Ix[length(Ix)])
		xi.data <- Ix
	}
	colnames(xi.data) <- colnames(x.data)
	list(sum=ss,Ix.data=xi.data)
}
setGeneric("integ",def=integ.default)

#' @describeIn integ.default integrates a \code{ts}
integ.ts <- function(x.data, dt=NA) {
	multi.trace <- is.mts(x.data)

	if ( is.na(dt) )
		dt <- deltat(x.data)
	start <- start(x.data)[1]
	ss <- NULL
	if ( multi.trace == TRUE ) {
		x.len <- dim(x.data)[1]
		xi.data = NULL
		for ( cn in 1:dim(x.data)[2] ) {
			ok <- ! is.na(x.data[,cn])
			Ix <- cumsum(x.data[ok,cn]) * dt
			ss <- c(ss, Ix[length(Ix)])
			if ( is.null(xi.data) )
				xi.data <- ts(Ix, start = start, deltat = dt)
			else
				xi.data <- ts(xi.data, Ix, start = start, deltat = dt)
		}
	} else {
		x.len <- length(x.data)
		ok <- ! is.na(x.data)
		Ix <- cumsum(x.data[ok]) * dt
		ss <- c(ss, Ix[length(Ix)])
		xi.data <- Ix
	}
	names(xi.data) <- names(x.data)
	list(sum=ss,Ix.data=xi.data)
}
setMethod("integ","ts",integ.ts)

#' @describeIn integ.default integrates a \code{signalSeries}
integ.signalSeries <- function(x.data, dt=NA) {
	multi.trace <- ! is.null(dim(x.data))

	if ( is.na(dt) )
		dt <- deltat(x.data)
	start <- x.data@positions@from

	units <- x.data@units
	if ( grep("cm/s",units) )
		new.units <- "cm"
	else if ( grep("cm/s/s",units) )
		new.units <- "cm/s"
	else if ( grep("ft/s",units) )
		new.units <- "ft"
	else if ( grep("ft/s/s",units) )
		new.units <- "ft/s"
	else
		new.units <- NULL

	ss <- NULL
	if ( multi.trace == TRUE ) {
		x.len <- dim(x.data)[1]
		xi.data = NULL
		for ( cn in 1:dim(x.data)[2] ) {
			ok <- ! is.na(x.data[,cn]@data)
			Ix <- cumsum(x.data[ok,cn])
			Ix@data <- Ix@data * dt
			ss <- c(ss, Ix@data[length(Ix)])
			if ( is.null(xi.data) )
				xi.data <- signalSeries(Ix, from = start, by = dt, units = new.units)
			else
				xi.data <- signalSeries(xi.data, Ix, from = start, by = dt, units = new.units)
		}
	} else {
		x.len <- length(x.data)
		ok <- ! is.na(x.data@data)
		Ix <- cumsum(x.data[ok])
		Ix@data <- Ix@data * dt
		ss <- c(ss, Ix@data[length(Ix)])
		xi.data <- signalSeries(Ix, from = start, by = dt, units = new.units)
	}
	names(xi.data) <- names(x.data)
	list(sum=ss,Ix.data=xi.data)
}
setMethod("integ","signalSeries",integ.signalSeries)

