#' Apply a Generalized Hann (Tukey) Window to a Time Series
#'
#' @description
#' \code{hanning} takes an input time series and multiplies it with a
#' generalized Hann (Tukey) window. The input time series can be univariate
#' or multivariate.
#'
#' @param x.data Equally-sampled univariate or multivariate series. Supported
#' types include a numeric \code{\link{vector}}, \code{\link{matrix}},
#' \code{\link{data.frame}}, \code{\link{ts}}, or \code{\link{signalSeries}}.
#' @param pct Percentage of data window to apply a taper. Must be
#' between 0 and 50. Default is 0 (no taper).
#' @param demean Should windowed data be demeaned? Default is F.
#' @return The windowed time series.
#' @details For factor > 1, this code augments by zero-padding the Fourier
#' transform of the input series at frequencies greater than the Nyquist,
#' and then taking the inverse transform.
#' @seealso \code{\link{windowTs.default}}
#' @keywords ts

hanning <- function(x.data, pct = NA, demean = NA) {

  if ( is.na(demean) )
    demean = FALSE

	multi.trace <- is.mts(x.data) || is.matrix(x.data) || length(dim(x.data)) > 1 ||
			( is(x.data, "signalSeries") && ! is.null(dim(x.data)) )

	if ( multi.trace == TRUE ) {
    if ( is(x.data, "signalSeries") )
      x.len <- dim(x.data)[1]
    else
      x.len <- dim(x.data)[1]
	} else
		x.len <- length(x.data)

	if ( demean && is.numeric(x.data) ) {
		if ( multi.trace == TRUE ) {
			for ( cn in 1:dim(x.data)[2] ) {
			  if ( is(x.data, "signalSeries") )
			    x.data@data[,cn] <- x.data@data[,cn] - mean(x.data@data[,cn])
			  else
			    x.data[,cn] <- x.data[,cn] - mean(x.data[,cn])
			}
		} else {
		  if ( is(x.data, "signalSeries") )
		    x.data@data <- x.data@data - mean(x.data)
		  else
		    x.data <- x.data - mean(x.data)
		}
	}

	if ( is.na(pct) || pct == 50. ) {
		# simple Hann window
		# taper <- 0.5*(1 - cos(seq(0,2*pi*(x.len-1)/x.len,length=x.len)))
		taper <- 0.5*(1 - cos(seq(0,2*pi,length=x.len)))
	} else {
		# Tukey window
		min.pct <- 100. / x.len # require at least 1 point to taper
		pct <- max(min.pct,min(pct,50.))
		tap.len <- as.integer(pct * x.len / 100.)
		taper <- 0.5 * c(
			# 1 - cos(seq(0, pi*(tap.len-1)/tap.len, length=tap.len)),
			# 1 - cos(seq(0, pi*(tap.len-.5)/tap.len, length=tap.len)),
			1 - cos(seq(0, pi, length=tap.len)),
			rep(2, x.len - 2 * tap.len),
			# 1 - cos(seq(pi, pi*(2*tap.len-1)/tap.len, length=tap.len))
			# 1 - cos(seq(pi*(tap.len+.5)/tap.len, 2 * pi, length=tap.len))
			1 - cos(seq(pi, 2 * pi, length=tap.len))
		)
	}

	if ( is(x.data, "signalSeries") )
	  x.data@data <- taper * x.data@data
	else
	  x.data <- taper * x.data

	return(x.data)
}
