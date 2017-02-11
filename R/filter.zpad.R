#' Convolve a filter with an input timeseries, after zero-padding both
#'
#' @description
#' \code{filter.zpad} takes an input time series and a filter, zero-pads
#' both of them, and then convolves the two series. This avoids wrao-around
#' effects that otherwise may occur.
#'
#' @param x Equally-sampled input series. Must convert to a numeric vector.
#' @param filt Equally-sampled filter series. Must convert to a numeric vector.
#' @return The convolved time series.
#' @details For factor > 1, this code augments by zero-padding the Fourier
#' transform of the input series at frequencies greater than the Nyquist,
#' and then taking the inverse transform.
#' @seealso \code{\link{filter}}
#' @keywords ts

filter.zpad <- function(x, f, truncate = F) {
	x.len <- length(x)
	f.len <- length(f)
	aug.len <- x.len + f.len - 1 # actual length of filtered time series
	if ( truncate )
		out.len <- x.len # truncate output to input data length
	else
		out.len <- aug.len # full output
	filter(zero.pad(x,aug.len), zero.pad(f,aug.len), method = "convolution", sides = 1, circular = T)[1:out.len]
}
