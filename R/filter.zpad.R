#' Convolve a filter with a time series, after zero-padding both
#'
#' @description
#' \code{filter.zpad} takes a univariate time series of length n and a univariate
#' filter of length m, zero-pads both of them to length n + m - 1, and then filters
#' the two zero-padded series. This avoids wrap-around effects that otherwise may occur.
#'
#' @param x Equally-sampled univariate time series. Must convert to a numeric vector.
#' @param filt Equally-sampled univariate filter series. Must convert to a numeric vector.
#' @param truncate If true, truncate output to input data length (def. FALSE).
#' @return The filtered time series.
#' @seealso \code{\link{filter}}
#' @keywords ts

filter.zpad <- function(x, filt, truncate = F) {
	x.len <- length(x)
	filt.len <- length(filt)
	aug.len <- x.len + filt.len - 1 # actual length of filtered time series
	if ( truncate )
		out.len <- x.len # truncate output to input data length
	else
		out.len <- aug.len # full output
	filter(zero.pad(x,aug.len), zero.pad(filt,aug.len), method = "convolution", sides = 1, circular = T)[1:out.len]
}
