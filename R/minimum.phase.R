#' Minimum-phase representation of an input time series
#'
#' \code{minimum.phase} returns the minimum-phase representation of an
#' input univariate time series
#'
#' @param x.data Equally-sampled input series. Must convert to numeric vector.
#' @param zero.pad Zero-pad the input using \code{\link{zero.pad}} before
#' transforming. Default is TRUE.
#'
#' @details This function uses a spectral factorization derived from the
#' Kramers-Kronig relation to determine the minimum-phase spectrum from
#' the amplitude spectrum of the data. The Kramers-Kronig relation shows that
#' the real and imaginary parts of an analytic spectrum form a Hilbert
#' transform pair, and therefore are not independent. The minimum phase spectrum
#' therefore can be determined from the Hilbert transform of the log of the
#' amplitude spectrum. See Buttkus (2000), eqn. 13.95, or Robinson and Trietel
#' (2000) eqn 9.43.
#' @return The minimum-phase representation.
#' @seealso \code{\link{zero.pad}}, \code{\link{hilbert}}
#' \itemize{
#' \item \href{https://en.wikipedia.org/wiki/Kramers-Kronig_relations}{Wikipedia}
#' \item \href{http://www.springer.com/us/book/9783540626749}{Buttkus (2000)}
#' Spectral Analysis and Filter Theory in Applied Geophysics.
#' \item \href{http://library.seg.org/doi/book/10.1190/1.9781560802327}{Robinson and Tritel (2000)}
#' Geophysical Signal Analysis.
#' }
#' @keywords ts

minimum.phase <- function(x.data, zero.pad=TRUE) {
	if ( zero.pad )
		x <- zero.pad(x.data)
	else
		x <- as.vector(x.data)
	x.fft <- fft(x) / length(x)
	x.aspec <- Mod(x.fft)
	x.min.phase <- hilbert(log(x.aspec), "hilbert", zero.pad)
	x.min.phase.fft <- x.aspec * exp(-1i * x.min.phase)
	# tsplot(x.min.phase)
	x.min.phase.ifft <- fft(x.min.phase.fft, inverse=TRUE)
	x.min.phase <- Re(x.min.phase.ifft)
	return (x.min.phase)
}
