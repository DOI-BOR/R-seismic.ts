#' Ricker Wavelet
#'
#' @description \code{ricker} computes a Ricker wavelet.
#' @param amp Amplitude of the wavelet.
#' @param dur Duration of the signal, in the same units as \code{[dt]}.
#' @param dt Sample interval.
#' @param sig Sigma value for for \code{ricker.fn}.
#' @param t Time value at which to compute \code{ricker.fn}.
#' @return A vector of Ricker function values.
#' @details Ricker uses Ricker.fn function parameter sig / 16.
#' @keywords ts

ricker <- function(amp,dur,dt) {
	half.len = 2 * 4
	sig <- 0.5 * dur / half.len
	len <- dur / dt
	a0 <- amp / ricker.fn(sig, 0)
	a0 * ricker.fn(sig, seq(-half.len*sig, half.len*sig, length=len))
}

ricker.fn <- function(sig,t) {
  t2 <- (t/sig)^2
  (2./(pi^.25)) * (1. - t2) * exp(-0.5 * t2) / sqrt(3.*sig)
}
