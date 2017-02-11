#' Ricker Wavelet
#'
#' @description \code{ricker} computes a Ricker wavelet for a specified duration
#' and amplitude.
#' @param amp Amplitude of the wavelet.
#' @param dur Duration of the signal, in the same units as \code{[dt]}.
#' @param dt Sample interval.
#' @return A vector of Ricker function values.
#' @details \code{ricker} uses \code{ricker.fn} with \code{sig=dur/16} for
#' the time vector \code{t <- seq(0,dur,length=dur/dt)}.
#' @keywords ts

ricker <- function(amp,dur,dt) {
	half.len = 2 * 4
	sig <- 0.5 * dur / half.len
	len <- dur / dt
	a0 <- amp / ricker.fn(sig, 0)
	a0 * ricker.fn(sig, seq(-half.len*sig, half.len*sig, length=len))
}

#' @describeIn ricker computes the Ricker function
#' \eqn{2/\pi^{1/4} (1 - {t/\sigma}^2) e^{-0.5 {t/\sigma}^2} / \sqrt{3 \sigma}}.
#' @param sig \eqn{\sigma} value for for \code{ricker.fn}.
#' @param t Time value or vector at which to compute \code{ricker.fn}.
ricker.fn <- function(sig,t) {
  t2 <- (t/sig)^2
  (2./(pi^.25)) * (1. - t2) * exp(-0.5 * t2) / sqrt(3.*sig)
}
