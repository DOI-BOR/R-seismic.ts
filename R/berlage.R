#' Berlage Wavelet
#'
#' @description \code{berlage} computes a Berlage wavelet.
#' @param amp Amplitude of the wavelet.
#' @param f Frequency of sine wave component, in the same units as \code{1/[dt]}.
#' @param dur Duration of the signal, in the same units as \code{[dt]}.
#' @param dt Sample interval.
#' @param A time exponent parameter for \code{berlage.fn}.
#' @param B time exponential coefficient for \code{berlage.fn}.
#' @param t Time value at which to compute \code{berlage.fn}.
#' @return A vector of Berlage function values.
#' @details Berlage uses Berlage function parameters A=0.5, and B = -3 * f.
#' @examples
#' amp <- 20 # cm/s
#' units <- "cm/s"
#' f0 <- 3.1 # Hz
#' dur <- 15.0 # sec
#' dt <- 0.01 # sec
#' zpad <- min(0.1 * dur, 3) # 10% of the duration, up to 3 second max
#' w.ber <- c(rep(0, zpad / dt), berlage(amp, f0, dur - zpad, dt))
#' wvlt <- signalSeries(w.ber,from=0,by=dt,units=units,units.position="seconds")
#' plot(wvlt)
#' @keywords ts

berlage <- function(amp,f,dur,dt) {
	len <- dur / dt
	A <- .5
	B <- -3 * f
	amp * berlage.fn(A, B, f, seq(0, dur, length=len))
}

#' @describeIn berlage \code{berlage.fn} computes the Berlage function
#' t^A * exp(B*t) * sin(2*pi*f*t) for a single time value
berlage.fn <- function(A,B,f,t) {
  t^A * exp(B*t) * sin(2*pi*f*t)
}
