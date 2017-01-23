ricker.fn <- function(sig,t) {
	t2 <- (t/sig)^2
	(2./(pi^.25)) * (1. - t2) * exp(-0.5 * t2) / sqrt(3.*sig)
}

ricker <- function(amp,dur,dt) {
	half.len = 2 * 4
	sig <- 0.5 * dur / half.len
	len <- dur / dt
	a0 <- amp / ricker.fn(sig, 0)
	a0 * ricker.fn(sig, seq(-half.len*sig, half.len*sig, length=len))
}
