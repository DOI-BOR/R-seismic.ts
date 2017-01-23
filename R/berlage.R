berlage.fn <- function(A,B,f,t) {
	t^A * exp(B*t) * sin(2*pi*f*t)
}

berlage <- function(amp,f,dur,dt) {
	len <- dur / dt
	A <- .5
	B <- -3 * f
	amp * berlage.fn(A, B, f, seq(0, dur, length=len))
}
