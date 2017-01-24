tsplot <- function(x, dt=.01, ...) {
	x <- as.vector(x)
	t <- seq(length(x)) * dt
	plot(t, x, type="l", ...)
}