# zero-pad input and filter to avoid wrap-around effects
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
