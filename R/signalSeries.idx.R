setMethod("[","signalSeries",
function (x, i, j, ..., drop = TRUE)
{
	nr <- nrow(x@data)
	nc <- ncol(x@data)
	if ( is.null(nr) ) {
		nr = length(x@data)
		nc = 1
	}
	if ( is.null(nc) )
		nc = 1
	if (missing(i))
		i <- seq(nr)
	if (missing(j))
		j <- if ( nc > 1 ) seq(nc) else 1
	if (i < 0 || i > nr)
		stop("index i is out of range")
	if (j < 0 || j > nc)
		stop("index j is out of range")
	pos <- as(positions(x), "numeric")[i]
	data <- if ( nc > 1 ) x@data[i,j] else x@data[i]
	z <- signalSeries(data = data, positions. = as(pos, "numericSequence"),
										units = x@units, units.position = x@units.position)
	z@title <- x@title
	z@documentation <- x@documentation
	z
}
)
