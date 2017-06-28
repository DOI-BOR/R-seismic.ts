#' Plot a signalSeries Object
#'
#' This function plots \code{\link{signalSeries}} onjects from
#' package \pkg{splus2R}, and provides additional functionality for
#' multivariate time series.
#'
#' @param x,y \code{\link{signalSeries}} object(s) from which to extract element(s)
#' or in which to replace element(s).
#' @param ... arguments to \code{\link{plot}}.
setMethod("plot","signalSeries",
function (x, y, ..., main = NULL, ylab = x@units[1], xlab = x@units.position,
					top.ticks = FALSE, right.ticks = FALSE, reference.grid = TRUE,
					merge.args = list(pos = "union", how = "interp"), x.axis.args = list(),
					y.axis.args = list(), plot.args = list(), log.axes = "",
					complex.convert = Mod, dB = FALSE, frame = sys.nframe(),
					col = NULL, lty = 1, lwd = 1, type = "l", cex.main = 1, legend = NA,
					legend.pos = NA)
{
	mergeArgList <- function(x, y) {
		if (!is.list(y))
			stop("y must be a list")
		if (is.null(x))
			return(y)
		if (!is.list(x))
			stop("x must be a list")
		x[names(y)] <- y
		x
	}
	y <- c(list(x), if (!missing(y)) list(y), list(...))
	ny <- length(y)
	ylim <- range(unlist(lapply(y, function(x) range(seriesData(x)))))
	xlim <- range(unlist(lapply(y, function(x) range(as(positions(x),
																											"numeric")))))
	plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab, log = log.axes)
	if (is.null(main))
		main <- paste(unlist(lapply(y, function(x) x@title)), sep = ",")
	leg.lty <- NULL
	leg.col <- NULL
	leg.lwd <- NULL
	leg.type <- NULL
	for (i in seq(along = y)) {
		xdata <- as(positions(y[[i]]), "numeric")
		ns <- ncol(y[[i]]@data)
		if ( is.null(ns) ) {
			mts <- FALSE
			ns <- 1
		} else {
			mts <- TRUE
		}
		if (is.null(col))
			col <- seq(ns)
		else
			col <- rep(col, length = ns)
		lty <- rep(lty, length = ns)
		lwd <- rep(lwd, length = ns)
		type <- rep(type, length = ns)
		for ( j in seq(ns) ) {
			ydata <- if ( mts ) y[[i]]@data[,j] else y[[i]]@data
			do.call("lines", mergeArgList(list(x = xdata, y = ydata,
																			 col = col[j], lty = lty[j], lwd = lwd[j],
																			 type = type[j]), plot.args))
			leg.lty <- c(leg.lty, lty[j])
			leg.col <- c(leg.col, col[j])
			leg.lwd <- c(leg.lwd, lwd[j])
			leg.type <- c(leg.type, type[j])
		}
	}
	title(main = main, cex.main = cex.main)
	if ( ! is.null(ns) && anyNA(legend) )
	  legend <- colnames(x@data)

	if ( ! is.null(legend) && ! anyNA(legend) ) {
		if ( is.na(legend.pos) )
			legend.pos <- "topleft"
		legend(x=legend.pos, legend=legend, lty=leg.lty, col=leg.col,
					 lwd=leg.lwd, plot.args, cex=0.8)
	}
	invisible(NULL)
}
)

