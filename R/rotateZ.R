#' Rotate 3-component data about the Z axis
#'
#' \code{rotateZ} rotates a 3-component time-series about the positive Z axis.
#'
#' @param xt,yt,zt equally-spaced time series in the X, Y, and Z directions. At
#' a minimum, the X and Y components must be present. If inputs \code{xt},
#' \code{yt}, and (optionally) \code{zt} are present, then they must be
#' equal-length univariate time series of the coordinates.
#' Alternatively, if \code{xt} is a 2- or 3-component multivariate time series,
#' including \code{\link{matrix}}, \code{\link{data.frame}}, \code{\link{ts}},
#' \code{\link{mts}} or \code{\link{signalSeries}}, then the Y and Z components
#' are taken from \code{xt}, and \code{yt} and \code{zt} are not used.
#' @param phi rotation angle about positive Z axis (from positive X axis towards
#' positive Y axis), in radians
#' @param cnames vector of column names to use for the rotated components.
#' Default is to append the rotation angle, in degrees, to the original
#' component name (or XYZ if the original component names are not set).
#' @param dt Sample interval, in seconds. Not used unless input is a
#' \code{\link{ts}} or \code{\link{signalSeries}} object, and the time
#' step of the object is not set. Default is 0.01 seconds.
#'
#' @details Computes \code{R = X * cos(phi) + Y * sin(phi)},
#' and \code{T = -X * sin(phi) + Y * cos(phi)}. The \code{XYZ} directions for
#' the 3-component seismogram must form a right-handed coordinate system.
#' @return The rotated time series.
#' @seealso \code{\link{ts}}, \code{\link{signalSeries}}
#' @examples
#' deg2rad <- pi / 180
#' # vector example
#' xt <- (1:4) / sqrt(2)
#' yt <- (1:4) / sqrt(2)
#' zt <- rep(1,4)
#' rotateZ(xt, yt, zt, phi=45*deg2rad)
#' #
#' # matrix example
#' data <- c(rep(sqrt(2),8), rep(1,4))
#' dim(data) <- c(4,3)
#' colnames(data) <- c("my.X","my.Y","my.Z")
#' rotateZ(xt=data, phi=45*deg2rad)
#' rotateZ(xt=data, phi=45*deg2rad, cnames=c("XR","YR","Z"))
#' #
#' # ts/mts examples
#' rotateZ(ts(xt),ts(yt), phi=45*deg2rad)
#' xt.mts <- ts(cbind(my.X=xt,my.Y=yt,my.Z=zt),start=0,deltat=0.01)
#' rotateZ(xt=xt.mts, phi=45*deg2rad)
#' #
#' # signalSeries example
#' xt.ss <- signalSeries(cbind(xt,yt,zt),
#'        from=0, by=0.01, units="cm", units.position="sec")
#' rotateZ(xt.ss, phi=45*deg2rad, cnames=c("XR","YR","Z"))
#'
#' @keywords ts
deg2rad <- pi / 180
rad2deg <- 180 / pi

#' @describeIn rotateZ.default rotates \code{vector} univariate
#' time series, or a \code{matrix} multivariate time series.
rotateZ.default <- function(xt, yt, zt, phi, cnames=NA, dt=NA) {
  if ( missing(xt) || missing(phi) )
    stop("Must provide input xt and phi")

  multi.trace <- is.matrix(xt) || length(dim(xt)) > 1

  # split out the components
  Z <- NULL
  cnames.new <- NULL
  if ( multi.trace ) {
    xt.len <- dim(xt)[1]
    n.traces <- dim(xt)[2]
    if ( n.traces < 2 )
      stop("multi-trace time series must have at least 2 components")
    X <- as.vector(xt[,1])
    Y <- as.vector(xt[,2])
    if ( n.traces >= 3 )
    Z <- as.vector(xt[,3])
    cnames.def <- colnames(xt)
    if ( any(is.null(cnames.def)) ) {
      if ( n.traces >= 3 )
        cnames.def <- c("X","Y","Z")
      else
        cnames.def <- c("X","Y")
    }
  } else {
    if ( missing(yt) )
      stop("Must provide input xt, yt and phi")
    else if ( length(xt) != length(yt) )
      stop("xt and yt must be the same length")
    X <- as.vector(xt)
    Y <- as.vector(yt)
    cnames.def <- c("X","Y")
    if ( ! missing(zt) && ! is.null(zt) ) {
      if ( length(xt) != length(zt) )
        stop("xt and zt must be the same length")
      Z <- as.vector(zt)
      cnames.def <- c("X","Y","Z")
    }
  }

  # get the rotated components
  RX <- X * cos(phi) + Y * sin(phi)
  RY <- -X * sin(phi) + Y * cos(phi)

  # bind rotated components into a matrix object and return
  if ( is.null(Z) )
    RXYZ <- cbind(RX, RY)
  else
    RXYZ <- cbind(RX, RY, Z)

  if ( anyNA(cnames) ) {
    # get phi as a non-negative rotation, in degrees, and append to default cnames
    phi.deg <- as.integer(round(phi * 180 / pi)) %% 360
    if ( phi.deg < 0 )
      phi.deg = phi.deg + 360
    phi.s <- sprintf("%d",phi.deg)
    cnames <- c(paste(cnames.def[1:2], phi.s, sep="."),cnames.def[3])
  }
  colnames(RXYZ) <- cnames

  return(RXYZ)
}
setGeneric("rotateZ",def=rotateZ.default)

#' @describeIn rotateZ.default rotates \code{ts} univariate or \code{mts}
#' multivariate time series
rotateZ.ts <- function(xt, yt, zt, phi, cnames=NA, dt=NA) {
  if ( missing(xt) || missing(phi) )
    stop("Must provide input xt and phi")

  multi.trace <- is.mts(xt)
  dt <- deltat(xt)
  if ( is.na(dt) )
    dt <- 0.01
  start <- start(xt)[1]

  # split out the components
  Z <- NULL
  cnames.def <- NULL
  if ( multi.trace ) {
    xt.len <- dim(xt)[1]
    n.traces <- dim(xt)[2]
    if ( n.traces < 2 )
      stop("multi-trace time series must have at least 2 components")
    X <- xt[,1]
    Y <- xt[,2]
    if ( n.traces >= 3 )
      Z <- xt[,3]
    cnames.def <- colnames(xt)
  } else {
    if ( missing(yt) )
      stop("Must provide input xt, yt and phi")
    else if ( length(xt) != length(yt) )
      stop("xt and yt must be the same length")
    X <- xt
    Y <- yt
    cnames.def <- c(colnames(xt), colnames(yt))
    if ( ! missing(zt) ) {
      Z <- zt
      cnames.def <- c(cnames.def, colnames(zt))
    }
  }
  if ( anyNA(cnames) && ! any(is.null(cnames.def)) ) {
    # get phi as a non-negative rotation, in degrees, and append to default cnames
    phi.deg <- as.integer(round(phi * 180 / pi)) %% 360
    if ( phi.deg < 0 )
      phi.deg = phi.deg + 360
    phi.s <- sprintf("%d",phi.deg)
    cnames <- c(paste(cnames.def[1:2], phi.s, sep="."),cnames.def[3])
  }

  # put the rotated components into a ts object
  RTZ <- ts(rotateZ.default(X,Y,Z,phi,cnames), start=start, deltat=dt)

  return(RTZ)
}
setMethod("rotateZ","ts",rotateZ.ts)

#' @describeIn rotateZ.default rotates a \code{signalSeries} univariate or
#' multivariate time series
rotateZ.signalSeries <- function(xt, yt, zt, phi, cnames=NA, dt=NA) {
  if ( missing(xt) || missing(phi) )
    stop("Must provide input xt and phi")

  multi.trace <- ! is.null(dim(xt))
  if ( is.na(dt) )
    dt <- deltat(xt)
  start <- xt@positions@from
  units <- xt@units

  # split out the X and Y components
  Z <- NULL
  cnames.def <- NULL
  if ( multi.trace ) {
    xt.len <- dim(xt)[1]
    n.traces <- dim(xt)[2]
    if ( n.traces < 2 )
      stop("multi-trace time series must have at least 2 components")
    X <- xt[,1]@data
    Y <- xt[,2]@data
    if ( n.traces >= 3 )
      Z <- xt[,3]@data
    cnames.def <- colnames(xt@data)
  } else {
    if ( missing(yt) )
      stop("Must provide input xt, yt and phi")
    else if ( length(xt) != length(yt) )
      stop("xt and yt must be the same length")
    X <- xt@data
    Y <- yt@data
    cnames.def <- c(colnames(xt@data), colnames(yt@data))
    if ( ! missing(zt) ) {
      Z <- zt@data
      cnames.def <- c(cnames.def, colnames(zt@data))
    }
  }
  if ( anyNA(cnames) && ! any(is.null(cnames.def)) ) {
    # get phi as a non-negative rotation, in degrees, and append to default cnames
    phi.deg <- as.integer(round(phi * 180 / pi)) %% 360
    if ( phi.deg < 0 )
      phi.deg = phi.deg + 360
    phi.s <- sprintf("%d",phi.deg)
    cnames <- c(paste(cnames.def[1:2], phi.s, sep="."),cnames.def[3])
  }

  # put the rotated components into a signalSeries object
  RTZ <- signalSeries(rotateZ.default(X,Y,Z,phi,cnames),
                      from=start, by=dt, units=units)

  return(RTZ)
}
setMethod("rotateZ","signalSeries",rotateZ.signalSeries)


