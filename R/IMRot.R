#' Orientetation-independent response spectra, spectrum intensities, and peak motions
#'
#' @description
#' \code{IMRot} computes the orientation-independent intensity measures of
#' 2-component ground motion time histories, including peak values, spectrum
#' intensity, and response spectra.
#'
#' @param xt,yt equally-spaced time series in the X and Y directions. If inputs
#' \code{xt} and \code{yt} are both present, then they must be
#' equal-length univariate time series.
#' Alternatively, if \code{xt} is a 2- or 3-component multivariate time series,
#' including \code{\link{matrix}}, \code{\link{data.frame}}, \code{\link{ts}},
#' \code{\link{mts}} or \code{\link{signalSeries}}, then the X and Y components
#' are taken from \code{xt}, and \code{yt} is not used. In this case, the X and Y
#' components are taken as the first and second components, respectively, unless
#' the column names use the SEED naming convention (e.g., BHZ, BHN, BHE), in which
#' case the X and Y components are assigned from the components having E and N names.
#' The Z component is ignored.
#' @param dt Sample interval, in seconds. Not used for
#' \code{\link{ts}} or \code{\link{signalSeries}} inputs, unless the time
#' step of the object is not set. Default is 0.01 seconds.
#' @param units.ts Units of time series. One of \code{c("cgs", "mks", "SI")}.
#' Default is \code{"cgs"}, or the units specified for \code{\link{signalSeries}} objects.
#' @param ts.type String indicating input type: \code{c("acceleration", "velocity",
#' "displacement")}. Can be abbreviated to just one character, and is case independent.
#' Default is \code{"velocity"}.
#' @param rs.type Response spectrum type: \code{c("acceleration", "velocity",
#' "displacement")} for pseudo-absolute acceleration, psuedo-relative velocity, or
#' relative displacement response spectra, respectively. Default is
#' \code{"acceleration"}.
#' @param pct Desired percentile for results. Default is 50 (median).
#' @param ptap Percentage of data window to apply a Hanning taper. Must be
#' between 0 and 50. Default is 0.
#' @param damp Damping value to use. Default is 0.05
#' @param tau.range Vector of length 2 defining the range of response
#' spectral periods, in seconds, to consider. Default is \code{c(0,10)}.
#' @param periods Vector of response spectral periods, in seconds, to
#' use. If not defined, a set of periods in the range tau.range is generated. The
#' periods generated are a superset of the 2008 NGA periods.
#' @param rs.meth Numeric value indicating the method to use for computing response
#' spectra. Currently implemented methods are: \tabular{rl}{
#'  0 \tab recursive IIR (default) \cr
#'  1 \tab bilinear Z-transform \cr
#'  2 \tab original method of Nigam and Jennings (1969) \cr
#'  3 \tab circular convolution (slow) \cr
#'  4 \tab FFT method \cr
#' }
#' These methods generally give similar results, except for periods < ~5 * dt,
#' when the results may substantially differ.
#' @param tau.si.range Period range to compute spectrum intensity. Default
#' is \code{c(0.1,0.5)} seconds for acceleration, \code{c(0.1,2.5)} seconds
#' for velocity, and \code{c(0.1,max_period)} seconds for displacement.
#'
#' @return A \code{list} with the complete orientation-independent response spectra
#' information, including:
#' response spectra (GMRotDpp, GMRotIpp, RotDpp, and RotIpp for the specified
#' percentile pp), peak amplitude and spectrum intensity for each of these, and
#' supporting info such as the input damping value, the input or default
#' response spectral periods, response-spectrum method, integration limits
#' for SI, and string labels for response spectra type, units, etc. The
#' specific list returned is:
#' \code{list(damping, periods, pct, GMpeakD, IMpeakD, GMpeakI, IMpeakI,
#' GMRotD, RotD, GMRotI, RotI, rs.type, rs.units, rs.method,
#' GMsiD, IMsiD, GMsiI, IMsiI, SI.per.range, SI.units,
#' IMpenalty, GMpenalty, GMrs.D, IMrs.D, GMAngleD, IMAngleD)}
#'
#' @details Boore et al (2006) and Boore (2010) define several
#' orientation-independent measures of ground motion, including GMRotD50,
#' GMRotI50, RotD50, and RotI50. Empirical relations between several of
#' these measures are provided by Boore and Kishida (2017). See the cited
#' references for a complete description of these measures.
#'
#' Ground motion measures such as peak motion or spectral acceleration generally
#' depend on the specific orientation of the recorded data components. To obtain
#' orientation-independent measures, the original components are mathematically
#' rotated through a range of non-repeating angles, and the ground measure is
#' computed at each of these angles. Assuming a uniform distribution of the angles,
#' an empirical distribution of the resulting ground measures is obtained. The
#' result for a specified percentile from this distribution is then used as the
#' orientation-independent measure. The 50th percentile (median) measure
#' commonly is used. In contrast, the 0th and 100th percentile values provide
#' minimum and maximum values.
#'
#' The "GM" measures (e.g., GMRotI50) are based on the geometric mean of the two
#' horizontal components, whereas the non-GM measures (e.g., RotI50 or IMpeakI50) are
#' based on combining the two horizontal components into a single time series that
#' depends on the rotation angle.
#'
#' For response spectra, the computed percentile-values generally will depend on
#' the response spectral period. The computed "D" measures (e.g., GMRotD50) are
#' these period-dependent measures. It may be useful to obtain period-independent
#' measures. This is done by finding the rotation angle that minimizes the spread
#' of a "D" measure. The "I" measures (e.g., GMRotI50 or RotI50) are period-independent.
#' @seealso
#' \itemize{
#' \item  \code{\link{rspec}}
#' \item  \href{https://pubs.geoscienceworld.org/ssa/bssa/article/59/2/909/116747/calculation-of-response-spectra-from-strong-motion}{Nigam and Jennings (1969)}
#' Nigam, N.C., and Jennings, P.C., 1969, Calculation of
#' response spectra from strong-motion earthquake records: Bulletin of the Seismological
#' Society of America, v. 59, no. 2, p. 909-922.
#' \item  \href{https://doi.org/10.1785/0120050209}{Boore et al (2006)} Boore,
#' D.M., Watson, L.J., and Abrahamson, N.A., 2006, Orientation-independent measures
#' of ground motion: Bulletin of the Seismological Society of America, v. 96, no. 4a,
#' p. 1502-1511.
#' \item  \href{https://doi.org/10.1785/0120090400}{Boore (2010)} Boore, D.M., 2010,
#' Orientation-Independent, Nongeometric-Mean Measures of Seismic Intensity from
#' Two Horizontal Components of Motion: Bulletin of the Seismological Society of
#' America, v. 100, no. 4, p. 1830-1835.
#' \item  \href{https://doi.org/10.1785/0120160250}{Boore and Kishida (2017)} Boore,
#' D.M., and Kishida, T., 2017, Relations between Some Horizontal‐Component
#' Ground‐Motion Intensity Measures Used in Practice: Bulletin of the Seismological
#' Society of America, v. 107, no. 1, p. 334.
#' }
#' @keywords ts

IMRot.default <- function(xt, yt, dt=NA, units.ts=NA, ts.type="vel", rs.type="acc",
                   pct=50, ptap=NA, damp=NA, tau.range=NA, periods=NA,
                   rs.meth=NA, tau.si.range=NA) {
  if ( missing(xt) )
    stop("Must provide input xt")

  if ( is.na(dt) )
    dt = 0.01

  pct <- min(max(0,pct),100)

  # this method can use only one damping value, unlike rspec
  if ( ! is.na(damp) && length(damp) > 1 )
    damp <- damp[1]

  multi.trace <- is.matrix(xt) || is.mts(xt) ||
      ( ! is.null(dim(xt)) && length(dim(xt)) > 1 && dim(xt)[2] > 1 )

  # split out the X and Y components
  if ( multi.trace ) {
    if ( dim(xt)[2] < 2 )
      stop("multi-trace time series must have at least 2 components")
    cnames <- colnames(xt@data)
    have.SEED.traces = FALSE
    if ( ! is.null(cnames) ) {
      # try to get traces assuming SEED naming convention
      x.ind <- grep("E$", cnames)
      y.ind <- grep("N$", cnames)
      if ( length(x.ind) > 0 && length(y.ind) > 0 ) {
        X <- xt[,x.ind]
        Y <- xt[,y.ind]
        have.SEED.traces = TRUE
      }
    }
    if ( ! have.SEED.traces ) {
      # no matching SEED names, so just use the first 2 traces
      X <- xt[,1]
      Y <- xt[,2]
    }
  } else {
    if ( missing(yt) )
      stop("Must provide input xt and yt")
    else if ( length(xt) != length(yt) )
      stop("xt and yt must be the same length")
    X <- as.vector(xt)
    Y <- as.vector(yt)
  }

  # get a set of non-redundant angles (0-90 deg). For those IMs that
  # need 0-180 deg, we'll just add 90 deg to a separate case in the 0-90 loop
  ddeg <- 1
  angles <- seq(from=0, to=(90 - ddeg), by=ddeg )

  ######### get the period-dependent results
  periods <- NULL
  GMrs.D <- NULL
  GMpeak.D <- NULL
  GMsi.D <- NULL
  IMrs.D <- NULL
  IMpeak.D <- NULL
  IMsi.D <- NULL
  deg2rad <- pi / 180.
  angles.names <- NULL
  angles1.names <- NULL
  for ( phi in angles ) {
    angles.names <- c(angles.names, sprintf("%.1f deg",phi))
    angles1.names <- c(angles1.names, sprintf("%.1f deg",phi),
                       sprintf("%.1f deg",phi + 90))

    # get the rotated time histories
    RX <- X * cos(phi * deg2rad) + Y * sin(phi * deg2rad) # 0-90 deg
    RY <- -X * sin(phi * deg2rad) + Y * cos(phi * deg2rad) # 0-90 deg
    RX.1 <- X * cos((phi + 90) * deg2rad) + Y * sin((phi + 90) * deg2rad) # 90-180 deg

    # add the peak values to list
    GMpeak.D <- c(GMpeak.D, sqrt(max(abs(RX)) * max(abs(RY)))) # 0-90 deg
    # GMpeakD <- c(GMpeakD, sqrt(max(abs(RX * RY)))) # 0-90 deg
    IMpeak.D <- c(IMpeak.D, max(abs(RX)), max(abs(RX.1))) # 0-90, 90-180 deg

    # get the response spectra of the rotated time histories
    rs.x <- rspec(RX, dt=dt, units.ts=units.ts, ts.type=ts.type, rs.type=rs.type,
                  ptap=ptap, damp=damp, tau.range=tau.range, periods=periods,
                  rs.meth=rs.meth, tau.si.range=tau.si.range) # 0-90 deg
    rs.y <- rspec(RY, dt=dt, units.ts=units.ts, ts.type=ts.type, rs.type=rs.type,
                  ptap=ptap, damp=damp, tau.range=tau.range, periods=periods,
                  rs.meth=rs.meth, tau.si.range=tau.si.range) # 0-90 deg
    rs.x1 <- rspec(RX.1, dt=dt, units.ts=units.ts, ts.type=ts.type, rs.type=rs.type,
                   ptap=ptap, damp=damp, tau.range=tau.range, periods=periods,
                   rs.meth=rs.meth, tau.si.range=tau.si.range) # 90-180 deg

    # add response spectra for rotated time histories, by column, to the matrix
    GMrs.D <- cbind(GMrs.D, sqrt(rs.x$rspect * rs.y$rspect)) # 0-90 deg
    IMrs.D <- cbind(IMrs.D, rs.x$rspect, rs.x1$rspect) # 0-90, 90-180 deg

    # add spectrum intensities to list
    GMsi.D <- c(GMsi.D, sqrt(rs.x$SI * rs.y$SI)) # 0-90 deg
    IMsi.D <- c(IMsi.D, rs.x$SI, rs.x1$SI) # 0-90, 90-180 deg

    # get list of response periods and other info the first time through
    if ( is.null(periods) ) {
      rs.damping <- rs.x$damping
      rs.periods <- rs.x$periods
      rs.type.out <- rs.x$rs.type
      rs.units <- rs.x$rs.units
      rs.method <- rs.x$rs.method
      SI.per.range <- rs.x$SI.per.range
      SI.units <- rs.x$SI.units
    }
  }
  colnames(GMrs.D) <- angles.names
  colnames(IMrs.D) <- angles1.names

  # get the period-dependent quantile results
  prob <- 0.01 * pct
  GMpeakD <- quantile(GMpeak.D, prob, names=FALSE)
  IMpeakD <- quantile(IMpeak.D, prob, names=FALSE)
  GMsiD <- quantile(GMsi.D, prob, names=FALSE)
  IMsiD <- quantile(IMsi.D, prob, names=FALSE)
  GMRotD <- NULL
  IMRotD <- NULL
  for ( ii in 1:length(rs.periods) ) {
    # for each response period, get the percentile value over all angles
    GMRotD <- c(GMRotD, quantile(GMrs.D[ii,], prob, names=FALSE))
    IMRotD <- c(IMRotD, quantile(IMrs.D[ii,], prob, names=FALSE))
  }

  # remove periods with zero-spectra (velocity and displacement spectra)
  nonzero.ind <- which(GMRotD > 0)
  if ( length(nonzero.ind) < length(rs.periods) ) {
    rs.periods <- rs.periods[nonzero.ind]
    GMrs.D <- GMrs.D[nonzero.ind,]
    IMrs.D <- IMrs.D[nonzero.ind,]
    GMRotD <- GMRotD[nonzero.ind]
    IMRotD <- IMRotD[nonzero.ind]
  }

  # get the angles for the percentile values of the response spectra, by period
  GMAngleD <- NULL
  IMAngleD <- NULL
  for ( ii in 1:length(rs.periods) ) {
    # for each response period, get the percentile value over all angles
    eps <- 1e-2
    ang.ind <- which.min( abs(GMrs.D[ii,] - GMRotD[ii]) )
    GMAngleD <- c(GMAngleD, angles[ang.ind])
    # for IMs over 0-180 deg, we computed alternating pairs of values at phi,
    # and phi + 90, so to recover the actual angle from the index we can use
    # the fact that odd indexes are for phi, and even indexes are for phi + 90
    IM.ang.ind <- which.min( abs(IMrs.D[ii,] - IMRotD[ii]) )
    ang.ind <- (IM.ang.ind + 1) %/% 2
    phi <- angles[ang.ind]
    if ( IM.ang.ind %% 2 == 0 )
      phi <- phi + 90
    IMAngleD <- c(IMAngleD, phi)
  }


  ######### get the period-independent results

  # compute the penalty functions
  GMpenalty <- NULL
  IMpenalty <- NULL
  for ( jj in 1:length(angles) ) {
    # for each non-redundant angle, normalize the rotated IMs
    # by the percentile result, and compute the penalty function

    # Note: these GMs go from 0-90
    GMrs.norm <- GMrs.D[,jj] / GMRotD - 1 # vector over periods, for this phi
    GMpenalty.phi <- sum( GMrs.norm * GMrs.norm ) / length(rs.periods)
    GMpenalty <- c(GMpenalty, GMpenalty.phi) # cat the result for this angle

    # Note: these IMs go from 0-180, which we do in pairs at phi and phi + 90
    IMrs.norm <- IMrs.D[,(2 * jj - 1)] / IMRotD - 1 # vector at phi, over periods
    IMpenalty.phi <- sum( IMrs.norm * IMrs.norm ) / length(rs.periods)
    IMrs.norm.1 <- IMrs.D[,(2 * jj)] / IMRotD - 1 # vector at phi + 90, over periods
    IMpenalty.phi.1 <- sum( IMrs.norm.1 * IMrs.norm.1 ) / length(rs.periods)
    IMpenalty <- c(IMpenalty, IMpenalty.phi, IMpenalty.phi.1) # cat the pair of results
  }

  # find the angles (and indices) which minimize the penalty functions
  GMmin.phi.ind <- which.min(GMpenalty)
  GMmin.phi <- angles[GMmin.phi.ind]
  IMmin.phi.ind <- which.min(IMpenalty)
  ang.ind <- (IMmin.phi.ind + 1) %/% 2
  IMmin.phi <- angles[ang.ind]
  if ( IMmin.phi.ind %% 2 == 0 )
    IMmin.phi <- IMmin.phi + 90

  # get the peak values at the minimum angle
  GMpeakI <- GMpeak.D[GMmin.phi.ind]
  IMpeakI <- IMpeak.D[IMmin.phi.ind]

  # get the response spectra at the minimum angle
  GMRotI <- GMrs.D[,GMmin.phi.ind]
  IMRotI <- IMrs.D[,IMmin.phi.ind]

  # get the spectrum intensities at the minimum angle
  GMsiI <- GMsi.D[GMmin.phi.ind]
  IMsiI <- GMsi.D[IMmin.phi.ind]

  # reorder IM results for output
  ind.180 <- c(seq(1, by=2, length.out=length(angles)),
               seq(2, by=2, length.out=length(angles)))
  IMpenalty <- IMpenalty[ind.180]
  IMrs.D <- IMrs.D[,ind.180]

  ret <- list(damping=rs.damping, periods=rs.periods, pct=pct,
              GMpeakD=GMpeakD, IMpeakD=IMpeakD,
              GMRotD=GMRotD, RotD=IMRotD,
              GMmin.phi=GMmin.phi, IMmin.phi=IMmin.phi,
              GMpeakI=GMpeakI, IMpeakI=IMpeakI,
              GMRotI=GMRotI, RotI=IMRotI,
              rs.type=rs.type.out, rs.units=rs.units, rs.method=rs.method,
              GMsiD=GMsiD, IMsiD=IMsiD, GMsiI=GMsiI, IMsiI=IMsiI,
              SI.per.range=SI.per.range, SI.units=SI.units,
              GMpenalty=GMpenalty, IMpenalty=IMpenalty,
              GMrs.D=GMrs.D, IMrs.D=IMrs.D,
              GMAngleD=GMAngleD, IMAngleD=IMAngleD)

}
setGeneric("IMRot",def=IMRot.default)


#' @describeIn IMRot.default computes orientation-independent values for \code{ts}
#' univariate or \code{mts} multivariate time series
IMRot.ts <- function(xt, yt, dt=NA, units.ts=NA, ts.type="vel", rs.type="acc",
                      pct=50, ptap=NA, damp=NA, tau.range=NA, periods=NA,
                      rs.meth=NA, tau.si.range=NA) {
  if ( missing(xt) )
    stop("Must provide input xt")

  multi.trace <- is.mts(xt)
  dt <- deltat(xt)
  if ( is.na(dt) )
    dt <- 0.01

  # split out the components
  if ( multi.trace ) {
    if ( dim(xt)[2] < 2 )
      stop("multi-trace time series must have at least 2 components")
    cnames <- colnames(xt@data)
    have.SEED.traces = FALSE
    if ( ! is.null(cnames) ) {
      # try to get traces assuming SEED naming convention
      x.ind <- grep("E$", cnames)
      y.ind <- grep("N$", cnames)
      if ( length(x.ind) > 0 && length(y.ind) > 0 ) {
        X <- xt[,x.ind]
        Y <- xt[,y.ind]
        have.SEED.traces = TRUE
      }
    }
    if ( ! have.SEED.traces ) {
      # no matching SEED names, so just use the first 2 traces
      X <- xt[,1]
      Y <- xt[,2]
    }
  } else {
    if ( missing(yt) )
      stop("Must provide input xt and yt")
    else if ( length(xt) != length(yt) )
      stop("xt and yt must be the same length")
    X <- xt
    Y <- yt
  }

  # call default method
  ret <- IMRot.default(X, Y, dt=dt, units.ts=units.ts, ts.type=ts.type,
                        rs.type=rs.type, pct=pct, ptap=ptap, damp=damp,
                        tau.range=tau.range, periods=periods,
                        rs.meth=rs.meth, tau.si.range=tau.si.range)
  return(ret)
}
setMethod("IMRot","ts",IMRot.ts)


#' @describeIn IMRot.default computes orientation-independent
#' values for a \code{signalSeries} univariate or multivariate time series
IMRot.signalSeries <- function(xt, yt, dt=NA, units.ts=NA, ts.type="vel", rs.type="acc",
                      pct=50, ptap=NA, damp=NA, tau.range=NA, periods=NA,
                      rs.meth=NA, tau.si.range=NA) {
  if ( missing(xt) )
    stop("Must provide input xt")

  multi.trace <- ! is.null(dim(xt))
  if ( is.na(dt) )
    dt <- deltat(xt)

  # get units from first time series
  units <- xt@units
  if ( is.na(units) || is.null(units) ) {
    units <- NA
  } else if ( units == "cm/s/s" || units == "cm/s^2" ||
              units == "cm/s" || units == "cm" ) {
    units <- "cgs"
  } else if ( units == "m/s/s" || units == "m/s^2" ||
              units == "m/s" || units == "m" ) {
    units <- "SI"
  } else
    units <- NA

  # split out the X and Y components
  if ( multi.trace ) {
    if ( dim(xt)[2] < 2 )
      stop("multi-trace time series must have at least 2 components")
    cnames <- colnames(xt@data)
    have.SEED.traces = FALSE
    if ( ! is.null(cnames) ) {
      # try to get traces assuming SEED naming convention
      x.ind <- grep("E$", cnames)
      y.ind <- grep("N$", cnames)
      if ( length(x.ind) > 0 && length(y.ind) > 0 ) {
        X <- xt@data[,x.ind]
        Y <- xt@data[,y.ind]
        have.SEED.traces = TRUE
      }
    }
    if ( ! have.SEED.traces ) {
      # no matching SEED names, so just use the first 2 traces
      X <- xt@data[,1]
      Y <- xt@data[,2]
    }
  } else {
    if ( missing(yt) )
      stop("Must provide input xt and yt")
    else if ( length(xt) != length(yt) )
      stop("xt and yt must be the same length")
    X <- xt@data
    Y <- yt@data
  }

  # call default method
  ret <- IMRot.default(X, Y, dt=dt, units.ts=units.ts, ts.type=ts.type,
                        rs.type=rs.type, pct=pct, ptap=ptap, damp=damp,
                        tau.range=tau.range, periods=periods,
                        rs.meth=rs.meth, tau.si.range=tau.si.range)
  return(ret)
}
setMethod("IMRot","signalSeries",IMRot.signalSeries)
