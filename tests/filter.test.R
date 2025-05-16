# test filter.llnl
library("splus2R")
library(seismic.ts)

# make a delta function, centered in a window of 2000 points
dirac <- c(rep(0,1000),1,rep(0,1000))

# choose time increment, and high- and low-cut frequencies
dt <- 0.01
f.lo <- 20 / (length(dirac) * dt)
f.hi <-  (2 / 5) * 1 / (2 * dt)

# filter the delta function
response <- filter.llnl(dirac, dt, order=8, pb.type="bp", filt.type="c2",
                           f.lo=f.lo, f.hi=f.hi, dir="zp", cheb.sb.atten=100,
                           cheb.tr.bw=0.4)

# convert filtered delta function to a ts object, and plot
response <- ts(response, start=0, deltat=dt)
tsplot(response)

# compute the power spectrum, and plot
filt.spec <- spec.pgram(response,plot=FALSE)
ymax <- 1.2 * 10 ^ round(log10(max(filt.spec$spec)),digits=0)
ymin <- 0.8 * 10 ^ (round(log10(max(filt.spec$spec)),digits=0) - 10)
plot(filt.spec,ylab="power spectrum",xlab="frequency, hz",ylim=c(ymin,ymax))
