# library("splus2R", lib.loc="~/R/win-library/3.3")

# test creating seismic test pulses using 1/2 Ricker, or Berlage wavelet
mk.ricker = TRUE

# for Roman
vs = 7880 * .3048 # 7880 ft/s -> m/s
rho <- 150 * (1/2.204623) * (1/.3048)^3 # 150 lb/ft^3 -> kg/m^3
amp <- rho * vs * 0.01 # velocity to stress conversion for a 10 cm/s amplitude, Pascal
units <- "Pa"

f0 <- 3.1
dur <- 20.0
zpad <- min(0.1 * dur,1) # 10% of the duration, up to 1 second max
zpad <- 5 # for Roman
dt <- .01
f.lo <- min(5 / dur, f0/5)
f.hi <- 13

if ( mk.ricker ) {
  # test pulse based on part of a Ricker wavelet
  rk.len <- round(8 / ( pi * f0 * dt ))
  if ( rk.len %% 2 == 0 ) rk.len = rk.len + 1
  rk.hlen <- rk.len %/% 2
  rk.dur <- rk.len * dt
  # get central pulse half-width, per Wang (2015)
  rk.taulen <- round(0.88521 * rk.dur / (16 * dt))
  syn.rck <- ricker(amp,rk.dur,dt)[(rk.hlen+1-rk.taulen):rk.len]
  tsplot(syn.rck)
  readline(prompt = "hit Enter to continue")
  syn.rcklen <- length(syn.rck)
  # make minimum-phase, smooth tail with a hanning window
  #syn.w <- minimum.phase(c(syn.rck,rep(0,1000-length(syn.rck)), F))
  #tsplot(syn.w[(syn.rcklen+1):200])
  #readline(prompt = "hit Enter to continue")
  zp1 <- max(10, syn.rcklen %/% 5)
  pct <- 100 * zp1 / ( 2 * zp1 + syn.rcklen )
  zp2 <- round((100 / pct - 1) * zp1 - syn.rcklen + 0.5)
  s.len <- zp1 + syn.rcklen + zp2
  f.pad <- rep(0,zp1)
  syn.w <- hanning(c(f.pad,syn.rck)[1:s.len],pct,F)
  tsplot(syn.w)
  readline(prompt = "hit Enter to continue")
  # zero-pad out to the specified duration
  syn.w <- c(rep(0,zpad/dt-zp1),syn.w,rep(0,dur/dt - s.len - zpad/dt + zp1))
  tsplot(syn.w)
  readline(prompt = "hit Enter to continue")
} else {
  # test pulse based on a Berlage wavelet
  syn.ber <- berlage(amp,f0,dur-zpad,dt)
  tsplot(syn.ber)
  readline(prompt = "hit Enter to continue")
  # make minimum-phase
  syn.w <- minimum.phase(syn.ber,F)
  tsplot(syn.w)
  readline(prompt = "hit Enter to continue")
  syn.w <- c(rep(0,zpad/dt),syn.w)
  tsplot(syn.w)
  readline(prompt = "hit Enter to continue")
}
syn.w <- as.double(syn.w[!is.na(syn.w)])
# convert unfiltered pulse to a signal series
syn.m.u <- signalSeries(syn.w,from=0,by=dt,units=units,units.position="seconds")
plot(syn.m.u)
readline(prompt = "hit Enter to continue")

# baseline correct
filt.w <- filter.llnl(syn.w, dt, order=8, pb.type="bp", filt.type="be",
                      f.lo=f.lo, f.hi=f.hi, dir="Forward")
tsplot(filt.w[400:800])
readline(prompt = "hit Enter to continue")

# convert filtered pulse to a signal series
syn.m <- signalSeries(filt.w,from=0,by=dt,units=units,units.position="seconds")
plot(syn.m)
readline(prompt = "hit Enter to continue")

# plot amplitude spectrum of wavelet
spec.u <- spec.pgram(signalSeries2ts(syn.m.u),plot=FALSE)
ymax <- 1.2 * 10 ^ round(log10(max(spec.u$spec)),digits=0)
ymin <- 0.8 * 10 ^ (round(log10(max(spec.u$spec)),digits=0) - 8)
plot(spec.u,ylab="power spectrum",xlab="frequency, hz",ylim=c(ymin,ymax))
readline(prompt = "hit Enter to continue")
spec.f <- spec.pgram(signalSeries2ts(syn.m),plot=FALSE)
plot(spec.f,ylab="power spectrum",xlab="frequency, hz",ylim=c(ymin,ymax))
readline(prompt = "hit Enter to continue")

#getwd()
#exportData(syn.m, type="ASCII", delim=",", file="syn_wavelet_new_clean.csv")
#setwd("C:/Documents and Settings/cwood/My Documents/Studies/BF Sisk/Loma Prieta")
#write.table(syn.m, eol="\n", sep=",", file="syn_wavelet_new_clean1.csv", row.names=FALSE)

tt <- seq(from=dt,to=length(filt.w)*dt,by=dt)
tx <- data.frame(time_sec=tt,stress_Pa=filt.w)
setwd("C:/Users/CWood/Documents/Studies/Research/Monticello")
write.table(tx, eol="\n", sep=",", file="syn_rk_wavelet.csv", row.names=FALSE)

