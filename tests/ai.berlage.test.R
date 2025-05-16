# test berlage and AI
library("splus2R")
library(seismic.ts)

amp <- 1500 # cm/s^2
dt <- .01 # sec
f0 <- 0.2 # Hz
dur <- 20 # sec
syn.ber <- ts(c(rep(0, 0.05 * dur/dt), berlage(amp,f0,dur,dt)),
              start=0, deltat=dt)
tsplot(syn.ber)
ai.ber <- ai(syn.ber,dt,in.units="cgs")
sprintf("Arias Intensity = %.5f %s, Arias Duration = %.3f %s",
        ai.ber$AI, ai.ber$AI.units, ai.ber$dur, ai.ber$dur.units)
tsplot(ts(ai.ber$ai.integral,start=0,deltat=dt))
