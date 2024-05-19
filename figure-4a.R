rm(list=ls(all=TRUE))
library(phaseR)
library(deSolve)
scaledModelmy=function(t, y, parms){
  #parms = k1, k2, k3, k4, b, direction
  u=y[1]; v=y[2]; 
  k1=parms[1]
  k2=parms[2]
  k3=parms[3]
  k4=parms[4]
  b=parms[5]
  direction=parms[6]
  du=(k3*v/(1+u+v)-1-k2*u)*b*u
  dv=(k4*u/(1+u+v)-k1)*b*v
  return(list(direction*c(du,dv)))
}

# FOR MU (k4)
k=20 #put k into the correct position & give baseline value
parms=c(10, 0.1, 3.75, k, 1)
lim=20 #set the plotting limit
eq0=findEquilibrium(scaledModelmy,y0=c(lim, lim),
                    parameters=c(parms,1), plot.it=FALSE)
eq_locations=matrix(c(eq0$ystar[1,], eq0$ystar[2,]), 1, 2)
step=0.2
for (i in 1:25) {
  k=k+step
  parms=c(10, 0.1, 3.75, k, 1)
  runningeq=findEquilibrium(scaledModelmy,y0=c(lim, lim),
                            parameters=c(parms,1), plot.it=FALSE)
  runningeqlocation=c(runningeq$ystar[1,], runningeq$ystar[2,])
  eq_locations=rbind(eq_locations, runningeqlocation)
}
mu=eq_locations

# FOR ETA (k3)
k=3.75 #put k into the correct position & give baseline value
parms=c(10, 0.1, k, 20, 1)
lim=20 #set the plotting limit
eq0=findEquilibrium(scaledModelmy,y0=c(lim, lim),
                    parameters=c(parms,1), plot.it=FALSE)
eq_locations=matrix(c(eq0$ystar[1,], eq0$ystar[2,]), 1, 2)
step=0.1
for (i in 1:50) {
  k=k+step
  parms=c(10, 0.1, k, 20, 1)
  runningeq=findEquilibrium(scaledModelmy,y0=c(lim, lim),
                            parameters=c(parms,1), plot.it=FALSE)
  runningeqlocation=c(runningeq$ystar[1,], runningeq$ystar[2,])
  eq_locations=rbind(eq_locations, runningeqlocation)
}
eta=eq_locations

# FOR C (k2)
k=0.1 #put k into the correct position & give baseline value
parms=c(10, k, 3.75, 20, 1)
lim=20 #set the plotting limit
eq0=findEquilibrium(scaledModelmy,y0=c(lim, lim),
                    parameters=c(parms,1), plot.it=FALSE)
eq_locations=matrix(c(eq0$ystar[1,], eq0$ystar[2,]), 1, 2)
step=0.01
for (i in 1:10) {
  k=k-step
  parms=c(10, k, 3.75, 20, 1)
  runningeq=findEquilibrium(scaledModelmy,y0=c(lim, lim),
                            parameters=c(parms,1), plot.it=FALSE)
  if (runningeq$classification == "Stable node") {
    runningeqlocation=c(runningeq$ystar[1,], runningeq$ystar[2,])
    eq_locations=rbind(eq_locations, runningeqlocation)
  } 
  else {
    runningeq=findEquilibrium(scaledModelmy,y0=c(2*lim, 2*lim),
                              parameters=c(parms,1), plot.it=FALSE)
    if (runningeq$classification == "Stable node") {
      runningeqlocation=c(runningeq$ystar[1,], runningeq$ystar[2,])
      eq_locations=rbind(eq_locations, runningeqlocation)
    } 
  }
}
c=eq_locations

# FOR D (k1)
k=10 #put k into the correct position & give baseline value
parms=c(k, 0.1, 3.75, 20, 1)
lim=20 #set the plotting limit
eq0=findEquilibrium(scaledModelmy,y0=c(lim, lim),
                    parameters=c(parms,1), plot.it=FALSE)
eq_locations=matrix(c(eq0$ystar[1,], eq0$ystar[2,]), 1, 2)
step=0.2
for (i in 1:40) {
  k=k-step
  parms=c(k, 0.1, 3.75, 20, 1)
  runningeq=findEquilibrium(scaledModelmy,y0=c(lim, lim),
                            parameters=c(parms,1), plot.it=FALSE)
  if (runningeq$classification == "Stable node" & runningeq$ystar[1,]>0) {
    runningeqlocation=c(runningeq$ystar[1,], runningeq$ystar[2,])
    eq_locations=rbind(eq_locations, runningeqlocation)
  } 
  else {
    runningeq=findEquilibrium(scaledModelmy,y0=c(100*lim, 100*lim),
                              parameters=c(parms,1), plot.it=FALSE)
    if (runningeq$classification == "Stable node"& runningeq$ystar[1,]>0) {
      runningeqlocation=c(runningeq$ystar[1,], runningeq$ystar[2,])
      eq_locations=rbind(eq_locations, runningeqlocation)
    } 
  }
}
d=eq_locations

plot(mu[,1], mu[,2], type="p",pch=17, cex=0.5, xlim=c(0,25), ylim=c(0,25), 
     ylab="scaled pollinator population density (v)", 
     xlab="scaled plant population density (u)", 
     main="Coexistence equilibrium locations for varying parameters", col="darkgreen")
lines(eta[,1], eta[,2], type="p",pch=17, cex=0.5, col="blue")
lines(c[,1], c[,2], type="p",pch=17, cex=0.5, col="red")
lines(d[,1], d[,2], type="p",pch=17, cex=0.5, col="orange")

