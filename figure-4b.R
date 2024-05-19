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


parms=c(10, 0.1, 3.75, 20, 1)
lim=25
par(pty="s")
plot(c(0,0), type="n", xlim=c(0,lim), ylim=c(0,lim), 
     ylab="scaled pollinator population density", 
     xlab="scaled plant population density", main="Plant and pollinator densities nullclines")
u=flowField(deriv=scaledModelmy, points=20, xlim=c(0,lim), ylim=c(0,lim), col="lightblue", parameters=c(parms, 1))
nullclines = nullclines(scaledModelmy, xlim = c(0, lim), ylim = c(0, lim), parameters = c(parms,1), points = 250, add.legend=FALSE )
