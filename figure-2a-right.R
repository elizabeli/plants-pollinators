rm(list=ls(all=TRUE))

library(phaseR)
library(deSolve)

scaledModel=function(t,y,parms) {
  #parms=(k, b)
  u=y[1]; v=y[2]; 
  k=parms[1]
  b=parms[2]
  direction=parms[3]
  du=(3.75*v/(1+u+v)-1-0.1*u)*u*b
  dv=(k*u/(1+u+v)-10)*b*v
  return(list(direction*c(du, dv)))
}

addtraj=function(parms, y0){
  integration_time=200
  t=trajectory(scaledModel, tlim = c(0,integration_time),
               parameters=c(parms,1), y0, lwd=2)
  t=trajectory(scaledModel, tlim = c(0,integration_time),
               y0, parameters=c(parms,-1), lwd=2)
}
addtraj2=function(parms, y0, integration_time){
  t=trajectory(scaledModel, tlim = c(0,integration_time),
               parameters=c(parms,1), y0, lwd=2)
  t=trajectory(scaledModel, tlim = c(0,integration_time),
               y0, parameters=c(parms,-1), lwd=2)
}

parms=c(20, 0.1)
lim=7
par(pty="s")
plot(c(0,0), type="n", xlim=c(0,7), ylim=c(0,6), 
     ylab="scaled pollinator population density (v)", 
     xlab="scaled plant population density (u)")
u=flowField(deriv=scaledModel, points=20, xlim=c(0,lim), ylim=c(0,lim), col="lightblue", parameters=c(parms, 1))
#nclines = nullclines(scaledModel, xlim = c(0, lim), ylim = c(0, lim), 
#                     parameters = c(parms,1), points = 250, add.legend=FALSE )
e2=findEquilibrium(scaledModel,y0=c(4, 3), 
                   parameters=c(parms,1), plot.it=TRUE)
e2=findEquilibrium(scaledModel,y0=c(10, 13), 
                   parameters=c(parms,1), plot.it=TRUE)
e2=findEquilibrium(scaledModel,y0=c(0.1, 0.1), 
                   parameters=c(parms,1), plot.it=TRUE)
#drawManifolds(scaledModel, y0=c(4,3), parameters=c(parms, 1))
addtraj(parms, y0=c(1.7, 4))
addtraj(parms, y0=c(3, 1))
addtraj2(parms, y0=c(4.5, 1), integration_time=400)
addtraj(parms, y0=c(4.9, 1))
addtraj2(parms, y0=c(3.5, 5), integration_time=400)
addtraj2(parms, y0=c(3.2, 5), integration_time=400)
addtraj(parms, y0=c(5.3, 3))
addtraj(parms, y0=c(4.9, 4.9))
addtraj(parms, y0=c(5.7, 6))
addtraj(parms, y0=c(6.3, 4))
addtraj(parms, y0=c(7.4, 3.5))
addtraj(parms, y0=c(7.4, 5.7))
