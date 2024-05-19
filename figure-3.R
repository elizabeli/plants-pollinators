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



addtraj=function(parms, y0){
  integration_time=200
  t=trajectory(scaledModel, tlim = c(0,integration_time),
               parameters=c(parms,1), y0, lwd=2)
  t=trajectory(scaledModel, tlim = c(0,integration_time),
               y0, parameters=c(parms,-1), lwd=2) 
}


  parms=c(10, 0.1, 3.75, 20, 1) #change parms here
  lim=25
  par(pty="s")
  plot(c(0,0), type="n", xlim=c(0,lim), ylim=c(0,lim), 
       ylab="scaled pollinator population density", 
       xlab="scaled plant population density", 
  )
  u=flowField(deriv=scaledModelmy, points=20, xlim=c(0,lim), ylim=c(0,lim), col="lightblue", parameters=c(parms, 1))
  nullclines = nullclines(scaledModelmy, xlim = c(0, lim), ylim = c(0, lim), parameters = c(parms,1), points = 250, add.legend=FALSE )
  mani=drawManifolds(scaledModelmy, y0=c(2.3, 2), parameters=c(parms, 1))
  
  x=c(rev(mani$stable.1[,2]), mani$stable.2[,2])
  x=head(x,-1)
  #y1 is upper
  y1=rep(lim+1, length(x))
  #y2 is y-values for the stable manifold
  y2=c(rev(mani$stable.1[,3]), mani$stable.2[,3])
  y2=head(y2,-1)
  
  #find where y1 and y2 intersect
  y2wonan=y2[!is.na(y2)]
  indexiny2wonan=which(abs(y2wonan-(lim+1))==min(abs(y2wonan-(lim+1))))
  val=y2[indexiny2wonan]
  indexforxmin=which(y2==val)
  plot(x, y1, type = "l",
       ylim = c(0,lim), xlim=c(0, lim),ylab = "scaled pollinator density", 
       xlab="scaled plant density", main=paste("parms=", parms[1], parms[2], parms[3],
                                               parms[4]))
  lines(x, y2, type = "l", col = 2)
  
  # Min and max X values
  xmin <- floor(x[indexforxmin])
  xmax <- 2*lim
  
  #add shaded region
  polygon(c(x[x >= xmin & x <= xmax],
            rev(x[x >= xmin & x <= xmax])),
          c(y2[x >= xmin & x <= xmax],
            rev(y1[x >= xmin & x <= xmax])),
          col = "#6BD7AF")
  e5=findEquilibrium(scaledModelmy,y0=c(lim, lim), parameters=c(parms,1), plot.it=TRUE)
  e5=findEquilibrium(scaledModelmy,y0=c(0.1, 0.1), parameters=c(parms,1), plot.it=TRUE)
  
