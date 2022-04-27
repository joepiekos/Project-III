library(fda)

#this code creates figure 6.1, the covariance function of the 78 mortality curves.
#it also plots the contour plot associated with the surface to aid interpretation.

jpeg("mortality covariance.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.2,cex.lab=1.5,mfrow=c(1,2),mar = c(4.5,4.5,2,0))

mortality.covfd = var.fd(mortality.posfd) 
#this function above takes the set of smoothed and positively contrained mortality curves. 
#These are computed using mortality.pos in the "mortality smoother" file, and just make sure this is set to smooth all 78 curves before computing thecovariance function here
day.seq = seq(30,100,length.out = 50)
mortality.covfunc = eval.bifd(day.seq,day.seq,mortality.covfd)
contour(day.seq,day.seq,mortality.covfunc,xlab="Day",ylab="Day")
title("Contour plot of surface")
par(cex.axis=0.7,cex.lab=1) 
persp(day.seq,day.seq, mortality.covfunc,
      theta=-25, phi=20, r=3, expand = 0.5,
      ticktype="detailed",
      xlab="Day",
      ylab="Day",
      zlab="",col="orange", shade = 0.2,border="black")
title("Mortality covariance (days 30 to 100)")
dev.off()