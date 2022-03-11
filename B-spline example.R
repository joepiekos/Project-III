library(fda)
splinebasis <- create.bspline.basis(norder=3,c(0,5),breaks=c(0,1,2,3,4,5))
plot(splinebasis,lty=1) #plots full B-spline basis function

piecewise <- function(x){
  f <- NULL
  f[x<=1 & x>=4] <- 0*(x[x<=1 & x>4])
  f[x>=1 & x<2] <- (1/2)*x[x>=1 & x<2]^2 - x[x>=1 & x<2] + (1/2)
  f[x>=2 & x<3] <- (-1)*x[x>=2 & x<3]^2 + 5*x[x>=2 & x<3] - (11/2)
  f[x>=3 & x<4] <- (1/2)*x[x>=3 & x<4]^2 - 4*x[x>=1 & x<2]
  f
} #defines the piecewise B-spline B_{1,3}(t) which was computed using the recursive Cox-de Boor formula

x <- seq(0, 5, 0.01)
plot(x, piecewise(x), type="l", las=1,ylim=c(0,1),ylab="",xlab="t")

jpeg("order 3 B-spline basis.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.3,cex.lab=1.5,cex.pch=2,mfrow=c(1,2),mar = c(4,2.5,1,0.5))
plot(splinebasis,lty=1,xlab="t") #plots full B-spline basis function
plot(x, piecewise(x), type="l", las=1,ylim=c(0,1),ylab="",xlab="t")
dev.off()
