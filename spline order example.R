library(stats)
library(fda)
#plot(splinefun(c(0,pi/2,pi,(3*pi)/2,2*pi),c(0,1,0,-1,0))) 
basis2 = create.bspline.basis(c(0,2*pi), 5, 2)
basis3 = create.bspline.basis(c(0,2*pi), 6, 3)
basis4 = create.bspline.basis(c(0,2*pi), 7, 4)
splinefd = smooth.basis(c(0,pi/2,pi,(3*pi)/2,2*pi),c(0,1,0,-1,0),basis4)$fd
par(mfrow=c(3,1))
plot(splinefd)
create.bspline.basis()

n = 5
argvals = seq(0,(2*pi),len=n)
#  The true curve values are sine function values with period 1/2
x = sin(argvals)
#  Add independent Gaussian errors with std. dev. 0.2 to the true values
sigerr = 0.05
y = x 
#+ rnorm(x)*sigerr
#  When we ran this code, we got these values of y (rounded to two
#  decimals):
#  Set up a B-spline basis system of order 4 (piecewise cubic) and with
#  knots at 0, 0.1, ..., 0.9 and 1.0, and plot the basis functions
norder = 2
nbasis = n

basisobj = create.bspline.basis(c(0,2*pi),nbasis,norder)

ys = smooth.basis(argvals=argvals, y=y, fdParobj=basisobj)$fd

norder = 3
basisobj2 = create.bspline.basis(c(0,2*pi),nbasis,norder)
ys2 = smooth.basis(argvals=argvals, y=y, fdParobj=basisobj2)$fd

jpeg("spline order 2 vs order 3.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.3,cex.lab=1.5,cex.pch=2,mfrow=c(1,2),mar = c(4,2.5,1,0.5))
plot(ys,ylim=c(-1,1),xlab="t",ylab="",yaxt="n")
points(argvals,y)
plot(ys2,ylim=c(-1,1),xlab="t",ylab="",yaxt="n")
points(argvals,y)
dev.off()
