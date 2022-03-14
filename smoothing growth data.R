library(fda)
jpeg("smoothing vs interpolation height.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.3,cex.lab=1.5,cex.pch=2,mfrow=c(1,2),mar = c(4.75,2.5,2,0.7))
heightmat <- growth$hgtm[,1:5]
age <- growth$age
heightbasis = create.bspline.basis(c(1,18),35,6,age)
heightfdPar = fdPar(heightbasis, 2, 5)
heightfd = smooth.basis(age, heightmat,
                        heightfdPar)
create.bspline.basis()
# jpeg("10 boys growth curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
# par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)

plot(heightfd,xlab="Age (years)",ylab="Height (cm)")
lines(heightfd,lty=1,col="blue")
for(j in 1:5){
  points(age,heightmat[,j],col="black",pch=20)
}

dev.off()
#
jpeg("order 3 B-spline basis.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.3,cex.lab=1.5,cex.pch=2,mfrow=c(1,2),mar = c(4,2.5,1,0.5))

solve(crossprod(heightbasismat))

heightbasismat = eval.basis(age, heightbasis)
