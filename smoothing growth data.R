library(fda)
heightmat <- growth$hgtm[,1:20]
age <- growth$age
heightbasis = create.bspline.basis(c(1,18),35,6,age)
heightfdPar = fdPar(heightbasis, 2, 0.01)
heightfd = smooth.basis(age, heightmat,
                        heightfdPar)

jpeg("10 boys growth curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)

plot(heightfd,xlab="Age (years)",ylab="Height (cm)")
lines(heightfd,lty=1,col="black")
for(j in 1:20){
  points(age,heightmat[,j],col="black",pch=20)
}
dev.off()

