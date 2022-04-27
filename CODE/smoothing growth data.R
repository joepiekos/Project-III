library(fda)
#growth data is already stored in the fda package, so does not need reading in like the other data sets.
#run the whole code in order to produce figure 1.1, and store it in a jpeg in R's current working directory

jpeg("smoothing height example.jpeg", units = "in", width = 7.5, height = 5, res = 400)
par(cex.axis=1.3,cex.lab=1.5,mfrow=c(1,1),mar = c(4.75,4.5,2,0.7))
heightmat <- growth$hgtm[,1:10] #takes the height data for the first 10 males and puts it into a matrix
age <- growth$age #puts the sampling points (ages between 1 and 18) into a vector

heightbasis = create.bspline.basis(c(1,18),35,6,age)
heightfdPar = fdPar(heightbasis, 2, 0.001)
heightfd = smooth.basis(age, heightmat,
                        heightfdPar)

plot(heightfd,xlab="Age (years)",ylab="Height (cm)",type="l",lty=1)
colours = c(1:10)
lines(heightfd,lty=1,col=colours)
dev.off()

