library(fda)
library(refund)
mortality <- read.csv("/Users/josephpiekos/Desktop/project III/data/english_mortality_final.csv")

counties <- names(mortality) #contains a list of all counties
#The index in this vector is what should be put into smoothing functions for "start.region" and "end.region", to choose the number of mortality curves that will be plotted.


###################################################################################################################################################################################################

#mortality.plot() plots the mortality curves for a subset of the 78 counties we are considering, from day 0 to any day up to day 180.
#day 0 = 15/02/20
#day 180 = 13/08/20

mortality.plot <- function(start.region, end.region, total.days, dataset=mortality, nbasis = 180, norder = 4,lambda=100){

  number.regions = (end.region - start.region) + 1 #number of replicates i.e. number of regions to plot on top of each other
  mortality.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions
  
  for (i in 1:total.days){ #creates a matrix containing the data which can be used with fda commands
    for (j in start.region:end.region){
      if(number.regions > 1){
        mortality.matrix[i,j-(start.region-1)] <- dataset[i,j]
      }
      else {
        mortality.matrix[i,1] <- dataset[i,j]
      }
    }
  }
  
  mortalitybasis = create.bspline.basis(c(0,total.days-1), nbasis, norder)
  mortalityfdPar = fdPar(mortalitybasis,2,lambda) 
  mortalityfd = smooth.basis((0:(total.days-1)),mortality.matrix,mortalityfdPar)$fd
  plot(mortalityfd,xlab="Day number",ylab="Mortality (Deaths per 100,000)",ylim=c(-0.01,2.2)) #max(eval.fd((0:(total.days-1)),mortalityfd)))
  if (number.regions > 1){
    title(paste("Mortality curves for",number.regions, "counties"))
    lines(mortalityfd, lty=1, col="black")
    lines(mean.fd(mortalityfd),lty = 2,lwd=3,col="orange")
  }
  else{
    lines(mortalityfd, lty=1, col = "blue")
  }
  abline(v = 37,lty = 5, col = "red") #adds red vertical line at day 37
  axis(1, at=37,labels=37,cex.axis = 0.7, tck = -0.06) #label on axis at x=37
  abline(v = 30,lty = 3, col = "orange") #adds orange vertical line at day 37
  axis(1, at=30,labels=30, cex.axis = 0.7, tck = -0.06) #label on axis at x = 30
  legend(85, 2.2, legend = c("No non-essential contact and travel","First Lockdown"), col = c("orange","red"), lty = c(3,5), cex = 0.65, text.font = 2, title = "Key Announcements from Government")
}

mortality.plot(start.region = 7,end.region = 8, total.days = 180, dataset=mortality,nbasis = 180 ,norder = 4,lambda=800) #lambda is the smoothing parameter and has been selected using the GCV methods (see GCV code file)


################################################################################################################################################################################################################################################################################################

#mortality.pos smooths mortality with the constraint that it is positive everywhere. (constrained function as described in chapter 4.4 of my report)
#it is an iterative method, so may take a bit of time if the full 78 curves are desired.
mortality.pos <- function(start.region,end.region,total.days, dataset, nbasis = 19, norder = 4,lambda=100){
  
  number.regions = (end.region - start.region) + 1 #number of replicates i.e. number of regions to plot on top of each other
  mortality.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions
  
  for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
    for (j in start.region:end.region){ #there are 3 columns, each column j is a region
      if(number.regions > 1){
        mortality.matrix[i,j-(start.region-1)] <- dataset[i,j]
      }
      else {
        mortality.matrix[i,1] <- dataset[i,j]
      }
    }
  }
  
  mortalitybasis = create.bspline.basis(c(0,total.days-1), nbasis, norder)
  mortalityfdPar = fdPar(mortalitybasis,2,lambda) 
  mortalityfd = smooth.pos((0:(total.days-1)),as.matrix(mortality.matrix[,1]),mortalityfdPar) #smooth.pos is a function from the fda package which returns the unconstrained smooth logarithm function W
  Wfd = mortalityfd$Wfdobj #the unconstrained smoothed logarithm function we require for positive constrained smoothing.
  
  posfit = matrix(0,180,number.regions) #creates an empty matrix, where each column will contain the final positive function, by evaluating our W function at a fine mesh of points
  posfit[,1] = exp(eval.fd((0:(total.days-1)), Wfd)) 
  plot(1:total.days,mortality[1:total.days,1],col="white",pch=20,ylim=c(-0.01,2.2),xlab="Day number",ylab="Mortality (Deaths per 100,000)")
  lines((0:(total.days-1)), posfit[,1],lwd=2) 
  
  for(k in 2:number.regions){ #this for loop creates the rest of the curves other than the first
  mortalityfd = smooth.pos((0:(total.days-1)),as.matrix(mortality.matrix[,k]),mortalityfdPar)
  Wfd = mortalityfd$Wfdobj
  mortality.matrix
  posfit[,k] = exp(eval.fd((0:(total.days-1)), Wfd))
  lines((0:(total.days-1)), posfit[,k],lwd=2) 
  }
  #up to this point the function will simply plot the now positive constrained functions
  #however if we want to perform further fda with the positive functions, we need to return them as a functional data object out of the function
  mortalitybasis = create.bspline.basis(c(0,total.days), 100, 4) #sets up a basis to resmooth the values of the now positive mortality function
  mortalityfdPar = fdPar(mortalitybasis,2,0.001) #a functional parameter object, containing the derivative to penalise in the roughness penalty and smoothing parameter to use.
  mortality.posfd = smooth.basis((0:(total.days-1)),posfit,mortalityfdPar)$fd 
  lines(mean.fd(mortality.posfd),lty=2,col="orange",lwd=2) #adds the mean of the positive functions to the plot, in orange.
  abline(v = 37,lty = 5, col = "red")
  axis(1, at=37,labels=37,cex.axis = 1.2, tck = -0.03) #label on axis at x=37
  abline(v = 30,lty = 3, col = "orange")
  axis(1, at=30,labels=30, cex.axis = 1.2, tck = -0.03) #label on axis at x = 30
  abline(h = 0,lty = 5, col = "black")
  legend(85, 2.2, legend = c("No non-essential contact and travel","First Lockdown"), col = c("orange","red"), lty = c(3,5), cex = 0.75, text.font = 2, title = "Key Announcements from Government")
  return(mortality.posfd) #returns a functional data object for the positive mortality functions that can be used in other FDA methods.
}
mortality.posfd = mortality.pos(1,3,180,mortality,180,4,5) #runs the positive smoothing function and stores the resulting set of curves in a variable

####################
#creates Figure 1.2, and saves it as a jpeg in R's current working directory
jpeg("10 nonnegative mortality curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=2,cex.axis=2,mar=c(5,6,2,1)+.1)
mortality.pos(1,10,180,mortality,180,4,5) #if using mortality.regions, make sure not to choose more than 9 regions
dev.off()
###################
jpeg("aligned vs not aligned 20 mortality curves.jpeg", units = "in", width = 10, height = 4 , res = 400)
par(cex.axis=1.3,cex.lab=1.5,mfrow=c(1,2),mar = c(4,2.5,1,0.5))
mortality.plot(start.region = 1,end.region = 20, total.days = 180, dataset=mortality,nbasis = 90 ,norder = 4,lambda=794)
continuous.reg(1,20,180)
#lines(mean.fd(accelfdCR),lty=2, col = "orange",lwd = 4)
dev.off()


jpeg("50 mortality curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)
mortality.plot(1,50,90,19,4)
dev.off()

