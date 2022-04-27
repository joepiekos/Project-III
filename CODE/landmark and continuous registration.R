library(fda)
library(refund)
mortality <- read.csv("/Users/josephpiekos/Desktop/project III/data/english_mortality_final.csv")

########################################################################################################################
########################################################################################################################
####### this first block of code until line 56 sets up landmark registration by smoothing mortality derivative curves, then identifying 4 landmarks and storing the xcoordinates. ####
total.days = 179
nbasis = 180
norder = 4
lambda = 800

mortality.matrix <- (matrix(0,total.days,78)) #creates empty matrix, number of rows = time period, number of columns = number of regions

for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
  for (j in 1:78){ #there are 3 columns, each column j is a region
    if(number.regions > 1){
      mortality.matrix[i,j] <- mortality[i,j]
    }
    else {
      mortality.matrix[i,1] <- mortality[i,j]
    }
  }
}

mortalitybasis = create.bspline.basis(c(0,total.days), nbasis, norder)
mortalityfdPar = fdPar(mortalitybasis,2,lambda) 
mortalityfd = smooth.basis((1:(total.days-1)),mortality.matrix[1:178,1:2],mortalityfdPar)
mortalityaccfd = deriv.fd(mortalityfd$fd, 1)

landmark.i  = matrix(0,78,4) #creates an empty 78x4 matrix that will store the x coordinates of the 4 landmarks for the 78 curves.
#add columns to this matrix if you want to align more than four landmarks
days = seq(1,178,len=178) #a fine mesh of time argument values to evaluate the curves at when identifying the landmarks
par(mfrow=c(1,1), ask=TRUE) #graphical parameter to ensure the curves are plotted one by one and can be cycled through manually


for (ilandmark in 1:ncol(landmark.i)){
  for (icase in 1:78){
  mortalityveci = predict.fd(mortalityaccfd[icase], days)
  plot(days,mortalityveci,"l", ylim=c(-0.11,0.18),
       xlab="day", ylab="Mortality derivative",
       main=paste("County",icase))
  lines(c(0,180),c(0,0),lty=2)
  landmark.i[icase,ilandmark] = locator(1)$x
  }}

landmark.i.mean = colMeans(landmark.i)
ncol(landmark.i)
wbasisLM = create.bspline.basis(c(0,179), 3+ncol(landmark.i), 3,
                                c(0,landmark.i.mean,179))
WfdLM    = fd(matrix(0,3+ncol(landmark.i),1),wbasisLM)
WfdParLM = fdPar(WfdLM,1,0.00000000001)
regListLM = landmarkreg(mortalityaccfd, landmark.i,
                        landmark.i.mean, WfdParLM, TRUE)

mortality.warpfunc = regListLM$warpfd #stores the set of time warping functions in a functional data object

########################################################################################################################
########################################################################################################################
#this section creates figure 5.4 by plotting the four time warping functions for the aligned curves in figure 5.3, and saving into a jpeg to R's current working directory

jpeg("warping function for 4 counties.jpeg", units = "in", width = 6, height = 6 , res = 400)
par(cex.lab=1.8,cex.axis=1.5,mar=c(5,6,2,1)+.1)
plot(mortality.warpfunc[1:4],xlab="original time (day number)",ylab="transformed time (day number)",asp=1)
lines(mortality.warpfunc,lty=1,lwd=1.5)
title("first 4 time warping functions")
abline(0,1,lty=3)
dev.off()

########################################################################################################################
########################################################################################################################
#the next block of code creates figure 5.3 ,which is the four landmark registered, mortality derivative curves.

jpeg("landmark registration example.jpeg", units = "in", width = 12, height = 5, res = 400)
par(cex.axis=1.3,cex.lab=1.5,mfrow=c(1,2),mar = c(4,4.5,1,0.5))
plot(mortalityaccfd,xlab="day",ylab="Derivative of mortality",xlim=c(0,180))
lines(mortalityaccfd,lty=1,lwd=2)
title("unaligned")


colours = 1:78
registered.func = matrix(0,1771,4) #this matrix will store the values of x(h(t)) because we have to manually evaluate the registered function
#normally we should be able to plot it directly as an output from the landmarkreg() function but it doesn't seem to work

for(i in 1:4){ #for each of the 4 time warping functions we evaluate the original time values to produce the transformed time values
new.timepoints = eval.fd(seq(1,178,0.1),mortality.warpfunc[i]) #the full time interval of 0 to 180 days is not used here, as it seemed to produce numerical issues.
new.timepoints
registered.func[,i] = eval.fd(new.timepoints,mortalityaccfd[i]) #the original unaligned mortality curves arer then evaluated at the new transformed time points
}

#this then plots the values 
plot(seq(1,178,0.1),registered.func[,4],col="white",xlab="day",ylab="",xlim=c(0,180),ylim=c(-0.11,0.18))
title("aligned with landmark registration")
for(i in 1:78){
 lines(seq(1,178,0.1),registered.func[,i],lwd=1,col = colours[i])
}
abline(h=0,lty=2)
dev.off()

########################################################################################################################
########################################################################################################################

#this next section is slightly more complicated and uses the above code to produce figure 5.6 in the end.
#firstly the first block of code must be run again, but this time changed to identify 4 landmarks for all 78 curves

days2 = seq(0,179,length.out=1771) #
mortality.matrix2 = cbind(days2,registered.func) #registered.func here is the matrix of values of 78 landmark registered functions using 4 landmarks.
#this is created using the landmark registration code at the top of the page, changing it to use 78 curves instead the original 4 curves.
mortalitybasis2 = create.bspline.basis(c(0,total.days), nbasis, norder) #these next three lines resmooth the landmark registered functions into a functional data object, which can be inputted into register.fd for continuous registration
mortalityfdPar2 = fdPar(mortalitybasis2,2,0.0001) 
mortalityfd2 = smooth.basis(days2,mortality.matrix2[1:1771,2:79],mortalityfdPar2)$fd
#mortalityfd2 is the functional data object containing 78 landmark registered curves.

#now perform continuous registration on the already landmark registered functions 

land.then.cont = register.fd(mean.fd(mortalityfd2),mortalityfd2) 

#this block of code then creates figure 5.6, which is two plots side by side, one of the 78 original unregistered curves and one of the 78 landmark/continuously registered curves

jpeg("landmark then continuous reg of 78 curves.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.3,cex.lab=1.5,cex.pch=2,mfrow=c(1,2),mar = c(4.1,4.5,1,0.3))

plot(mortalityaccfd,xlim=c(0,180),asp=1,ylim=c(-0.11,0.18),xlab="Day Number",ylab="Rate of Change of Mortality") #change all to mortalityaccfd if we want to look at the derivative
lines(mortalityaccfd,lty=1,col = colours) #unaligned mortality derivative for 78 counties
title("Unaligned")

plot(land.then.cont$regfd,xlim=c(0,180),ylim=c(-0.11,0.18),xlab="Day Number",ylab="")
lines(land.then.cont$regfd,lty=1,lwd=1, col = colours)
title("Landmark then Continuous registration")
dev.off()

########################################################################################################################
########################################################################################################################

##### this function continuously aligns the 78 mortality curves and their derivative curves, and then I use it to create figure 5.7 which shows mean functions of the two.

continuous.reg2 <- function(start.region, end.region, total.days){

  number.regions = (end.region - start.region) + 1

  mortality.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions

  for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
    for (j in 1:(end.region-start.region+1)){ #there are 3 columns, each column j is a region
      if(number.regions > 1){
        mortality.matrix[i,j] <- mortality[i,j]
      }
      else {
        mortality.matrix[i,1] <- mortality[i,j]
      }
    }
  }

  mortalitybasis = create.bspline.basis(c(0,total.days), 19, 4)
  mortalityfdPar = fdPar(mortalitybasis,2,794)
  mortalityfd = smooth.basis((0:(total.days-1)),mortality.matrix,mortalityfdPar)$fd
  mortalityaccfd = deriv.fd(mortalityfd, 1)

  wbasisCR = create.bspline.basis(c(0,total.days-1), 15, 5) #warping function basis for continuous registration (CR)
  Wfd0CR   = fd(matrix(0,15,number.regions),wbasisCR)
  WfdParCR = fdPar(Wfd0CR, 2, 1)
  regList  = register.fd(mean.fd(mortality.posfd), #put in here what you want to continuously align
                         mortality.posfd, WfdParCR)
  mortalityaccfd.CR = regList$regfd #shifted mortality curves as a functional data object
  plot(mean.fd(mortalityaccfd.CR),xlab="Day number",ylab="Mortality Derivative (Deaths per 100,000)",ylim = c(-0.05,0.07))
  title(paste(""))
  #lines(mortalityaccfd.CR,lty=1,col = "black")
  lines(mean.fd(mortalityaccfd.CR),lty=1, col = "black",lwd = 3)
  lines(mean.fd(mortalityaccfd),lty=3,col="black",lwd = 3)
  abline(v = 37,lty = 5, col = "red")
  abline(v = 30,lty = 3, col = "orange")
  #mortality.plot(start.region,end.region,total.days,19,4)
  return(mortalityaccfd.CR)
}

#these two variables contain functional data objects for the 78 aligned mortality curves and derivative curves respectively.
mortality.contreg = continuous.reg2(1,78,180) #change register.fd in function back to non differentiated curves
mortalityderiv.contreg = continuous.reg2(1,78,180) #make sure you set the register.fd to align the derivatives

#this remaining section of code plots the mean curves for the two sets of aligned curves, and saves into a jpeg for figure 5.7
jpeg("aligned mean mortality curves.jpeg", units = "in", width = 10, height = 4.1, res = 400)
par(cex.axis=1.3,cex.lab=1.5,mfrow=c(1,2),mar = c(4.2,5,2.1,0.7))

plot(mean.fd(mortality.contreg),asp=1,xlab="Day Number",ylab="Mean mortality value")
lines(mean.fd(mortality.contreg),lty=1,lwd=3)
abline(v = 37,lty = 5, col = "red")
axis(1, at=37,labels=37,cex.axis = 0.7, tck = -0.06) #label on axis at x=37
abline(v = 30,lty = 3, col = "orange")
axis(1, at=30,labels=30, cex.axis = 0.7, tck = -0.06) #label on axis at x = 30

plot(mean.fd(mortalityderiv.contreg),asp=1,xlab="Day Number",ylab="Mean rate of change of mortality")
lines(mean.fd(mortalityderiv.contreg),lty=1,lwd=3)
abline(v = 37,lty = 5, col = "red")
axis(1, at=37,labels=37,cex.axis = 0.7, tck = -0.06) #label on axis at x=37
abline(v = 30,lty = 3, col = "orange")
axis(1, at=30,labels=30, cex.axis = 0.7, tck = -0.06) #label on axis at x = 30

dev.off()

##################