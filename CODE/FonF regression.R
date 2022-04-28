library(fda)
library(refund)

######### Code in this section is used to create figure 7.1 #######################

mortality.matrix2 <- t(as.matrix(mortality)) #we must use the transpose of the data matrices when using the refund package for fitting functional models
positivity.matrix2 <- t(as.matrix(positivity[,1:78]))
positivity.matrix2
mobility <- read.csv("/Users/josephpiekos/Desktop/project III/data/mobilitymatrix.csv")

#this line fits the total functional model for mortality ~ positivity
#it takes about a minute to run sometimes so don't worry if nothing happens straight away
m3_ff <- pffr(mortality.matrix2 ~ ff(positivity.matrix2,xind=seq(1,181,l=181))) #creates regression fit

#model summary if interested
summary(m3_ff)

#this line stores all the information about the coefficient surface for the positivity regression model. psi_plot3$fit will contain the z-values of the surface
psi_plot3 <- plot(m3_ff, select = 3, pers=TRUE)[[2]]

#to fit the other total functional model of mortality ~ mobility, I directly smooth the coefficient surface using a double basis. As I describe in section 7.3 of my report.

# the code up to line 75 smooths the three variables again into a functional data objects of 78 curves, for model fitting.
number.regions = 78
total.days = 180
mortality.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions

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
mortalityfdPar = fdPar(mortalitybasis,2,800) #lambda = 800 selected using GCV
mortalityfd = smooth.basis((0:(total.days-1)),mortality.matrix,mortalityfdPar)$fd


mobility.matrix <- (matrix(0,total.days,number.regions))

for(i in 1:total.days){ #there are 60 rows, as you go down i represents another day
  for (j in 1:78){ #there are 3 columns, each column j is a region
    if(number.regions > 1){
      mobility.matrix[i,j] <- mobility[i,j]
    }
    else {
      mobility.matrix[i,1] <- mobility[i,j]
    }
  }
}

mobilitybasis = create.bspline.basis(c(0,total.days), nbasis=90, norder = 4)
mobilityfdPar = fdPar(mobilitybasis,2,1) #lambda = 1 chosen already using GCV
mobilityfd = smooth.basis((0:(total.days-1)),mobility.matrix,mobilityfdPar)$fd

positivity.matrix <- (matrix(0,total.days,number.regions))

for(i in 1:total.days){ 
  for (j in 1:78){
    if(number.regions > 1){
      positivity.matrix[i,j] <- positivity[i,j]
    }
    else {
      positivity.matrix[i,1] <- positivity[i,j]
    }
  }
}

positivitybasis = create.bspline.basis(c(0,total.days), nbasis=90, norder = 4)
positivityfdPar = fdPar(positivitybasis,2,5) #lambda = 1 chosen already using GCV
positivityfd = smooth.basis((0:(total.days-1)),positivity.matrix,positivityfdPar)$fd


# I now begin the set up for fitting the mobility model using the double basis expansion

beta0Par  = fdPar(create.bspline.basis(c(0,180),10,4) , 2, 0.1) #set up the basis for the intercept function in the model

sbasis.n = 16 # the number of basis functions for the phi(s) basis of the coefficient surface. truncate this to avoid overfitting
tbasis.n = 15 # the number of basis functions for the theta(t) basis of the coefficient surface. trunacte this to ensure predictions are smooth

#sets up a bivariate functional data object containing the two bases we will use to expand the coefficient surface
bifd.st <- bifd(matrix(0,sbasis.n,tbasis.n),create.bspline.basis(c(0,180),sbasis.n),create.bspline.basis(c(0,180),tbasis.n))

bifdpar.st <- bifdPar(bifd.st)

betaList  = list(beta0Par, bifdpar.st) #puts all the bases into one list for use in the model fitting function

linmodSmooth = linmod(mobilityfd, mortalityfd,betaList) #this function fits the total linear model for mortality ~ mobility.
#both functional objects here can have registration applied to them first before putting into the model

bifd_mat  = eval.bifd(seq(0,180,length.out = 60),seq(0,180,length.out = 60), linmodSmooth$beta1estbifd) #the linmod smooth function returns a bivariate functional data object for the coefficient surface
#eval.bifd() evaluates the coefficient surface at a grid of points for eventual plotting.

#this block of code up to line 88 then plots the two coefficient surfaces seen in figure 7.1
jpeg("coefficient surfaces for positivity and mobility.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=0.7,cex.lab=1,mfrow=c(1,2),mar = c(1.5,2.5,0,0.5))

#this persp() function plots the coefficient surface for the positivity model, using the information from psi_plot3
persp(psi_plot3$x, psi_plot3$y, matrix(psi_plot3$fit, 40, 40),
      theta=-45, phi=25,r=2.9,expand = 0.5,
      ticktype="detailed",col="orange", shade = 0.2,border="black",xlab="Day (Positivity)",ylab="Day (Mortality)",zlab="")

#this persp() function plots the coefficient surface for the mobility model, using the information from bifd_mat
persp(seq(0,180,length.out =60),seq(0,180,length.out = 60), bifd_mat,
      theta=-45, phi=25, r=2.9,expand = 0.5,
      ticktype="detailed",col="orange", shade = 0.2,border="black",xlab="Day (Mobility)",ylab="Day (Mortality)",zlab="")

dev.off()



########################################################################################################################################################################################
### This section of code (line 98 to line 131) creates figure 7.3: Plots of the R^2 functions of both the mobility and positivity model ####

jpeg("r^2 functions for mobility then positivity.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.5,cex.lab=1.5,mfrow=c(1,2),mar = c(4.2,4.5,2,0.5))
y = mortality.contreg #continuously aligned mortality functions from the curve registration code
ybar = mean.fd(y)
year = seq(0,180,by=0.1)
mortality.finemat = eval.fd(year,y) #evaluates the aligned mortality curves at a fine mesh of points up to 180
mortality.meanvec = eval.fd(year,ybar) #similarly evaluates the mean function
mortality.hatfinemat = eval.fd(year,linmodSmooth$yhatfdobj)  #similarly evaluates the fitted functions from the mobility model
resmat = mortality.finemat - mortality.hatfinemat #the next few lines of code represent the terms in the function for R^2 as seen in expression (7.8) on page 49 of my report. resmat is named as a shortened version of "residual matrix"
resmat0 = mortality.finemat - mortality.meanvec %*% matrix(1,1,78)
SSE0 = apply((resmat0)^2,1,sum)
SSE1 = apply(resmat^2,1,sum)
Rsqr = (SSE0-SSE1)/SSE0 #calculates the R
plot(Rsqr,ylim=c(-0.2,1),col="white",xaxt="n",ylab=expression(R^2),xlab="Day",main=expression(paste(R^2," function (mobility as predictor)")))
axis(1,at=seq(0,180/0.1,by=50/0.1),labels=c("0","50","100","150"))
abline(h=0,lty=3)
lines(Rsqr,lwd=2)

linmodSmooth2 = linmod(positivityfd, mortalityfd,betaList) 
#fit the positivity ~ mortality model again, using linmod, so that it returns a functional data object of the fitted functions. 
#The pffr function before, did the same thing, but plot a better surface but doesn't return the fitted functions as nicely

#repeat the process above to plot the second R^2 function
y = mortality.contreg
ybar = mean.fd(y)
year = seq(0,180,by=0.1)
mortality.finemat = eval.fd(year,y)
mortality.meanvec = eval.fd(year,ybar)
mortality.hatfinemat = eval.fd(year,linmodSmooth2$yhatfdobj)
resmat = mortality.finemat - mortality.hatfinemat
resmat0 = mortality.finemat - mortality.meanvec %*% matrix(1,1,78)
SSE0 = apply((resmat0)^2,1,sum)
SSE1 = apply(resmat^2,1,sum)
Rsqr = (SSE0-SSE1)/SSE0
plot(Rsqr,ylim=c(-0.2,1),col="white",xaxt="n",ylab=expression(R^2),xlab="Day",main=expression(paste(R^2," function (positivity as predictor)")))
axis(1,at=seq(0,180/0.1,by=50/0.1),labels=c("0","50","100","150"))
abline(h=0,lty=3)
lines(Rsqr,lwd=2)
dev.off()

#################################################################################################################################################################
#################################################################################################################################################################

#### The code up to line 173 plots two of the predicted mortality curves from the mobility model, each on top of the actual observed mortality curve from the data#

jpeg("Predicted mortality curves.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.2,cex.lab=1.5,mfrow=c(1,2),mar = c(4.5,2.5,2,0.5))
names = c("County Durham", "Kent")
for(i in c(15,30)){
mortality.plot(i,i,180,mortality,lambda=700)
title(paste(names[i/15],"mortality curve prediction"))
lines(linmodSmooth$yhatfdobj[i],lty=2) #ylim=c(-0.01,2.2)
legend(80, 2.1, legend = c("Predicted Curve","Actual Curve"), col = c("black","blue"), lty = c(2,1), cex = 1, text.font = 2)
}
dev.off()

#################################################################################################################################################################
#################################################################################################################################################################





