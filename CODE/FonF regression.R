library(fda)
library(refund)
#code in this file is used to create figure 7.1

mortality.matrix2 <- t(as.matrix(mortality)) #we must use the transpose of the data matrices when using the refund package for fitting functional models
positivity.matrix2 <- t(as.matrix(positivity[,1:78]))
positivity.matrix2
mobility2 <- read.csv("/Users/josephpiekos/Desktop/project III/data/mobility.csv")
mobility.matrix2 <- t(as.matrix(mobility2))

#this line fits the total functional model for mortality ~ positivity
#it takes about a minute to run sometimes so don't worry if nothing happens straight away
m3_ff <- pffr(mortality.matrix2 ~ ff(positivity.matrix2,xind=seq(1,181,l=181))) #creates regression fit

#model summary if interested
summary(m3_ff)

#this line stores all the information about the coefficient surface for the positivity regression model. psi_plot3$fit will contain the z-values of the surface
psi_plot3 <- plot(m3_ff, select = 3, pers=TRUE)[[2]]

#to fit the other total functional model of mortality ~ mobility, I directly smooth the coefficient surface using a double basis. As I describe in section 7.3 of my report.

#firstly the code up to line  smooths the mortality data again into a functional data object of 78 curves
number.regions = 78
total.days = 180
mortality.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions
mobility.matrix <- (matrix(0,total.days,number.regions))

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
mobility.matrix

mobilitybasis = create.bspline.basis(c(0,total.days), nbasis=90, norder = 4)
mobilityfdPar = fdPar(mobilitybasis,2,1) #lambda = 1 chosen already using GCV
mobilityfd = smooth.basis((0:(total.days-1)),mobility.matrix,mobilityfdPar)$fd

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




### R^2 function ####
jpeg("r^2 functions for mobility then positivity.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.5,cex.lab=1.5,mfrow=c(1,2),mar = c(4.2,4.5,2,0.5))
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
plot(Rsqr,ylim=c(-0.2,1),col="white",xaxt="n",ylab=expression(R^2),xlab="Day",main=expression(paste(R^2," function (mobility as predictor)")))
axis(1,at=seq(0,180/0.1,by=50/0.1),labels=c("0","50","100","150"))
abline(h=0,lty=3)
lines(Rsqr,lwd=2)

y = mortality.contreg
ybar = mean.fd(y)
year = seq(0,180,by=0.1)
mortality.finemat = eval.fd(year,y)
mortality.meanvec = eval.fd(year,ybar)
mortality.hatfinemat = eval.fd(year,linmodSmooth3$yhatfdobj)
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


#######################
plot(y)
lines(ybar)


persp(seq(0,180,length.out =60),seq(0,180,length.out = 60), bifd_mat,
      theta=-25, phi=25, r=2,expand = 0.4,
      ticktype="detailed",col="orange", shade = 0.2,border="black",xlab="Day (Mobility)",ylab="Day (Mortality)",zlab="")
contour(seq(0,180,length.out = 60),seq(0,180,length.out = 60), bifd_mat)


#####concurrent model######
betabasis = create.bspline.basis(c(0,180),180,4) 
xfdlist = list(rep(1,78),mobilityfd)
betafdPar = fdPar(betabasis,2,lambda=100)
betalist = list(betafdPar,betafdPar)
fRegress.mod = fRegress(mortalityfd,xfdlist,betalist)
##########################

####checking prediction of the models#
plot(fRegress.mod$yhatfdobj)
plot(linmodSmooth$yhatfdobj)
mortality.plot(1,78,180,mortality,lambda=800)

mobility.plot(start.region = 1,end.region = 78, total.days = 180, K = 100, b.order = 4,lambda = 100)

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
####################################

jpeg("mortality boxplot.jpeg", units = "in", width = 10, height = 6.67 , res = 400)
par(cex.lab=1.8,cex.axis=1.7,mar=c(5,5,1,1)+.1,cex.main = 2.5)
mortalityfd.boxplot = boxplot.fd(mortalityfd,xlab="Day",ylab="Mortality (Deaths per 100,000 people)") 
mortalityfd.boxplot
legend(100, 3.1, legend = c("Outliers"), col = c("red"), lty = c(2),lwd = c(1.8), cex = 1.9, text.font = 2)
dev.off()

mortalityfd.boxplot$medcurve
mortality.matrix2

mortality.pca = pca.fd(mortality.contreg,2)
mortality.rotpca = varmx.pca.fd(mortality.pca)
plot.pca.fd(mortality.rotpca)

?cor.fd
mortality.covfd = cor.fd(mortality.posfd)
day.seq = seq(30,100,length.out = 50)
mortality.covfunc = eval.bifd(day.seq,day.seq,mortality.covfd)
min(mortality.covfunc)
jpeg("mortality covariance.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.2,cex.lab=1.5,mfrow=c(1,2),mar = c(4.5,4.5,2,0))
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



