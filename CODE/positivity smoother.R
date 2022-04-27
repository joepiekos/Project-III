library(fda)
library(refund)

positivity <- read.csv("/Users/josephpiekos/Desktop/project III/data/english_positivity_final.csv")

counties <- names(positivity)
###############

#positivity.plot plots the positivity curves for a subset of the 78 counties we are considering, with up to 181 days of data able to plot.
#day 0 = 15/02/20
#day 180 = 13/08/20

positivity.plot <- function(start.region,end.region,total.days, dataset, nbasis = 180, norder = 4,lambda= 5){
  
  number.regions = (end.region - start.region) + 1 #number of replicates i.e. number of regions to plot on top of each other
  positivity.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions
  
  for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
    for (j in start.region:end.region){ #there are 3 columns, each column j is a region
      if(number.regions > 1){
        positivity.matrix[i,j-(start.region-1)] <- dataset[i,j]
      }
      else {
        positivity.matrix[i,1] <- dataset[i,j]
      }
    }  
  }
  positivity.matrix
  
  positivitybasis = create.bspline.basis(c(0,total.days-1), nbasis, norder)
  positivityfdPar = fdPar(positivitybasis,2,lambda) 
  positivityfd = smooth.basis((0:(total.days-1)),positivity.matrix,positivityfdPar)$fd
  plot(positivityfd,xlab="Day number",ylab="positivity (% of positive tests)",ylim=c(-1,100)) #max(eval.fd((0:(total.days-1)),positivityfd)))
  if (number.regions > 1){
    if(ncol(dataset)<11){
      title(paste("positivity curves for",number.regions,"regions")) #change back to this: positivity curves for",number.regions,"regions
    }
    else{
      title(paste("Positivity curves of metropolitan counties")) #change back to this: positivity curves for",number.regions,"counties
    }
    lines(positivityfd, lty=1, col="black",lwd=2)
    #lines(mean.fd(positivityfd),lty = 2,lwd=3,col="orange",ylim=c(-1,100))
  }
  else if(number.regions == 1 & ncol(dataset)<11){
    points(1:total.days,dataset[1:total.days,j],col="black",pch=20)
    lines(positivityfd, lty=1,lwd=3, col = "blue")
  }
  else{
    #lines(positivityfd, lty=1, col = "blue")
    points(1:total.days,dataset[1:total.days,j],col="black",pch=1,cex=0.8)
  }
  abline(v = 37,lty = 5, col = "red")
  axis(1, at=37,labels=37,cex.axis = 0.7, tck = -0.06) #label on axis at x=37
  abline(v = 30,lty = 3, col = "orange")
  axis(1, at=30,labels=30, cex.axis = 0.7, tck = -0.06) #label on axis at x = 30
  #legend(0, max(positivity.matrix)-1.6, legend = c("No non-essential contact and travel","First Lockdown"), col = c("orange","red"), lty = c(3,5), cex = 1, text.font = 2, title = "Key Announcements from Government")
}

positivity.plot(start.region = 72,end.region = 78, total.days = 181, dataset=positivity,nbasis = 90 ,norder = 4,lambda=5) #lambda = 0.2 is best according to gcv, but 5 is still low on the gcv scale and gives slightly nicer curves.
#positivity.plot

################################################################################################################################################################################################################################################################################################

#continuous.reg() is the same continuous registration function from the separate landmark & continuous registration R file.
#It will be used in creating figure 2.2, which contains some unaligned and aligned positivity curves
#NOTE: ensure when continuous.reg is run that the dataset parameter is set to "dataset = positivity"
#Continuous registration is a numerical method, and can take a long time if we want to align a lot of curves.
#The orange dotted line on the plot when this is run, is the mean curve of the aligned curves.

continuous.reg <- function(start.region, end.region, total.days, dataset = positivity){
  
  number.regions = (end.region - start.region) + 1  
  
  positivity.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions
  
  for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
    for (j in start.region:end.region){ #there are 3 columns, each column j is a region
      if(number.regions > 1){
        positivity.matrix[i,j-(start.region-1)] <- dataset[i,j]
      }
      else {
        positivity.matrix[i,1] <- dataset[i,j]
      }
    }
  }
  
  positivitybasis = create.bspline.basis(c(0,total.days-1), 19, 4)
  positivityfdPar = fdPar(positivitybasis,2,794) 
  positivityfd = smooth.basis((0:(total.days-1)),positivity.matrix,positivityfdPar)$fd
  
  wbasisCR = create.bspline.basis(c(0,total.days-1), 15, 5) #warping function basis for continuous registration (CR)
  Wfd0CR   = fd(matrix(0,15,number.regions),wbasisCR) 
  WfdParCR = fdPar(Wfd0CR, 2, 1)
  regList  = register.fd(mean.fd(positivityfd),
                         positivityfd, WfdParCR)
  positivityfd.CR = regList$regfd #shifted positivity curves as a functional data object
  plot(positivityfd.CR,xlab="Day number",ylab="positivity (% of positive tests)")
  title(paste("Aligned"))
  lines(positivityfd.CR,lty=1,col = "black")
  lines(mean.fd(positivityfd.CR),lty=2, col = "orange",lwd = 2.2)
  abline(v = 37,lty = 5, col = "red")
  abline(v = 30,lty = 3, col = "orange")
  #positivity.plot(start.region,end.region,total.days,19,4)
  return(positivityfd.CR)
}

#aligns all 78 positivity curves and stores it in a variable for easy plotting. So that we don't have to run the alignment again which takes a long time.
#note that this will take a while to align, up to 10 minutes depending on the PC used.
positivityfd.CR = continuous.reg(1,78,180)

#this code creates two figure 2.2, with the positivity curves of the metropolitan counties on the left and full set of 78 aligned curves on the right

jpeg("some positivity curves.jpeg", units = "in", width = 11, height = 4 , res = 400)
par(cex.axis=1.3,cex.lab=1.5,mfrow=c(1,2),mar = c(4,4.5,1,0.5))
positivity.plot(72,78,180,positivity,180,4,5)
plot(positivityfd.CR,xlab="Day number",ylab="")
title(paste("Full set of aligned positivity curves"))
lines(positivityfd.CR,lty=1,col = "black")
lines(mean.fd(positivityfd.CR),lty=1, col = "orange",lwd = 3.5)
abline(v = 37,lty = 5, col = "red")
abline(v = 30,lty = 3, col = "orange")
axis(1, at=37,labels=37,cex.axis = 0.7, tck = -0.06) #label on axis at x=37
axis(1, at=30,labels=30, cex.axis = 0.7, tck = -0.06) #label on axis at x = 30
dev.off()

################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################

#positivity.pos works in the same way as mortality.pos() from the mortality smoother.R code and plots a set of positively constrained positivity curves.


positivity.pos <- function(start.region,end.region,total.days, dataset, nbasis = 19, norder = 4,lambda=100){
  
  number.regions = (end.region - start.region) + 1 #number of replicates i.e. number of regions to plot on top of each other
  positivity.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions
  
  for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
    for (j in start.region:end.region){ #there are 3 columns, each column j is a region
      if(number.regions > 1){
        positivity.matrix[i,j] <- dataset[i,j]
      }
      else {
        positivity.matrix[i,1] <- dataset[i,j]
      }
    }
  }
  positivity.matrix
  
  positivitybasis = create.bspline.basis(c(0,total.days-1), nbasis, norder)
  positivityfdPar = fdPar(positivitybasis,2,lambda) 
  positivityfd = smooth.pos((0:(total.days-1)),as.matrix(positivity.matrix[,1]),positivityfdPar)
  Wfd = positivityfd$Wfdobj
  positivity.matrix
  precfit = exp(eval.fd((0:(total.days-1)), Wfd))
  plot(1:total.days,positivity[1:total.days,1],col="white",pch=20,ylim=c(0,100),xlab="Day number",ylab="positivity (Deaths per 100,000)")
  lines((0:(total.days-1)), precfit) 
  for(k in 2:number.regions){ #this for loop creates the rest of the curves other than the first
    positivityfd = smooth.pos((0:(total.days-1)),as.matrix(positivity.matrix[,k]),positivityfdPar)
    Wfd = positivityfd$Wfdobj
    positivity.matrix
    precfit = exp(eval.fd((0:(total.days-1)), Wfd))
    #plot(1:total.days,positivity[1:total.days,1],col="white",pch=20,ylim=c(-0.01,max(precfit)),xlab="Day number",ylab="positivity (Deaths per 100,000)")
    lines((0:(total.days-1)), precfit) 
  }
  abline(v = 37,lty = 5, col = "red")
  axis(1, at=37,labels=37,cex.axis = 0.7, tck = -0.03) #label on axis at x=37
  abline(v = 30,lty = 3, col = "orange")
  axis(1, at=30,labels=30, cex.axis = 0.7, tck = -0.03) #label on axis at x = 30
  abline(h = 0,lty = 5, col = "black")
  legend(80, 80, legend = c("No non-essential contact and travel","First Lockdown"), col = c("orange","red"), lty = c(3,5), cex = 0.9, text.font = 2, title = "Key Announcements from Government")
}

dev.off() #clears previous graphical parameters
positivity.pos(1,2,180,positivity,180,4,5)

