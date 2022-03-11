library(fda)
library(refund)
regions <- read.csv("/Users/josephpiekos/Desktop/project III/data/english authorities.csv")
mortality <- read.csv("/Users/josephpiekos/Desktop/project III/data/english_mortality_final.csv")
mortality.regions <- read.csv("/Users/josephpiekos/Desktop/project III/data/regions_mortality_final.csv")
CanadianWeather
###############
dev.off()
#mortality.plot plots the mortality curves for a subset of the 78 counties we are considering, with up to 181 days of data able to plot.
#day 0 = 15/02/20
#day 180 = 13/08/20

mortality.plot <- function(start.region,end.region,total.days, dataset, nbasis = 19, norder = 4,lambda=100){

  number.regions = (end.region - start.region) + 1 #number of replicates i.e. number of regions to plot on top of each other
  mortality.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions
  
  for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
    for (j in start.region:end.region){ #there are 3 columns, each column j is a region
      if(number.regions > 1){
        mortality.matrix[i,j] <- dataset[i,j]
      }
      else {
        mortality.matrix[i,1] <- dataset[i,j]
      }
    }
  }
  mortality.matrix
  
  mortalitybasis = create.bspline.basis(c(0,total.days-1), nbasis, norder)
  mortalityfdPar = fdPar(mortalitybasis,2,lambda) 
  mortalityfd = smooth.basis((0:(total.days-1)),mortality.matrix,mortalityfdPar)$fd
  plot(mortalityfd,xlab="Day number",ylab="Mortality (Deaths per 100,000)",ylim=c(-0.01,2.2)) #max(eval.fd((0:(total.days-1)),mortalityfd)))
  if (number.regions > 1){
    if(ncol(dataset)<11){
      title(paste("Mortality curves for",number.regions,"regions"))
    }
    else{
      title(paste("Mortality curves for",number.regions,"counties"))
    }
    lines(mortalityfd, lty=1, col="black")
    #lines(mean.fd(mortalityfd),lty = 2,lwd=3,col="orange")
  }
  else if(number.regions == 1 & ncol(dataset)<11){
    points(1:total.days,dataset[1:total.days,j],col="black",pch=20)
    lines(mortalityfd, lty=1,lwd=3, col = "blue")
  }
  else{
    lines(mortalityfd, lty=1, col = "blue")
  }
  abline(v = 37,lty = 5, col = "red")
  axis(1, at=37,labels=37,cex.axis = 0.7, tck = -0.06) #label on axis at x=37
  abline(v = 30,lty = 3, col = "orange")
  axis(1, at=30,labels=30, cex.axis = 0.7, tck = -0.06) #label on axis at x = 30
  #legend(0, max(mortality.matrix)-1.6, legend = c("No non-essential contact and travel","First Lockdown"), col = c("orange","red"), lty = c(3,5), cex = 1, text.font = 2, title = "Key Announcements from Government")
}

mortality.plot(start.region = 1,end.region = 9, total.days = 180, dataset=mortality,nbasis = 80 ,norder = 4,lambda=100)
#mortality.plot

################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################

#single.region.mortality plots the mortality curve for a single region, the desired region is named in the function parameters
match("GREATER LONDON",region)
single.region.mortality <- function(total.days, nbasis = 90 , norder = 4,region.name,lambda){
  start.region = match(region.name, region)
  end.region = start.region
  mortality.plot(start.region,end.region,total.days,dataset=mortality, nbasis, norder,lambda)
  title(main = paste("Mortality Curve for", region[[start.region]]))#CHANGE DEPENDING ON ACTIVITY
}

single.region.mortality(90,90,4,"Kent",100)

################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################

#continuous.reg plots mortality curves, and also plots shifted mortality curves so we can compare.

continuous.reg <- function(start.region, end.region, total.days){

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

mortalitybasis = create.bspline.basis(c(0,total.days-1), 19, 4)
mortalityfdPar = fdPar(mortalitybasis,2,10) 
mortalityfd = smooth.basis((0:(total.days-1)),mortality.matrix,mortalityfdPar)$fd

wbasisCR = create.bspline.basis(c(0,total.days-1), 15, 5) #warping function basis for continuous registration (CR)
Wfd0CR   = fd(matrix(0,15,number.regions),wbasisCR) 
WfdParCR = fdPar(Wfd0CR, 2, 1)
regList  = register.fd(mean.fd(mortalityfd),
                       mortalityfd, WfdParCR)
mortalityfd.CR = regList$regfd #shifted mortality curves as a functional data object
plot(mortalityfd.CR,xlab="Day number",ylab="Mortality (Deaths per 100,000)")

lines(mortalityfd.CR,lty=1,col = "black")
lines(mean.fd(mortalityfd.CR),lty=2, col = "orange",lwd = 7)
abline(v = 37,lty = 5, col = "red")
abline(v = 30,lty = 3, col = "orange")
#mortality.plot(start.region,end.region,total.days,19,4)
}

################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################

mortality.pos <- function(start.region,end.region,total.days, dataset, nbasis = 19, norder = 4,lambda=100){
  
  number.regions = (end.region - start.region) + 1 #number of replicates i.e. number of regions to plot on top of each other
  mortality.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions
  
  for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
    for (j in start.region:end.region){ #there are 3 columns, each column j is a region
      if(number.regions > 1){
        mortality.matrix[i,j] <- dataset[i,j]
      }
      else {
        mortality.matrix[i,1] <- dataset[i,j]
      }
    }
  }
  mortality.matrix
  
  mortalitybasis = create.bspline.basis(c(0,total.days-1), nbasis, norder)
  mortalityfdPar = fdPar(mortalitybasis,2,lambda) 
  mortalityfd = smooth.pos((0:(total.days-1)),as.matrix(mortality.matrix[,1]),mortalityfdPar)
  Wfd = mortalityfd$Wfdobj
  mortality.matrix
  precfit = exp(eval.fd((0:(total.days-1)), Wfd))
  plot(1:total.days,mortality[1:total.days,1],col="white",pch=20,ylim=c(-0.01,2.2),xlab="Day number",ylab="Mortality (Deaths per 100,000)")
  lines((0:(total.days-1)), precfit) 
  for(k in 2:number.regions){ #this for loop creates the rest of the curves other than the first
  mortalityfd = smooth.pos((0:(total.days-1)),as.matrix(mortality.matrix[,k]),mortalityfdPar)
  Wfd = mortalityfd$Wfdobj
  mortality.matrix
  precfit = exp(eval.fd((0:(total.days-1)), Wfd))
  #plot(1:total.days,mortality[1:total.days,1],col="white",pch=20,ylim=c(-0.01,max(precfit)),xlab="Day number",ylab="Mortality (Deaths per 100,000)")
  lines((0:(total.days-1)), precfit) 
  }
  
  # for(k in 2:number.regions){
  #   
  # }
  # mortalityfd = smooth.basis((0:(total.days-1)),mortality.matrix,mortalityfdPar)$fd
  # plot(mortalityfd,xlab="Day number",ylab="Mortality (Deaths per 100,000)",ylim=c(-0.05,max(eval.fd((0:(total.days-1)),mortalityfd))))
  # if (number.regions > 1){
  #   if(ncol(dataset)<11){
  #     title(paste("Mortality curves for",number.regions,"regions"))
  #   }
  #   else{
  #     title(paste("Mortality curves for",number.regions,"counties"))
  #   }
  #   lines(mortalityfd, lty=1, col="black")
  #   lines(mean.fd(mortalityfd),lty = 2,lwd=3,col="orange")
  # }
  # else if(number.regions == 1 & ncol(dataset)<11){
  #   points(1:total.days,dataset[1:total.days,j],col="black",pch=20)
  #   lines(mortalityfd, lty=1,lwd=3, col = "blue")
  # }
  # else{
  #   lines(mortalityfd, lty=1, col = "blue")
  # }
  abline(v = 37,lty = 5, col = "red")
  axis(1, at=37,labels=37,cex.axis = 0.7, tck = -0.03) #label on axis at x=37
  abline(v = 30,lty = 3, col = "orange")
  axis(1, at=30,labels=30, cex.axis = 0.7, tck = -0.03) #label on axis at x = 30
  abline(h = 0,lty = 5, col = "black")
  legend(80, 2.1, legend = c("No non-essential contact and travel","First Lockdown"), col = c("orange","red"), lty = c(3,5), cex = 1.3, text.font = 2, title = "Key Announcements from Government")
}
plot.new()
jpeg("10 nonnegative mortality curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)
mortality.pos(1,10,180,mortality,90,4,5) #if using mortality.regions, make sure not to choose more than 9 regions
dev.off()

mortality.plot(1,9,180,mortality,90,4,5)
 
?axis


#jpeg("50 aligned mortality curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
#par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)
continuous.reg(1,50,90)
#lines(mean.fd(accelfdCR),lty=2, col = "orange",lwd = 4)
#dev.off()


jpeg("50 mortality curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)
mortality.plot(1,50,90,19,4)
dev.off()
#plot()


mortalityfd = smooth.pos((0:(total.days-1)),mortality.matrix,mortalityfdPar)
Wfd = mortalityfd$Wfdobj
mortality.matrix
precfit = exp(eval.fd((0:(total.days-1)), Wfd))
plot(1:total.days,mortality[1:total.days,1],col="black",pch=20)
lines((0:(total.days-1)), precfit,lwd=2) 
