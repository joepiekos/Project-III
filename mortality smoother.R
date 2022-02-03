library(fda)
library(refund)
regions <- read.csv("/Users/josephpiekos/Desktop/project III/english authorities.csv")
mortality <- read.csv("/Users/josephpiekos/Desktop/project III/english_mortality_final.csv")


###############

#mortality.plot plots the mortality curves for a subset of the 78 counties we are considering, with up to 181 days of data able to plot.
#day 0 = 15/02/20
#day 180 = 13/08/20

mortality.plot <- function(start.region,end.region,total.days, nbasis = 30, norder = 4){
  
  number.regions = (end.region - start.region) + 1 #number of replicates i.e. number of regions to plot on top of each other
  
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
  mortality.matrix
  
  mortalitybasis = create.bspline.basis(c(0,total.days-1), nbasis, norder)
  mortalityfdPar = fdPar(mortalitybasis,2,10) 
  mortalityfd = smooth.basis((0:(total.days-1)),mortality.matrix,mortalityfdPar)$fd
  mortality.matrix
  plot(mortalityfd,xlab="Day number",ylab="Mortality (Deaths per 100,000)")
  if (number.regions > 1){
    #title("mortality Curves")
    lines(mortalityfd, lty=1)
    lines(mean.fd(mortalityfd),lty = 2)
  }
  else {
    lines(mortalityfd, lty=1, col = "blue")
  }
  abline(v = 37,lty = 5, col = "red")
  #axis(1, at=37,labels=37,cex.axis = 0.7, tck = -0.06) label on axis at x=37
  abline(v = 30,lty = 3, col = "orange")
  #axis(1, at=30,labels=30, cex.axis = 0.7, tck = -0.06) label on axis at x = 30
  #legend(total.days-45, max(mortality.matrix)-5, legend = c("No non-essential contact and travel","First Lockdown"), col = c("orange","red"), lty = c(3,5), cex = 1.3, text.font = 2, title = "Key Announcements from Government")
}

mortality.plot(start.region = 1,end.region = 50, total.days = 90, nbasis = 19 ,norder = 4)

##################

#single.region.mortality plots the mortality curve for a single region, the desired region is named in the function parameters

single.region.mortality <- function(total.days, nbasis = 65 , norder = 4,region.name){
  start.region = match(region.name, region)
  end.region = start.region
  mortality.plot(start.region,end.region,total.days, nbasis, norder)
  title(main = paste("Mortality Curve for", region[[start.region]]))#CHANGE DEPENDING ON ACTIVITY
}
single.region.mortality(40,19,4,"Bath and North East Somerset")
mortality[1,1]

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


#jpeg("50 aligned mortality curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
#par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)
continuous.reg(1,50,90)
#lines(mean.fd(accelfdCR),lty=2, col = "orange",lwd = 4)
#dev.off()


#jpeg("50 mortality curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
#par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)
mortality.plot(1,50,90,19,4)
#dev.off()
#plot()
