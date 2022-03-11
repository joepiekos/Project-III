library(fda)
library(refund)

regions <- read.csv("/Users/josephpiekos/Desktop/project III/data/english authorities.csv")
mobility <- read.csv("/Users/josephpiekos/Desktop/project III/data/english_mobility_final.csv")

region = vector("list",78) #list of regions we are looking at
for (i in 1:78){
  region[[ i]] = regions[i,]
}
region

mobility.plot <- function(start.region,end.region,total.days, mobility.category, K, b.order = 4,lambda){

number.regions = (end.region - start.region) + 1 #number of replicates i.e. number of regions to plot on top of each other
mobility.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions

for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
  for (j in start.region:end.region){ #there are 3 columns, each column j is a region
    if(number.regions > 1){
      mobility.matrix[i,j] <- mobility.category[181*(j-1)+i]
      }
    else {
      mobility.matrix[i,1] <- mobility.category[181*(j-1)+i]
      }
    }
  }
mobility.matrix

mobilitybasis = create.bspline.basis(c(0,total.days-1), K,b.order)
mobilityfdPar = fdPar(mobilitybasis,2,lambda)
mobilityfd = smooth.basis((0:(total.days-1)),mobility.matrix,mobilityfdPar)$fd
mobility.gcv = smooth.basis((0:(total.days-1)),mobility.matrix,mobilitybasis)$gcv

mobility.gcv
mobility.matrix
plot(mobilityfd,xlab="Day number",ylab="Mobility")
if (number.regions > 1){
  title("Mobility Curves")
  lines(mobilityfd, lty=1, col = "black")
  lines(mean.fd(mobilityfd), lty=2, col = "orange", lwd = 3)
}
else {
  lines(mobilityfd, lty=1, col = " dark green", lwd = 2)
}
abline(v = 37,lty = 5, col = "red")
axis(1, at=37,cex.axis = 0.7, tck = -0.03) #label on axis at x=37
abline(v = 30,lty = 3, col = "orange")
axis(1, at=30, cex.axis = 0.7, tck = -0.03) #label on axis at x = 30
legend(60, max(mobility.matrix)-5, legend = c("No non-essential contact and travel","First Lockdown"), col = c("orange","red"), lty = c(3,5), cex = 1.6, text.font = 2, title = "Key Announcements from Government")
}

jpeg("50 mobility curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)
mobility.plot(start.region = 1,end.region = 50, total.days = 90, mobility.category = mobility$grocery_and_pharmacy_percent_change_from_baseline,  K = 90, b.order = 4,lambda = 10)
dev.off()

single.region.mobility <- function(total.days, K, b.order = 4,lambda,region.name){
  start.region = match(region.name, region)
  end.region = start.region
  mobility.plot(start.region,end.region,total.days, mobility.category = mobility$grocery_and_pharmacy_percent_change_from_baseline, K, b.order,lambda)
  #title(main = paste("Mobility Curve for", region[[start.region]], "(Groceries and Pharmacy)"))#CHANGE DEPENDING ON ACTIVITY
}

jpeg("mobility curve for county durham.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=2,cex.axis=1.5,mar=c(4.4,5.3,1,0.5)+.1)
single.region.mobility(180,90,4,1,"County Durham")
points(x=0:180,y=mobility[2535:2715,6])
dev.off()

#######
continuous.reg <- function(start.region, end.region, total.days,mobility.category){
  
  number.regions = (end.region - start.region) + 1  
  
  mobility.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions
  
  for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
    for (j in start.region:end.region){ #there are 3 columns, each column j is a region
      if(number.regions > 1){
        mobility.matrix[i,j] <- mobility.category[181*(j-1)+i]
      }
      else {
        mobility.matrix[i,1] <- mobility.category[181*(j-1)+i]
      }
    }
  }
  
  mobilitybasis = create.bspline.basis(c(0,total.days-1), 19, 4)
  mobilityfdPar = fdPar(mobilitybasis,2,10) 
  mobilityfd = smooth.basis((0:(total.days-1)),mobility.matrix,mobilityfdPar)$fd
  
  wbasisCR = create.bspline.basis(c(0,total.days-1), 15, 5) #warping function basis for continuous registration (CR)
  Wfd0CR   = fd(matrix(0,15,number.regions),wbasisCR) 
  WfdParCR = fdPar(Wfd0CR, 2, 1)
  regList  = register.fd(mean.fd(mobilityfd),
                         mobilityfd, WfdParCR)
  mobilityfd.CR = regList$regfd #shifted mobility curves as a functional data object
  plot(mobilityfd.CR,xlab="Day number",ylab="Percentage change from Baseline")
  
  lines(mobilityfd.CR,lty=1, col ="black")
  lines(mean.fd(mobilityfd.CR), lty=2, col = "orange", lwd = 3)
  abline(v = 37,lty = 5, col = "red")
  abline(v = 30,lty = 3, col = "orange")
  #mobility.plot(start.region,end.region,total.days,19,4)
}

jpeg("50 aligned mobility curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)
continuous.reg(1,50,90,mobility$grocery_and_pharmacy_percent_change_from_baseline)
dev.off()



