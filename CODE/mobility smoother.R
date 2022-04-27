library(fda)
library(refund)

mobility <- read.csv("/Users/josephpiekos/Desktop/project III/data/mobilitymatrix.csv") #loads mobility data as a 181 x 78 data frame, where each column represents a different county.
counties <- names(mobility)

#mobility.plot() plots a desired number of mobility curves, smoothed using a roughness penalty.
#nbasis here is the number of basis functions, 
mobility.plot <- function(start.region,end.region,total.days, nbasis, norder = 4,lambda){

number.regions = (end.region - start.region) + 1 #number of replicates i.e. number of regions to plot on top of each other
mobility.matrix <- (matrix(0,total.days,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions

for (i in 1:total.days){ #there are 60 rows, as you go down i represents another day
  for (j in start.region:end.region){ #there are 3 columns, each column j is a region
    if(number.regions > 1){
      mobility.matrix[i,j-(start.region-1)] <- mobility[i,j]
      }
    else {
      mobility.matrix[i,1] <- mobility[i,j]
      }
    }
  }
mobility.matrix

mobilitybasis = create.bspline.basis(c(0,total.days-1), nbasis,norder)
mobilityfdPar = fdPar(mobilitybasis,2,lambda)
mobilityfd = smooth.basis((0:(total.days-1)),mobility.matrix,mobilityfdPar)$fd
plot(mobilityfd,xlab="Day",ylab="Mobility",ylim=c(-50,25))
if (number.regions > 1){
  title(paste(number.regions,"Mobility Curves"))
  lines(mobilityfd, lty=1, col = "black")
  #lines(mean.fd(mobilityfd), lty=2, col = "red", lwd = 3,ylim=c(-45,25))
}
else { #if only one curve is being plotted then this code runs, where it adds the data points as well as the curve
  lines(mobilityfd, lty=1, col = "black", lwd = 2)
  points(0:180,mobility[,start.region],cex=0.8)
}
abline(v = 37,lty = 5, col = "red")
axis(1, at=37,cex.axis = 0.9, tck = -0.03) #label on axis at x=37
abline(v = 30,lty = 3, col = "orange")
axis(1, at=30, cex.axis = 0.9, tck = -0.03) #label on axis at x = 30
legend(60, max(mobility.matrix)-5, legend = c("No non-essential contact and travel","First Lockdown"), col = c("orange","red"), lty = c(3,5), cex = 0.9, text.font = 2, title = "Key Announcements from Government")
}

#this code creates figure 2.1, and saves into a jpeg in R's working directory
jpeg("county durham mobility curve.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)
mobility.plot(start.region = 14,end.region = 15, total.days = 180, nbasis = 180, norder = 4,lambda = 1)
dev.off()


##############################################################################################################################
#least squares vs roughness penalty (for creating figure 4.1)
##############################################################################################################################

#smooth.squares() uses least squares smoothing instead of roughness penalty smoothing 
#and plots the mobility curve of durham using this method

smooth.squares <- function(){
number.regions = 78 #number of replicates i.e. number of regions to plot on top of each other
mobility.matrix <- (matrix(0,180,number.regions)) #creates empty matrix, number of rows = time period, number of columns = number of regions

for (i in 1:180){ #there are 60 rows, as you go down i represents another day
  for (j in 1:78){ #there are 3 columns, each column j is a region
    if(number.regions > 1){
      mobility.matrix[i,j] <- mobility[i,j]
    }
    else {
      mobility.matrix[i,1] <- mobility[i,j]
  }
}}

mobilitybasis = create.bspline.basis(c(0,180-1),nbasis=90,norder=4)
mobilitybasismat = eval.basis((0:(180-1)), mobilitybasis)
mobilitycoef = solve(crossprod(mobilitybasismat),
                   crossprod(mobilitybasismat,mobility.matrix)) #calculates the coefficient vector for the basis expansion using the matrix expression on page 16 in my report.
mobilityfd = fd(mobilitycoef, mobilitybasis,
              list("Day", "County", "Mobility"))
plot(mobilityfd, lty=1, lwd=2, col=1)
plotfit.fd(mobility.matrix[,13], (0:(180-1)), mobilityfd[13], #the index 13 represents we want the 13th county which is county durham in our variable "counties"
           lty=1, lwd=1,cex.pch=0.8,ylim=c(-65,25))
abline(v = 37,lty = 5, col = "red")
axis(1, at=37,cex.axis = 0.7, tck = -0.03) #label on axis at x=37
abline(v = 30,lty = 3, col = "orange")
axis(1, at=30, cex.axis = 0.7, tck = -0.03)
abline(h = 0,lty=3)
}

#this curve plots, side by side, a least squares smoothed and roughness penalty smoothed of the county durham mobility curve
#once again saves the results into a jpeg in R's current working directory
jpeg("least squares vs roughness penalty.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.3,cex.lab=1.5,mfrow=c(1,2),mar = c(4.2,4.3,1,0.5))
smooth.squares() #plots the mobility curve, instead using least squares smoothing with 90 basis functions
mobility.plot(start.region = 13,end.region = 13,total.days = 180, nbasis = 180, norder = 4, lambda = 1)
dev.off()




