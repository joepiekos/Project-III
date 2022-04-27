
#firstly sets up the functional data matrix for use in smoothing
mobility <- read.csv("/Users/josephpiekos/Desktop/project III/data/mobilitymatrix.csv") #loads mobility data as a 181 x 78 data frame, where each column represents a different county.
mobility.matrix <- (matrix(0,180,78)) #creates empty matrix, number of rows = time period, number of columns = number of regions

for (i in 1:180){ #there are 60 rows, as you go down i represents another day
  for (j in 1:78){ #there are 3 columns, each column j is a region
    if(number.regions > 1){
      mobility.matrix[i,j] <- mobility[i,j]
    }
    else {
      mobility.matrix[i,1] <- mobility[i,j]
    }
  }
}

#the next section of code is the GCV calculating section.
loglambda = seq(-1,3,0.1) #sets up the range of log(lambda) values to test using GCV. This range is predetermined by the user to capture a full range of variation.
nlambda = length(loglambda) #number of log(lambda) values we will be testing
gcv.values = rep(NA,nlambda) #sets up an emppty matrix for the GCV values of all the curves

total.days = 180
norder = 4 #order of B-splines
nbasis  = 180 #number of basis functions to use

mobilitybasis = create.bspline.basis(c(0,total.days-1), nbasis, norder)

#although this is all set up to use the mobility data, it can simply be changed to use the positivity/mortality data 
#by loading in the respective data and changing all "positivity" instances to the name of the dataset you want
#for each log(lambda) value in our test interval, smooth the mobility curves and calculate the total GCV.

for (lambda.i in 1:nlambda) {
  cat(paste("log10 lambda =",loglambda[lambda.i],"\n"))
  lambda = 10^loglambda[lambda.i]
  mobilityfdPar = fdPar(mobilitybasis, 2, lambda)
  mobilityfd = smooth.basis((0:(total.days-1)), mobility.matrix,
                          mobilityfdPar)
  gcv.values[lambda.i] = sum(mobilityfd$gcv)
}

#this last block of code creates the plot of GCV values (figure 4.2) and saves it in a jpeg in R's current working directory
jpeg("Mobility GCV for lambda.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=2,cex.axis=1.5,mar=c(4.4,5.3,1,0.5)+.1)
plot(loglambda,gcv.values,ylab=expression("GCV(" ~ lambda ~")"),xlab=expression("log(" ~ lambda ~ ")"),cex=1.5)
lines(loglambda,gcv.values)
abline(v=loglambda[8],lty=3)
dev.off()

