
straight = cbind(seq(0,20,0.5),seq(0,20,0.5)) #a simulation of points in a line with gradient 1 (the points shown on the left of figure 5.5)
wiggly = cbind(seq(0,20,0.5),seq(0,20,0.5) + c(rnorm(41,0,2))) #adding some random error to these points to produce the set of points seen on the right of figure 5.5

plot2d.pca<-function(Z,add.proj=FALSE,...){ #this function plots the first principal component (2D) of a data matrix for two variables
  
  #"add.proj = FALSE" is a reference to part of the function, defined later,
  #which we may or may not want to call.
  #The dots mean any graphics parameters we enter when calling the function
  #will be passed along to the plot function.
  #Z will be the data set we are generating PCs for
  
  
  n<-dim(Z)[1]
  plot(Z,...)
  eigen.var<-eigen(var(Z))
  
  #Note that eigen() finds both eigenvectors and eigenvalues, hence subscript use below
  
  gamma1<-eigen.var[[2]][,1]
  lambda1<-eigen.var[[1]][1]
  pcline1a<-colMeans(Z)-2*sqrt(lambda1)*gamma1
  pcline1b<-colMeans(Z)+2*sqrt(lambda1)*gamma1
  
  #Here we are making an arbitrary decision
  #to draw lines with a length of four SDs of projected data
  
  segments(pcline1a[1],pcline1a[2],pcline1b[1],pcline1b[2],col=3,lwd=2)
  if(add.proj){
    t1<-t(gamma1)%*%(t(Z)-colMeans(Z))  
    for(i in 1:n){
      segments(Z[i,1],Z[i,2],colMeans(Z)[1]+t1[i]*gamma1[1],colMeans(Z)[2]+t1[i]*gamma1[2],col=3,lty=2)  
    }
  }
}

#this remaining code runs the pca plot function on the two simulated data sets at the top of this file. Producing figure 5.5
jpeg("continuous registration pca example.jpeg", units = "in", width = 9, height = 4, res = 400)
par(cex.axis=1.3,cex.lab=1.5,mfrow=c(1,2),mar = c(4,4.3,1,0.5))
plot2d.pca(straight,add.proj=TRUE,asp=1,xlim=c(-0.5,20.5),ylim=c(-1,21),xlab=expression("x"[0]*"(t)"),ylab="x[h(t)]")
plot2d.pca(wiggly,add.proj=TRUE,asp=1,xlim=c(-0.5,20.5),ylim=c(-1,21),xlab=expression("x"[0]*"(t)"),ylab="")
dev.off()

