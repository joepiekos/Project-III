#the code in this file is used to create figure 5.1
#it is a simple example of phase vs amplitude variation where sine functions are transform appropriately

#set up the argument values over the interval [-pi,3pi]
x <- seq(-pi,3*pi,length.out=100)
#evaluate three sine functions which differ only in amplitude
y1 <- c(sin(x),1.2*sin(x),1.4*sin(x))
#evaluate three sine functions varying only in phase
y2 <- c(sin(x),sin(x-(pi*0.1)),sin(x-(pi*0.2)))

#the remaining code plots each set of 3 functions side by side, to produce figure 5.1
jpeg("simple phase vs amplitude variation.jpeg", units = "in", width = 10, height = 4, res = 400)
par(cex.axis=1.3,cex.lab=1.5,cex.pch=2,mfrow=c(1,2),mar = c(4,2.5,1,0.5))


plot(x,y2[1:100],type="l",ylab="y",ylim=c(min(y1),max(y1)),xlim=c(-0.2,(2*pi+0.2)),main="Phase variation",xlab="t")
lines(x,y2[101:200],type="l")
lines(x,y2[201:300],type="l")
abline(h=0,lty=3)

plot(x,y1[1:100],type="l",ylab="y",ylim=c(min(y1),max(y1)),xlim=c(-0.2,(2*pi+0.2)),main="Amplitude variation",xlab="t")
lines(x,y1[101:200],type="l")
lines(x,y1[201:300],type="l")
abline(h=0,lty=3)
dev.off()

