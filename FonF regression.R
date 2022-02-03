library("plot3D")
data(DTI)
DTI
DTI1 <- DTI[DTI$visit==1 & complete.cases(DTI),]

par(mfrow=c(1,2))

# Fit model with linear functional term for CCA
fit.lf <- pfr(pasat ~ ff(cca, k=30, bs="ps"), data=DTI1)
plot(fit.lf, ylab=expression(paste(beta(t))), xlab="t")



# Y(t) = f(t) + \int X1(s)\beta(s,t)ds + eps
set.seed(2121)
data1 <- pffrSim(scenario="ff", n=40)
t <- attr(data1, "yindex")
s <- attr(data1, "xindex")
m1 <- pffr(Y ~ ff(X1, xind=s), yind=t, data=data1)
summary(m1)
plot(m1, scheme = 2,pages=1)


mortality.matrix2 <- t(as.matrix(mortality))
mobility2 <- read.csv("/Users/josephpiekos/Desktop/project III/mobility_matrix.csv")
mobility.matrix2 <- t(as.matrix(mobility2))

m_ff <- pffr(mortality.matrix2 ~ ff(mobility.matrix2,xind=seq(1,181,l=181))) #creates regression fit
summary(m_ff)
psi_plot <- plot(m_ff, select = 2, pers=TRUE)[[2]]
psi_plot
psi_plot$x
psi_plot$y
psi_plot$fit


jpeg("FonF regression plane.jpeg", units = "in", width = 12, height = 8 , res = 400)
par(cex.lab=1.5,cex.axis=1,mar=c(1,1,3,1)+.1,cex.main = 2.5,bg = "#EEF0F2")
persp(psi_plot$x, psi_plot$y, matrix(psi_plot$fit, 40, 40),
      xlab = "Mobility time", ylab = "Mortality time", r=3, phi=40, theta = 40, ticktype="detailed", main=expression(hat(beta)(t,s)), zlab = "", border= NA, col="orangered", shade = 0.4)
dev.off()
