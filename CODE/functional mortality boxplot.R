#this code creates figure 6.3, the functional boxplot of the 78 mortality curve

jpeg("mortality boxplot.jpeg", units = "in", width = 10, height = 6.67 , res = 400)
par(cex.lab=1.8,cex.axis=1.7,mar=c(5,5,1,1)+.1,cex.main = 2.5)
mortalityfd.boxplot = boxplot.fd(mortalityfd,xlab="Day",ylab="Mortality (Deaths per 100,000 people)") 
mortalityfd.boxplot
legend(100, 3.1, legend = c("Outliers"), col = c("red"), lty = c(2),lwd = c(1.8), cex = 1.9, text.font = 2)
dev.off()