library(fda)


regions <- read.csv("/Users/josephpiekos/Desktop/project III/english authorities.csv")
mobility <- read.csv("/Users/josephpiekos/Desktop/project III/english_mobility_final.csv")

region = vector("list",78) #list of regions we are looking at
for (i in 1:78){
  region[[ i]] = regions[i,]
}
region

mobility.plot <- function(start.region,end.region,total.days, mobility.category, K, b.order = 4){

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
mobilityfd = smooth.basis((0:(total.days-1)),mobility.matrix,mobilitybasis)$fd
mobility.gcv = smooth.basis((0:(total.days-1)),mobility.matrix,mobilitybasis)$gcv

mobility.gcv
mobility.matrix
plot(mobilityfd,xlab="Day number",ylab="Percentage change from Baseline")
if (number.regions > 1){
  title("Mobility Curves")
  lines(mobilityfd, lty=1, col = "orange")
}
else {
  lines(mobilityfd, lty=1, col = " dark green", lwd = 2)
}
abline(v = 37,lty = 5, col = "red")
#axis(1, at=37,labels=37,cex.axis = 0.7, tck = -0.06) label on axis at x=37
abline(v = 30,lty = 3, col = "orange")
#axis(1, at=30,labels=30, cex.axis = 0.7, tck = -0.06) label on axis at x = 30
legend(total.days-45, max(mobility.matrix)-5, legend = c("No non-essential contact and travel","First Lockdown"), col = c("orange","red"), lty = c(3,5), cex = 1.3, text.font = 2, title = "Key Announcements from Government")
}

#jpeg("mobility curves.jpeg", units = "in", width = 12, height = 8 , res = 400)
#par(cex.lab=2,cex.axis=1.5,mar=c(5,6,2,1)+.1)
mobility.plot(start.region = 1,end.region = 50, total.days = 90, mobility.category = mobility$grocery_and_pharmacy_percent_change_from_baseline,  K = 19, b.order = 4)
#dev.off()

single.region.mobility <- function(total.days, K, b.order = 4,region.name){
  start.region = match(region.name, region)
  end.region = start.region
  mobility.plot(start.region,end.region,total.days, mobility.category = mobility$grocery_and_pharmacy_percent_change_from_baseline, K, b.order)
  #title(main = paste("Mobility Curve for", region[[start.region]], "(Groceries and Pharmacy)"))#CHANGE DEPENDING ON ACTIVITY
}

#jpeg("mobility curve durham.jpeg", units = "in", width = 12, height = 8 , res = 400)
#par(cex.lab=2,cex.axis=1.5,cex.main = 1.5,mar=c(4.1,5,1,0.5)+.1)
single.region.mobility(180,19,4,"County Durham")
#dev.off()




