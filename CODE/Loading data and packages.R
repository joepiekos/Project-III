#the following two packages provide a large range of commands for working with functional data
install.packages("fda")
install.packages("refund")

library(fda)
library(refund)

#the next three commands load in the data sets for mortality, mobility and positivity
#the data can be downloaded from the same github respository, as .csv files. 
#File path urls will obviously need to be changed for each user.

mobility <- read.csv("/Users/josephpiekos/Desktop/project III/data/mobilitymatrix.csv")
mortality <- read.csv("/Users/josephpiekos/Desktop/project III/data/english_mortality_final.csv")
positivity <- read.csv("/Users/josephpiekos/Desktop/project III/data/english_positivity_final.csv")

