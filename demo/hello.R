library(FSToolboxR)
# Read CSV into R
datasetFile <- system.file("extdata", "test.csv", package="FSToolboxR")
dataset <- read.csv(datasetFile, header=TRUE, sep=",")
class=matrix(unlist(dataset[1]), ncol = 1, byrow = TRUE)
features=matrix(unlist(dataset[,2:ncol(dataset)]), ncol = ncol(dataset)-1, byrow = TRUE)
FSToolboxR("CMIM",10,features,NULL,class,0.0,0.0)