# source the 'REST' R file
source("REST.r")
# import data and specify target and suspect
# the data provided contains 10 suspects and 10 replicates of target sequentially 
data <- read.csv("ADE_otu.csv")
target = 11 #column number for the target/sink/evidence
suspect = 1:10 #column numbers for the sources/suspects 
data_target = data[,target]
data_suspect = data[,suspect]
# run REST for 1 target 1 time
pro <- REST(data_target,data_suspect,nbootstrap = 10)
# show proportions 
pro
# plot a graph to show the proportion
# the true proportion of the sample data set is 
prowithtruth <- rbind(c(0.5,0.3,0.2,0), pro)
rownames(prowithtruth) <- c("Truth", "REST")
barplot(prowithtruth, main="Simple Barchart to Compare the Propotions",
        col=c("blue","orange"),
        legend = rownames(prowithtruth), beside=TRUE)

