# mldist_plot.R plot the likelihood distribution of the best model and a set of alternative models with very low gene flow

# Set the path to the main directory
PATH <- "PATH_TO_MAIN_DIRECTORY" # Replace with the path to the runboot2.pl output directory
pair <- "D00_H00" # Replace with the name of the population pair of interest

setwd(PATH)
pdf("mldist.pdf", 15, 10)

maxL <- read.table(paste(pair, ".maxL.lhoods", sep=""), sep="\t")[,2]
A <- read.table(paste(pair, ".A.lhoods", sep=""), sep="\t")[,2]
B <- read.table(paste(pair, ".B.lhoods", sep=""), sep="\t")[,2]
C <- read.table(paste(pair, ".C.lhoods", sep=""), sep="\t")[,2]
  
if (min(maxL)<max(A,B,C)) {
  boxplot(maxL, A, B, C, main=paste(pair, " (overlap)"), names=c("maxL", "A", "B", "C"), 
          ylab="Log-likelihood", cex.lab=2, cex.axis=1.9, cex.main=2, xlab="Model")
} else {boxplot(maxL, A, B, C, main=pair, names=c("maxL", "A", "B", "C"), 
         ylab="Log-likelihood", cex.lab=2, cex.axis=1.8, cex.main=2, xlab="Model")}	

dev.off()
