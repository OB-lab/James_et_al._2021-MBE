# summaryplots_fsc.R visually summarises the performance of fastsimcoal2 runs across different models

# Set the path to the main directory
PATH <- "PATH_TO_MAIN_DIRECTORY" # Replace with the path to the results directory output by extract_ml.pl Perl script
pair <- "D00_H00" # Replace with the name of the population pair of interest

# Remove the first and last lines of the brent_lhoods files and save it in a new file
setwd(PATH)
for (k in 1:10) { # Do this action in all models per population pair. The number of models can be adjusted
  m <- read.table(paste(pair, "_model", k, ".brent_lhoods", sep=""), sep='\n')
  write.table(m[c(3:nrow(m)-1),], paste(pair, "_model", k, ".brent_lhoods.txt", sep=""), quote=F, row.names=F, col.names=F)
}

# Save the plots in a PDF file
pdf(paste("summary_fsc_", pair, ".pdf", sep=""), width=15, height=10)

# Summarise the likelihood values of all the runs per model as boxplots
allruns <- read.table(paste(pair, "_bestlhoods.txt", sep=""), sep="\t") # Verify all rows have equal number of columns, complete with NA if needed
boxplot(as.numeric(allruns[1,seq(2, ncol(allruns)-1, 1)]), # Exclude the frist and last columns because they are not likelihood values
        as.numeric(allruns[2,seq(2, ncol(allruns)-1, 1)]),
        as.numeric(allruns[3,seq(2, ncol(allruns)-1, 1)]),
        as.numeric(allruns[4,seq(2, ncol(allruns)-1, 1)]),
        as.numeric(allruns[5,seq(2, ncol(allruns)-1, 1)]),
        as.numeric(allruns[6,seq(2, ncol(allruns)-1, 1)]),
        as.numeric(allruns[7,seq(2, ncol(allruns)-1, 1)]),
        as.numeric(allruns[8,seq(2, ncol(allruns)-1, 1)]),
        as.numeric(allruns[9,seq(2, ncol(allruns)-1, 1)]),
        as.numeric(allruns[10,seq(2, ncol(allruns)-1, 1)]),
        main=pair, ylab="Likelihood", xlab="Demographic model",
        names=c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10")) # Adjust accoding to the number of models

# Plot the likelihood values of the best run per model along the optimisation cycles
library("ggplot2")
library("gridExtra")
m <- list()
for (w in 1:10) { # Number of models per pair
  m[[w]] <- read.table(paste(pair, "_model", w, ".brent_lhoods.txt", sep=""), sep='\t')
}
grid.arrange(qplot(c(1:nrow(m[[1]])), m[[1]][,length(m[[1]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair, " - M1", sep="")),
             qplot(c(1:nrow(m[[2]])), m[[2]][,length(m[[2]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair, " - M2", sep="")),
             qplot(c(1:nrow(m[[3]])), m[[3]][,length(m[[3]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair, " - M3", sep="")),
             qplot(c(1:nrow(m[[4]])), m[[4]][,length(m[[4]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair, " - M4", sep="")),
             qplot(c(1:nrow(m[[5]])), m[[5]][,length(m[[5]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair, " - M5", sep="")),
             qplot(c(1:nrow(m[[6]])), m[[6]][,length(m[[6]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair, " - M6", sep="")),
             qplot(c(1:nrow(m[[7]])), m[[7]][,length(m[[7]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair, " - M7", sep="")),
             qplot(c(1:nrow(m[[8]])), m[[8]][,length(m[[8]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair, " - M8", sep="")),
             qplot(c(1:nrow(m[[9]])), m[[9]][,length(m[[9]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair, " - M9", sep="")),
             qplot(c(1:nrow(m[[10]])), m[[10]][,length(m[[10]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair, " - M10", sep="")),
             nrow=3)


dev.off()
