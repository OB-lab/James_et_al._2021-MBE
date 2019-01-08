# summaryplots_fsc.R graphically summarises the performance of fastsimcoal runs across different models and populations

# Pairs names, which should match the corresponding directory names
pair <- c("D00_H00", "D01_H01", "D03_H02", "D04_H05", "D05_H06", "D12_H14", "D14_H15")
# Triads names, which should match the corresponding directory names
triad <- c("D32_H12_H12A")

# Remove the first and last lines of the brent_lhoods files and save it in a new file
for (i in 1:length(pair)) { # Do this action in all pairs
  setwd(paste("PATH", # Replace for the PATH to the directory that contains the pairs directories
              pair[i], "/results_", pair[i], sep="")) 
  for (k in 1:7) { # Do this action in all models per population pair. The number of models can be adjusted
    m <- read.table(paste(pair[i], "_model", k, ".brent_lhoods", sep=""), sep='\n')
    write.table(m[c(3:nrow(m)-1),], paste(pair[i], "_model", k, ".brent_lhoods.txt", sep=""), quote=F, row.names=F, col.names=F)
  }
}  
for (i in 1:length(triad)) { # Do this action for all triads
  setwd(paste("PATH", # Replace for the PATH to the directory that contains the triads directories
              triad[i], "/results_", triad[i], sep=""))
  for (k in 1:9) { # Do this action in all models per population triad. The number of models can be adjusted
    m <- read.table(paste(triad[i], "_model", k, ".brent_lhoods", sep=""), sep='\n')
    write.table(m1[c(3:nrow(m1)-1),], paste(triad[i], "_model", k, ".brent_lhoods.txt", sep=""), quote=F, row.names=F, col.names=F)
  }
}  

# Save the plots in a PDF file
setwd("PATH") # Replace for the PATH to the directory that contains the triads directories
pdf("SummaryPairRuns.pdf", width=10, height=12)

# Summarise the likelihood values of all the runs per model as boxplots
for (i in 1:length(pair)) { # Do this action for all the pairs
  setwd(paste("~/Dropbox/PhD_OB/SideProjects/PacBioMapping/", pair[i], "/results_", 
              pair[i], sep="")) # fastsimcoal output files should be contained in a directory named results_POP1_POP2 within the directory POP1_POP2
  allruns <- read.table(paste(pair[i], "_bestlhoods.txt", sep=""), sep="\t") #Verify all rows have equal number of columns, complete with NA if needed
  boxplot(as.numeric(allruns[1,seq(2, ncol(allruns)-1, 1)]), # Exclude the frist and last columns because they are not likelihood values
          as.numeric(allruns[2,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[3,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[4,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[5,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[6,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[7,seq(2, ncol(allruns)-1, 1)]),
          main=pair[i], ylab="Likelihood", xlab="Demographic model",
          names=c("M1","M2","M3","M4","M5","M6","M7"))
}
for (i in 1:length(triad)) { # Do this action for all the triads
  setwd(paste("~/Dropbox/PhD_OB/SideProjects/PacBioMapping/", triad[i], "/results_", 
              triad[i], sep="")) # fastsimcoal output files should be contained in a directory named results_POP1_POP2 within the directory POP1_POP2
  allruns <- read.table(paste(triad[i], "_bestlhoods.txt", sep=""), sep="\t") #Verify all rows have equal number of columns
  boxplot(as.numeric(allruns[1,seq(2, ncol(allruns)-1, 1)]), # Exclude the frist and last columns because they are not likelihood values
          as.numeric(allruns[2,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[3,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[4,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[5,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[6,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[7,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[8,seq(2, ncol(allruns)-1, 1)]),
          as.numeric(allruns[9,seq(2, ncol(allruns)-1, 1)]),
          main=pair[i], ylab="Likelihood", xlab="Demographic model",
          names=c("M1","M2","M3","M4","M5","M6","M7", "M8", "M9"))
}

# Plot the likelihood values of the best run per model along the optimisation cycles
library("ggplot2")
library("gridExtra")
for (i in 1:length(pair)) { # Do this action for all the pairs
  setwd(paste("~/Dropbox/PhD_OB/SideProjects/PacBioMapping/", pair[i], "/results_", pair[i], sep=""))
  m <- list()
  for (w in 1:7) { # Number of models per pair
    m[[w]] <- read.table(paste(pair[i], "_model", w, ".brent_lhoods.txt", sep=""), sep='\t')
  }
  grid.arrange(qplot(c(1:nrow(m[[1]])), m[[1]][,length(m[[1]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M1", sep="")),
               qplot(c(1:nrow(m[[2]])), m[[2]][,length(m[[2]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M2", sep="")),
               qplot(c(1:nrow(m[[3]])), m[[3]][,length(m[[3]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M3", sep="")),
               qplot(c(1:nrow(m[[4]])), m[[4]][,length(m[[4]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M4", sep="")),
               qplot(c(1:nrow(m[[5]])), m[[5]][,length(m[[5]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M5", sep="")),
               qplot(c(1:nrow(m[[6]])), m[[6]][,length(m[[6]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M6", sep="")),
               qplot(c(1:nrow(m[[7]])), m[[7]][,length(m[[7]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M7", sep="")),
               nrow=3)
}
for (i in 1:length(triad)) { # Do this action for all the triads
  setwd(paste("~/Dropbox/PhD_OB/SideProjects/PacBioMapping/", triad[i], "/results_", triad[i], sep=""))
  m <- list()
  for (w in 1:9) { # Number of models per triad
    m[[w]] <- read.table(paste(triad[i], "_model", w, ".brent_lhoods.txt", sep=""), sep='\t')
  }
  grid.arrange(qplot(c(1:nrow(m[[1]])), m[[1]][,length(m[[1]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M1", sep="")),
               qplot(c(1:nrow(m[[2]])), m[[2]][,length(m[[2]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M2", sep="")),
               qplot(c(1:nrow(m[[3]])), m[[3]][,length(m[[3]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M3", sep="")),
               qplot(c(1:nrow(m[[4]])), m[[4]][,length(m[[4]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M4", sep="")),
               qplot(c(1:nrow(m[[5]])), m[[5]][,length(m[[5]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M5", sep="")),
               qplot(c(1:nrow(m[[6]])), m[[6]][,length(m[[6]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M6", sep="")),
               qplot(c(1:nrow(m[[7]])), m[[7]][,length(m[[7]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M7", sep="")),
               qplot(c(1:nrow(m[[8]])), m[[8]][,length(m[[8]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M8", sep="")),
               qplot(c(1:nrow(m[[9]])), m[[9]][,length(m[[9]])], geom="path", xlab="EMC cycle", ylab="Likelihood", main=paste(pair[i], " - M9", sep="")),
               nrow=3)
}

dev.off()
setwd("~/Dropbox/PhD_OB/SideProjects/PacBioMapping/") # Go back to the main directory
