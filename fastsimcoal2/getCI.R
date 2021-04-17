# getCI.R compute confidence intervals from a file summarising bootstrap runs of fastsimcoal2
# This script that the best model was bidirectional secondary contact (model 5) as specified in James et al. 2021

# Set the path to the main directory
PATH <- "PATH_TO_MAIN_DIRECTORY" # Replace with the path to the runboot2.pl output directory
pair <- "D00_H00" # Replace with the name of the population pair of interest

library(rcompanion)
setwd(PATH)

# Read runboot2.pl output table
Table <- read.table(paste(pair, "_boot.txt", sep=""))
ci <- c()

# Extract the lower and upper 95% confidence interval of each parameter
ci[1] <- as.numeric(groupwiseMean(V1 ~ 1, data=Table, conf=0.95, digits=6)[5])
ci[2] <- as.numeric(groupwiseMean(V1 ~ 1, data=Table, conf=0.95, digits=6)[6])
ci[3] <- as.numeric(groupwiseMean(V2 ~ 1, data=Table, conf=0.95, digits=6)[5])
ci[4] <- as.numeric(groupwiseMean(V2 ~ 1, data=Table, conf=0.95, digits=6)[6])
ci[5] <- as.numeric(groupwiseMean(V3 ~ 1, data=Table, conf=0.95, digits=6)[5])
ci[6] <- as.numeric(groupwiseMean(V3 ~ 1, data=Table, conf=0.95, digits=6)[6])
ci[7] <- as.numeric(groupwiseMean(V4 ~ 1, data=Table, conf=0.95, digits=6)[5])
ci[8] <- as.numeric(groupwiseMean(V4 ~ 1, data=Table, conf=0.95, digits=6)[6])
ci[9] <- as.numeric(groupwiseMean(V5 ~ 1, data=Table, conf=0.95, digits=6)[5])
ci[10] <- as.numeric(groupwiseMean(V5 ~ 1, data=Table, conf=0.95, digits=6)[6])
ci[11] <- as.numeric(groupwiseMean(V6 ~ 1, data=Table, conf=0.95, digits=6)[5])
ci[12] <- as.numeric(groupwiseMean(V6 ~ 1, data=Table, conf=0.95, digits=6)[6])
ci[13] <- as.numeric(groupwiseMean(V7 ~ 1, data=Table, conf=0.95, digits=6)[5])
ci[14] <- as.numeric(groupwiseMean(V7 ~ 1, data=Table, conf=0.95, digits=6)[6])

# Parameters named in the same order that specified in the parameters file
colnames(ci) <- c("ANCSIZE_min", "ANCSIZE_max", "H_min", "H_max", "D_min", "D_max", "TDIV_min", "TDIV_max", 
                  "TSEC_min", "TSEC_max", "H2D_min", "H2D_max", "D2H_min", "D2H_max")
write.table(ci, "CI_model5.txt", sep="\t", quote=F, row.names=F)
