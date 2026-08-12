# This script reruns the simulation inference with updated parameters and settings.
# following reviewer comments


# load libraries
library(MRGN)
library(MRPC)
source("adapted_GMAC_func/GMAC_moded/R/GMAC.R")


root <- getwd()
setwd("./bioinfo_revision/")
#load in simulated datasets:
sim_data <- loadRData(file = "./simulated_data/simulated_trios.RData")


# This function applies the MRGN method to a single dataset
apply.mrgn <- function(data) {

}





#return to root directory
setwd(root)