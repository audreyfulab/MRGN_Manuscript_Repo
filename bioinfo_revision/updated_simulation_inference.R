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


# a function to bootstrap edge probabilities for a single trio
boostrap_edge_probabilities <- function(trio, number_of_samples=1000) {
    results <- replicate(number_of_samples, {
        sampled_trio <- trio[sample(nrow(trio), replace=TRUE), ]
        result <- MRGN::infer.trio(
            trio = sampled_trio,
            use.perm = FALSE,
        )
        return(result$edge.probabilities)
    })
}

# This function applies the MRGN method to a single dataset
apply.mrgn <- function(data, use.perm=FALSE) {
    result <- MRGN::infer.trio(
        trio = data,
        use.perm = use.perm,

    )
}





#return to root directory
setwd(root)