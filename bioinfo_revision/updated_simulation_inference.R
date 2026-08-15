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
    probabilities <- vector("list", length=number_of_samples)
    for (i in 1:number_of_samples) {
        sampled_trio <- trio[sample(nrow(trio), replace=TRUE), ]
        result <- MRGN::infer.trio(
            trio = sampled_trio,
            use.perm = FALSE
        )
        probabilities[[i]] <- MRGN::get.adj(result)
    }
    return(Reduce("+", probabilities) / number_of_samples)
}

# This function applies the MRGN method to a single dataset
apply.mrgn <- function(data, use.perm=FALSE) {
    result <- MRGN::infer.trio(
        trio = data,
        use.perm = use.perm,
        verbose = TRUE
    )
    inferred_model <- result$Inferred.Model
    inferred_adj <- MRGN::get.adj(result)
    final_result <- list(model=inferred_model, 
                         adj=inferred_adj,
                         edge.probabilities=boostrap_edge_probabilities(
                            data, 
                            number_of_samples=1000
                            ),
                            details = result
                         )
    return(final_result)
}





#return to root directory
setwd(root)