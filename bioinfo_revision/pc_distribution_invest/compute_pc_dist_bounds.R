library(MRGN)
library(ggplot2)

root <- getwd()
data.with.pcs <- loadRData("./GTEx/data/data.with.PCs.WholeBlood.RData")

setwd("./bioinfo_revision/pc_distribution_invest")
# Column layout of each element of data.with.pcs (see GTEx/data/PC_LRNA_PC_Selection_manu.R):
#   1        : SNP genotype (V)
#   2        : cis gene expression (T1)
#   3        : trans gene expression (T2)
#   4:6      : known covariates (pcr, platform, sex)
#   7:ncol   : the PCs selected for that trio (named "PC<k>")

# minimum SD a covariate column must have to be kept in the regression. PC670 (= PC_n,
# the null direction left over after centering 670 samples) is retained by the PC
# selection step in 431 trios with SD ~1e-13. It is not exactly constant, so lm() does
# not alias it away; instead its slope explodes to ~1e11 and dominates everything else.
SD.TOL <- 1e-08


get.conf.effects <- function(data, target, standardize = TRUE,
                             pcs.only = TRUE, verbose = FALSE) {
    # Get the effects of confounders (PCs) on the specified target variable
    if( target == "genotype") {
        Y <- data[,1]
    } else if (target == "cis") {
        Y <- data[,2]
    } else if (target == "trans") {
        Y <- data[,3]
    } else {
        stop("target must be one of 'genotype', 'cis', 'trans'")
    }

    covs <- data[, -c(1:3), drop = FALSE]

    # keep only the PC columns: pcr/platform/sex are binary known covariates, so their
    # slopes live on a completely different scale than the PC scores (SD ~ 12-25) and
    # pooling them into the same histogram is comparing apples to oranges
    if (pcs.only) {
        covs <- covs[, grepl("^PC[0-9]+$", colnames(covs)), drop = FALSE]
    }

    # drop degenerate (numerically constant) covariate columns
    cov.sd <- apply(covs, 2, sd)
    keep <- is.finite(cov.sd) & cov.sd > SD.TOL
    if (verbose && any(!keep)) {
        message("  dropping degenerate covariate(s): ",
                paste(colnames(covs)[!keep], collapse = ", "))
    }
    covs <- covs[, keep, drop = FALSE]
    cov.sd <- cov.sd[keep]

    if (ncol(covs) == 0 || sd(Y) == 0) {
        return(numeric(0))
    }

    # Fit a linear model with the target variable as the response and the confounders (PCs) as predictors
    model <- lm( Y ~ ., data = covs)
    conf.effects <- coef(model)[-1]
    # aliased predictors (rank deficiency) come back as NA
    conf.effects <- conf.effects[!is.na(conf.effects)]

    # Put every trio on a common scale. The raw slope is in units of
    # (response units) / (PC score units), and the responses are residualized
    # expression whose SD ranges over 8 orders of magnitude across trios
    # (4e-04 to 3.8e+04), so raw slopes are not comparable trio-to-trio.
    # b * sd(PC) / sd(Y) is the standardized (partial-correlation scale) effect.
    if (standardize) {
        conf.effects <- conf.effects * cov.sd[names(conf.effects)] / sd(Y)
    }

    if (verbose) {
        print(summary(model))
    }
    return(conf.effects)
}

collect.pc.effects <- function(data.with.pcs, target, standardize = TRUE,
                               pcs.only = TRUE, verbose = FALSE) {
    # Collect the effects of confounders (PCs) on the specified target variable across all datasets
    pc.effects <- list()
    for (i in 1:length(data.with.pcs)) {
        # Get the effects of confounders on the target variable for each dataset
        pc.effects[[i]] <- get.conf.effects(data = data.with.pcs[[i]], target = target,
                                            standardize = standardize,
                                            pcs.only = pcs.only, verbose = verbose)
    }
    # Combine the effects into a single vector
    return(unlist(pc.effects, use.names = FALSE))
}


summarize.effects <- function(effects, target) {
    # numeric summary of the bounds we are after
    cat("\n--- ", target, ": ", length(effects), " PC effects from ",
        length(data.with.pcs), " trios ---\n", sep = "")
    print(round(quantile(effects, c(0, 0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995, 1)), 4))
    cat("mean = ", round(mean(effects), 4),
        ", sd = ", round(sd(effects), 4),
        ", max |effect| = ", round(max(abs(effects)), 4), "\n", sep = "")
}


plot.distribution <- function(effects, target, save=FALSE, filename=NULL) {
    # Create a data frame for plotting
    effects.df <- data.frame(effects = effects)

    # Create the histogram with density overlay
    p <- ggplot(effects.df, aes(x = effects)) +
        geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "lightblue", color = "black") +
        geom_density(color = "red", linewidth = 1) +
        labs(title = paste("Distribution of PC Effects on", target),
             subtitle = paste0("standardized slopes, b * sd(PC) / sd(Y); ",
                               length(effects), " effects over ", length(data.with.pcs), " trios"),
             x = "Standardized Effect Size",
             y = "Density") +
        theme_minimal()

    # Save the plot if requested
    if (save) {
        if (is.null(filename)) {
            # If no filename is provided, create a default one
            path = getwd()
            filename = paste0(path, "/PC_effects_distribution_", target, ".png")
        }
        # Save the plot to the specified filename
        ggsave(filename, plot = p, width = 7, height = 5, dpi = 150)
    }
    return(p)
}


# run for each target variable: cis, trans

# cis gene distribution of confounder effects:
print("Collecting PC effects for cis genes...")
cis.effects <- collect.pc.effects(data.with.pcs, target = "cis")
summarize.effects(cis.effects, "cis")
plot.distribution(cis.effects, target = "cis", save = TRUE, filename = "PC_effects_distribution_cis.png")

print("Collecting PC effects for trans genes...")
trans.effects <- collect.pc.effects(data.with.pcs, target = "trans")
summarize.effects(trans.effects, "trans")
plot.distribution(trans.effects, target = "trans", save = TRUE, filename = "PC_effects_distribution_trans.png")



setwd(root)
