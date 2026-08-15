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


# `standardize` selects the scale the effects are reported on:
#   FALSE  (default) raw slopes, in (response units) / (allele count)
#   TRUE             b * sd(V) / sd(Y), the partial-correlation scale
#   "both"           a two column data.frame with one scale per column, so the combined
#                    plot gets both without refitting every model twice
select.effect.scale <- function(raw, standardized, standardize) {
    if (identical(standardize, "both")) {
        return(data.frame(raw = unname(raw), standardized = unname(standardized)))
    }
    if (isTRUE(standardize)) {
        return(standardized)
    }
    return(raw)
}

empty.effects <- function(standardize) {
    # skipped trios have to come back in the same shape select.effect.scale() would
    # have returned, otherwise the rbind()/unlist() in the collector chokes on them
    if (identical(standardize, "both")) {
        return(data.frame(raw = numeric(0), standardized = numeric(0)))
    }
    return(numeric(0))
}


get.snp.effects <- function(data, target, standardize = FALSE, adjust = TRUE,
                            pcs.only = TRUE, verbose = FALSE) {
    # Get the effect of the SNP (column 1) on the specified target gene.
    # Unlike the PC version this returns a single effect per trio, since there is only
    # one genotype column.
    if (target == "cis") {
        Y <- data[, 2]
    } else if (target == "trans") {
        Y <- data[, 3]
    } else {
        stop("target must be one of 'cis', 'trans'")
    }

    V <- data[, 1]

    # NOTE: the other gene is deliberately NOT in the model. The simulation generates
    # T1 = b.snp * V + (confounding) + e and T2 = b.med * T1 + ..., so the parameter
    # these effects calibrate is the V -> gene effect adjusted for confounding only.
    # Conditioning on the other gene would adjust away a mediator (for T1 -> T2) or
    # condition on a collider, and estimates a different quantity than b.snp.
    # That conditional version is MRGN's b11/b21 test, not what is wanted here.
    covs <- data[, -c(1:3), drop = FALSE]

    # keep only the PC columns: pcr/platform/sex are binary known covariates, so their
    # slopes live on a completely different scale than the PC scores (SD ~ 12-25)
    if (pcs.only) {
        covs <- covs[, grepl("^PC[0-9]+$", colnames(covs)), drop = FALSE]
    }

    # drop degenerate (numerically constant) covariate columns. na.rm matters: without
    # it a single missing value makes sd() return NA, is.finite(NA) is FALSE, and a
    # perfectly good covariate gets silently dropped.
    cov.sd <- apply(covs, 2, sd, na.rm = TRUE)
    keep <- is.finite(cov.sd) & cov.sd > SD.TOL
    if (verbose && any(!keep)) {
        message("  dropping degenerate covariate(s): ",
                paste(colnames(covs)[!keep], collapse = ", "))
    }
    covs <- covs[, keep, drop = FALSE]

    # Assemble the design first and restrict to complete cases. The genotype column
    # carries missing values in some trios (MRGN's own helpers na.omit V before using
    # it), which makes sd() return NA and every comparison against it an error. lm()
    # would drop these rows anyway, so doing it here means the slope and the SDs used
    # to standardize it are computed on exactly the same rows.
    if (adjust && ncol(covs) > 0) {
        design <- data.frame(Y = Y, V = V, covs, check.names = FALSE)
    } else {
        design <- data.frame(Y = Y, V = V)
    }
    design <- design[complete.cases(design), , drop = FALSE]

    # a monomorphic genotype column carries no information about the gene, and it is
    # also collinear with the intercept, so lm() would alias it away
    if (nrow(design) < 3 || sd(design$V) <= SD.TOL || sd(design$Y) == 0) {
        if (verbose) {
            message("  skipping trio: degenerate genotype or response, or < 3 complete cases")
        }
        return(empty.effects(standardize))
    }

    # Fit the gene on the SNP, adjusting for the confounders (PCs)
    model <- lm(Y ~ ., data = design)

    # a saturated fit has no residual df, so the SE and hence the slope are unusable
    if (model$df.residual < 1) {
        if (verbose) {
            message("  skipping trio: no residual degrees of freedom")
        }
        return(empty.effects(standardize))
    }

    # rank deficiency comes back as NA. Index by name rather than position: if an
    # earlier column were aliased, coef() would still be named but a positional
    # index would silently pick up a different predictor.
    snp.effect <- coef(model)["V"]
    if (is.na(snp.effect)) {
        if (verbose) {
            message("  skipping trio: genotype aliased by rank deficiency")
        }
        return(empty.effects(standardize))
    }

    # The raw slope is in units of (response units) / (allele count), and the responses
    # are residualized expression whose SD ranges over 8 orders of magnitude across
    # trios (4e-04 to 3.8e+04), so raw slopes are not comparable trio-to-trio.
    # b * sd(V) / sd(Y) puts every trio on the common (partial-correlation) scale.
    # Both SDs come from `design`, i.e. the same complete cases the slope was fit on.
    # Both scales are computed here; `standardize` only decides which is returned.
    std.effect <- snp.effect * sd(design$V) / sd(design$Y)

    if (verbose) {
        print(summary(model))
    }
    return(select.effect.scale(snp.effect, std.effect, standardize))
}

collect.snp.effects <- function(data.with.pcs, target, standardize = FALSE,
                                adjust = TRUE, pcs.only = TRUE, verbose = FALSE) {
    # Collect the effect of the SNP on the specified target gene across all datasets
    snp.effects <- list()
    for (i in 1:length(data.with.pcs)) {
        # Get the effect of the SNP on the target gene for each dataset
        snp.effects[[i]] <- get.snp.effects(data = data.with.pcs[[i]], target = target,
                                            standardize = standardize,
                                            adjust = adjust,
                                            pcs.only = pcs.only, verbose = verbose)
    }
    # Combine the effects: a single vector for one scale, a two column data.frame when
    # both scales were requested
    if (identical(standardize, "both")) {
        return(do.call(rbind, snp.effects))
    }
    return(unlist(snp.effects, use.names = FALSE))
}


summarize.effects <- function(effects, target) {
    # numeric summary of the bounds we are after. One effect per trio here, so any
    # shortfall against length(data.with.pcs) is trios that were skipped above.
    # Recurse over the columns when both scales are present, since the two live on
    # completely different scales.
    if (is.data.frame(effects)) {
        summarize.effects(effects$raw, paste0(target, " (raw slopes)"))
        summarize.effects(effects$standardized, paste0(target, " (standardized)"))
        return(invisible(NULL))
    }
    cat("\n--- ", target, ": ", length(effects), " SNP effects from ",
        length(data.with.pcs), " trios ---\n", sep = "")
    if (length(effects) < length(data.with.pcs)) {
        cat("(", length(data.with.pcs) - length(effects),
            " trio(s) skipped: degenerate genotype, no residual df, or aliased)\n", sep = "")
    }
    print(round(quantile(effects, c(0, 0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995, 1)), 4))
    cat("mean = ", round(mean(effects), 4),
        ", sd = ", round(sd(effects), 4),
        ", max |effect| = ", round(max(abs(effects)), 4), "\n", sep = "")
}


fd.binwidth <- function(x) {
    # Freedman-Diaconis bin width, 2 * IQR / n^(1/3). Driven by the IQR rather than the
    # range, so a few extreme slopes cannot dictate the bin size for everything else.
    # geom_histogram() accepts a function here and calls it once per facet, which is
    # what we want: the two scales need completely different widths.
    rng <- diff(range(x, na.rm = TRUE))
    if (!is.finite(rng) || rng <= 0) {
        return(1)
    }
    bw <- 2 * IQR(x, na.rm = TRUE) / length(x)^(1/3)
    if (!is.finite(bw) || bw <= 0) {
        bw <- rng / 60
    }
    # keep the bin count in a readable 20-300, since FD asks for thousands of bins when
    # the window is still wide relative to the IQR, and for a handful when it is not
    return(min(max(bw, rng / 300), rng / 20))
}


plot.distribution <- function(effects, target, save = FALSE, filename = NULL,
                              trim = 0.005) {
    # `effects` is the two column data.frame from collect.snp.effects(standardize = "both")
    if (!is.data.frame(effects)) {
        stop("effects must be the data.frame returned by ",
             "collect.snp.effects(..., standardize = 'both')")
    }

    # Stack the two scales into one long frame with a panel label
    scale.labels <- c("raw slope, b", "standardized, b * sd(V) / sd(Y)")
    effects.df <- rbind(
        data.frame(effects = effects$raw, scale = scale.labels[1]),
        data.frame(effects = effects$standardized, scale = scale.labels[2])
    )
    effects.df$scale <- factor(effects.df$scale, levels = scale.labels)

    # Restrict each panel to its central quantile window before binning. Bin width
    # alone cannot rescue the raw panel: a few extreme slopes stretch the range so far
    # that every other effect lands in the first bin no matter how the bins are chosen.
    # Done per panel because the two scales have completely different tails, and the
    # count that falls outside is reported in the subtitle rather than hidden.
    effects.df <- do.call(rbind, lapply(split(effects.df, effects.df$scale), function(d) {
        if (trim <= 0) {
            d$outside <- FALSE
        } else {
            lims <- quantile(d$effects, c(trim, 1 - trim), na.rm = TRUE)
            d$outside <- d$effects < lims[1] | d$effects > lims[2]
        }
        return(d)
    }))
    n.outside <- sum(effects.df$outside)
    plot.df <- effects.df[!effects.df$outside, , drop = FALSE]

    trim.note <- ""
    if (n.outside > 0) {
        trim.note <- paste0("; central ", round(100 * (1 - 2 * trim), 1),
                            "% of each panel shown, ", n.outside, " outside the axes")
    }

    # Create the histogram with density overlay
    p <- ggplot(plot.df, aes(x = effects)) +
        geom_histogram(aes(y = after_stat(density)), binwidth = fd.binwidth,
                       fill = "lightblue", color = "black") +
        geom_density(color = "red", linewidth = 1) +
        # scales = "free": the raw slopes and the standardized effects differ by orders
        # of magnitude, so a shared axis would compress one panel into a single bar
        facet_wrap(~ scale, ncol = 1, scales = "free") +
        labs(title = paste("Distribution of SNP Effects on", target),
             subtitle = paste0(nrow(effects), " effects over ",
                               length(data.with.pcs), " trios", trim.note),
             x = "Effect Size",
             y = "Density") +
        theme_minimal()

    # Save the plot if requested
    if (save) {
        if (is.null(filename)) {
            # If no filename is provided, create a default one
            path = getwd()
            filename = paste0(path, "/SNP_effects_distribution_", target, ".png")
        }
        # Save the plot to the specified filename. Taller than the single panel version
        # so the stacked facets each get room.
        ggsave(filename, plot = p, width = 7, height = 8, dpi = 150)
    }
    return(p)
}


# run for each target variable: cis, trans
# standardize = "both" so the summary and the figure cover the raw and the standardized
# scale from one pass of the regressions

# cis gene distribution of SNP effects:
print("Collecting SNP effects for cis genes...")
cis.effects <- collect.snp.effects(data.with.pcs, target = "cis", standardize = "both")
summarize.effects(cis.effects, "cis")
plot.distribution(cis.effects, target = "cis", save = TRUE, filename = "SNP_effects_distribution_cis.png")

print("Collecting SNP effects for trans genes...")
trans.effects <- collect.snp.effects(data.with.pcs, target = "trans", standardize = "both")
summarize.effects(trans.effects, "trans")
plot.distribution(trans.effects, target = "trans", save = TRUE, filename = "SNP_effects_distribution_trans.png")



setwd(root)
