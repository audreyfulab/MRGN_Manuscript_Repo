library(MRGN)
library(ggplot2)
library(patchwork)

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
#   FALSE  (default) raw slopes, in (response units) / (PC score units)
#   TRUE             b * sd(PC) / sd(Y), which here IS cor(PC, Y): see the note at
#                    std.effects below
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


get.conf.effects <- function(data, target, standardize = FALSE,
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
        return(empty.effects(standardize))
    }

    # Fit a linear model with the target variable as the response and the confounders (PCs) as predictors
    model <- lm( Y ~ ., data = covs)
    conf.effects <- coef(model)[-1]
    # aliased predictors (rank deficiency) come back as NA
    conf.effects <- conf.effects[!is.na(conf.effects)]

    # The raw slope is in units of (response units) / (PC score units), and the
    # responses are residualized expression whose SD ranges over 8 orders of magnitude
    # across trios (4e-04 to 3.8e+04), so raw slopes are not comparable trio-to-trio.
    # b * sd(PC) / sd(Y) puts every trio on a common, unit-free scale.
    #
    # That scale is exactly the Pearson correlation cor(PC_j, Y). A standardized
    # regression coefficient equals the marginal correlation whenever the predictors
    # are mutually uncorrelated, and these predictors are principal components, so
    # they are orthogonal by construction (checked: the largest off-diagonal
    # correlation among the retained PC columns of a trio is ~5e-15, and the
    # standardized coefficients match cor(PC_j, Y) to ~2e-15). Note this is the
    # marginal correlation, NOT the partial correlation, which differs by up to ~0.05.
    # Both are computed here; `standardize` only decides which is returned.
    std.effects <- conf.effects * cov.sd[names(conf.effects)] / sd(Y)

    if (verbose) {
        print(summary(model))
    }
    return(select.effect.scale(conf.effects, std.effects, standardize))
}

collect.pc.effects <- function(data.with.pcs, target, standardize = FALSE,
                               pcs.only = TRUE, verbose = FALSE) {
    # Collect the effects of confounders (PCs) on the specified target variable across all datasets
    pc.effects <- list()
    for (i in 1:length(data.with.pcs)) {
        # Get the effects of confounders on the target variable for each dataset
        pc.effects[[i]] <- get.conf.effects(data = data.with.pcs[[i]], target = target,
                                            standardize = standardize,
                                            pcs.only = pcs.only, verbose = verbose)
    }
    # Combine the effects: a single vector for one scale, a two column data.frame when
    # both scales were requested
    if (identical(standardize, "both")) {
        return(do.call(rbind, pc.effects))
    }
    return(unlist(pc.effects, use.names = FALSE))
}


summarize.effects <- function(effects, target) {
    # numeric summary of the bounds we are after. Recurse over the columns when both
    # scales are present, since the two live on completely different scales.
    if (is.data.frame(effects)) {
        summarize.effects(effects$raw, paste0(target, " (raw slopes)"))
        summarize.effects(effects$standardized, paste0(target, " (correlation)"))
        return(invisible(NULL))
    }
    cat("\n--- ", target, ": ", length(effects), " PC effects from ",
        length(data.with.pcs), " trios ---\n", sep = "")
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


fmt.effect <- function(v) {
    # the raw slopes run around 1e-04, so the default fixed notation prints them all as
    # "0.00" and the labels say nothing
    return(formatC(v, format = "e", digits = 2))
}


box.summary <- function(x) {
    # The five numbers geom_boxplot() draws: lower whisker, Q1, median, Q3, upper
    # whisker, with the whiskers pulled in to the most extreme observation within
    # 1.5 * IQR of the hinges. Taken from boxplot.stats() so the labels are guaranteed
    # to describe the box actually plotted rather than a separately computed quantile.
    # Computed on the full vector: nothing is dropped from the numbers, only from the
    # view (see the coord_cartesian() zoom in plot.raw.boxplot()).
    stats <- grDevices::boxplot.stats(x)$stats
    stat.names <- c("lower whisker", "Q1", "median", "Q3", "upper whisker")
    return(data.frame(
        stat = stat.names,
        value = stats,
        label = paste0(stat.names, "\n", fmt.effect(stats)),
        # stagger the labels onto two rows: Q1 and Q3 sit within a hinge width of the
        # median, so all five on one row would overprint each other
        y = c(0.5, -0.5, 0.5, -0.5, 0.5),
        stringsAsFactors = FALSE
    ))
}


plot.raw.boxplot <- function(raw, panel.label) {
    # The raw slopes cluster so tightly at zero relative to their extremes that a
    # histogram of them is a single spike (see the pre-revision figure). A boxplot with
    # its five-number summary printed as text carries the magnitudes that matter here.
    stats <- box.summary(raw)
    n.outside <- sum(raw < stats$value[1] | raw > stats$value[5], na.rm = TRUE)

    # Zoom to the whiskers. coord_cartesian() rather than xlim(): this clips the view
    # only, so the box itself is still computed from every effect.
    span <- diff(range(stats$value))
    if (!is.finite(span) || span <= 0) {
        span <- max(abs(stats$value[c(1, 5)]), 1)
    }
    pad <- 0.12 * span

    caption <- paste0("n = ", format(length(raw), big.mark = ","))
    if (n.outside > 0) {
        caption <- paste0(caption, "; ", format(n.outside, big.mark = ","),
                          " outlier(s) beyond the whiskers, not drawn")
    }

    return(ggplot(data.frame(effects = raw), aes(x = effects, y = 0, group = 1)) +
        # orientation is set explicitly: with a constant numeric y, ggplot's own guess
        # is ambiguous and warns about a continuous y aesthetic
        geom_boxplot(orientation = "y", outlier.shape = NA,
                     fill = "lightblue", color = "black", width = 0.5) +
        geom_text(data = stats, aes(x = value, y = y, label = label),
                  inherit.aes = FALSE, size = 3, lineheight = 0.9) +
        coord_cartesian(xlim = c(stats$value[1] - pad, stats$value[5] + pad),
                        ylim = c(-0.85, 0.85)) +
        scale_x_continuous(labels = scales::label_scientific(digits = 2)) +
        # the y axis only positions the box and its labels, so its numbers are noise
        scale_y_continuous(breaks = NULL) +
        labs(title = panel.label, x = "Effect Size", y = NULL, caption = caption) +
        theme_minimal() +
        theme(plot.title = element_text(size = 10, hjust = 0.5),
              plot.caption = element_text(size = 8, hjust = 0.5)))
}


plot.std.histogram <- function(standardized, panel.label, trim,
                               x.label = "Effect Size") {
    # Restrict to the central quantile window before binning, and report the count that
    # falls outside in the caption rather than hiding it.
    outside <- rep(FALSE, length(standardized))
    if (trim > 0) {
        lims <- quantile(standardized, c(trim, 1 - trim), na.rm = TRUE)
        outside <- standardized < lims[1] | standardized > lims[2]
    }
    n.outside <- sum(outside, na.rm = TRUE)

    caption <- paste0("n = ", format(length(standardized), big.mark = ","))
    if (n.outside > 0) {
        caption <- paste0(caption, "; central ", round(100 * (1 - 2 * trim), 1),
                          "% shown, ", format(n.outside, big.mark = ","),
                          " outside the axes")
    }

    return(ggplot(data.frame(effects = standardized[!outside]), aes(x = effects)) +
        geom_histogram(aes(y = after_stat(density)), binwidth = fd.binwidth,
                       fill = "lightblue", color = "black") +
        geom_density(color = "red", linewidth = 1) +
        labs(title = panel.label, x = x.label, y = "Density",
             caption = caption) +
        theme_minimal() +
        theme(plot.title = element_text(size = 10, hjust = 0.5),
              plot.caption = element_text(size = 8, hjust = 0.5)))
}


plot.distribution <- function(effects, target, save=FALSE, filename=NULL,
                              trim = 0.005) {
    # `effects` is the two column data.frame from collect.pc.effects(standardize = "both")
    if (!is.data.frame(effects)) {
        stop("effects must be the data.frame returned by ",
             "collect.pc.effects(..., standardize = 'both')")
    }

    # The two scales differ by orders of magnitude and now want different geoms, so they
    # are built as separate plots and stacked rather than faceted. The histogram gets
    # the taller share: the boxplot is one row of ink plus its labels.
    p <- plot.raw.boxplot(effects$raw, "raw slope, b") /
        plot.std.histogram(effects$standardized,
                           "correlation, r(PC, Y) = b * sd(PC) / sd(Y)", trim,
                           x.label = "Correlation") +
        plot_layout(heights = c(1, 2)) +
        plot_annotation(
            title = paste("Distribution of PC Effects on", target),
            subtitle = paste0(nrow(effects), " effects over ",
                              length(data.with.pcs), " trios")
            # deliberately no `theme =` here: patchwork 1.2.0 warns "annotation$theme
            # is not a valid theme" under ggplot2 4.x for anything passed to it,
            # theme() included, and then ignores it. Each panel carries its own
            # theme_minimal(), and the default title styling is what we want anyway.
        )

    # Save the plot if requested
    if (save) {
        if (is.null(filename)) {
            # If no filename is provided, create a default one
            path = getwd()
            filename = paste0(path, "/PC_effects_distribution_", target, ".png")
        }
        # Save the plot to the specified filename. Taller than the single panel version
        # so the stacked facets each get room.
        ggsave(filename, plot = p, width = 7, height = 8, dpi = 150)
    }
    return(p)
}


per.trio.r2 <- function(data.with.pcs) {
    # One row per trio: how many PCs it retained, and the proportion of each gene's
    # variance they explain. simulation/verify_simulation.R plots the simulated
    # R2-versus-confounder-count relationship against this.
    #
    # R2 is the sum of squared marginal correlations rather than an lm() fit. The two are
    # identical here because principal components are orthogonal (the same identity the
    # standardized effects rely on above), and it avoids fitting 3,248 models with up to
    # 51 predictors each.
    rows <- lapply(data.with.pcs, function(x) {
        covs <- x[, -c(1:3), drop = FALSE]
        covs <- covs[, grepl("^PC[0-9]+$", colnames(covs)), drop = FALSE]
        cov.sd <- apply(covs, 2, sd)
        covs <- covs[, is.finite(cov.sd) & cov.sd > SD.TOL, drop = FALSE]
        if (ncol(covs) == 0) {
            return(NULL)
        }
        P <- as.matrix(covs)
        data.frame(nPC = ncol(P),
                   R2.cis = sum(cor(P, x[, 2])^2),
                   R2.trans = sum(cor(P, x[, 3])^2))
    })
    return(do.call(rbind, rows))
}


cache.effect.pools <- function(cis.effects, trans.effects, path, trio.r2 = NULL) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    # Write the pooled STANDARDIZED effects to disk for simulation_utils.R to resample
    # from. Only the standardized column crosses over: it is the Pearson correlation
    # cor(PC, gene), which is unit-free, so a value of 0.117 means the same thing in a
    # GTEx trio and in a simulated one. The raw slopes cannot cross over -- sd(cis gene)
    # spans 0.042 to 3.75e+04 across trios, so a raw slope carries the units of one
    # particular trio's residualized expression and is meaningless outside it.
    #
    # The pools are conditional on selection: these PCs were retained at FDR 0.05 by
    # get.conf.trios(), which is why the distributions dip at zero. That is the relevant
    # distribution for choosing simulation parameters (see the README), but it means the
    # simulated U confounders represent confounders strong enough to have been found.
    pools <- list(cis = cis.effects$standardized,
                  trans = trans.effects$standardized,
                  # one row per trio: nPC, R2.cis, R2.trans. The reference curve for
                  # verify_simulation.R's R2-versus-confounder-count figure.
                  trio.r2 = trio.r2,
                  n.trios = length(data.with.pcs),
                  source = "GTEx/data/data.with.PCs.WholeBlood.RData",
                  scale = "Pearson correlation cor(PC, gene), = b * sd(PC) / sd(Y)")

    # a raw slope leaking into the pool would break the calibration silently, and the
    # cheapest tell is a value outside [-1, 1]
    for (nm in c("cis", "trans")) {
        if (any(abs(pools[[nm]]) > 1)) {
            stop("standardized ", nm, " effects outside [-1, 1]: the pool is not on the ",
                 "correlation scale")
        }
    }

    save(pools, file = path)
    cat("\ncached effect pools ->", path, "\n")
    for (nm in c("cis", "trans")) {
        cat("  ", nm, ": n = ", length(pools[[nm]]),
            ", sd = ", round(sd(pools[[nm]]), 4),
            ", E[r^2] = ", round(mean(pools[[nm]]^2), 5),
            ", max|r| = ", round(max(abs(pools[[nm]])), 4), "\n", sep = "")
    }
    if (!is.null(trio.r2)) {
        cat("  per-trio R2: ", nrow(trio.r2), " trios, median nPC = ",
            median(trio.r2$nPC), ", median R2 cis = ", round(median(trio.r2$R2.cis), 3),
            ", trans = ", round(median(trio.r2$R2.trans), 3), "\n", sep = "")
    }
    return(invisible(pools))
}


# run for each target variable: cis, trans
# standardize = "both" so the summary and the figure cover the raw and the standardized
# scale from one pass of the regressions

# cis gene distribution of confounder effects:
print("Collecting PC effects for cis genes...")
cis.effects <- collect.pc.effects(data.with.pcs, target = "cis", standardize = "both")
summarize.effects(cis.effects, "cis")
plot.distribution(cis.effects, target = "cis", save = TRUE, filename = "PC_effects_distribution_cis.png")

print("Collecting PC effects for trans genes...")
trans.effects <- collect.pc.effects(data.with.pcs, target = "trans", standardize = "both")
summarize.effects(trans.effects, "trans")
plot.distribution(trans.effects, target = "trans", save = TRUE, filename = "PC_effects_distribution_trans.png")

print("Collecting per-trio confounder counts and R2...")
trio.r2 <- per.trio.r2(data.with.pcs)

cache.effect.pools(cis.effects, trans.effects,
                   path = file.path(root, "bioinfo_revision", "pc_distribution_invest",
                                    "data", "real_pc_effect_pools.RData"),
                   trio.r2 = trio.r2)


setwd(root)
