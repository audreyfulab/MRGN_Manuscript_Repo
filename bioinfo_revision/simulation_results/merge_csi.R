# Fold the CS-i columns into the per-method master results.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/merge_csi.R
#
# Reads  data/inference_<method>_csi.RData   (written by apply_<method>_csi.R)
# Writes data/inference_<method>.RData/.csv  (the masters the results stage reads)
#
# The CS-i arm was run as its own pass so that the arms already on disk would not be
# recomputed -- see the header of apply_mrgn_csi.R. That leaves the CS-i results in a
# separate file, and the results stage (results_scripts/) reads one master per method. This
# joins them.
#
# ---------------------------------------------------------------------------------------
# WHAT IS COPIED, AND WHAT IS CHECKED
# ---------------------------------------------------------------------------------------
#
# ONLY the CS-i columns move: `CSi.*` (the selection scores) and `<method>.CSi.*` (the
# fit). Everything else in the master is left exactly as it was. The two files also share
# every id column and, for MRGN/MRPC/MR-GGI, the CS-q and CS-alpha selection scores --
# those are NOT copied, they are COMPARED, because they were computed from the same
# selection cache in both passes and any disagreement would mean the two runs saw different
# confounder sets. The script stops rather than writing if they differ.
#
# Rows are matched on `dataset`, not on position. The masters are sorted by dataset and so
# are the CS-i files, but a method whose CS-i pass covered fewer groups -- MRPC, which has
# no n = 670 or n = 1000 -- would silently shift every row if position were trusted.
#
# A master that already has CS-i columns is overwritten with the new values rather than
# gaining duplicates, so re-running this after re-running an arm is safe.
#
# Idempotent and non-destructive to the CS-i files: they stay on disk as provenance.

# MRGN for loadRData(), which inference_utils.R uses but does not define.
library(MRGN)
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

METHODS <- c("mrgn", "mrpc", "gmac", "mrggi")

read.results <- function(method) {
    path <- file.path(out.dir, sprintf("inference_%s.RData", method))
    if (!file.exists(path)) return(NULL)
    loadRData(path)
}

for (method in METHODS) {
    cat("\n===", method, "===\n")
    master <- read.results(method)
    csi    <- read.results(paste0(method, "_csi"))

    if (is.null(master)) { cat("  no master results -- skipping\n"); next }
    if (is.null(csi))    { cat("  no CS-i results yet -- skipping\n"); next }

    csi.cols <- grep(sprintf("^(CSi[.]|%s[.]CSi[.])", method), names(csi), value = TRUE)
    if (length(csi.cols) == 0) {
        stop(method, ": inference_", method, "_csi.RData has no CS-i columns")
    }

    # Integrity: the columns both passes ACTUALLY COMPUTED must agree. Those are the id
    # columns and the CS-q / CS-alpha selection scores, which both passes derived from the
    # same selection cache -- a difference there would mean the two runs saw different
    # confounder sets.
    #
    # The other arms' result columns are excluded, and must be: the CS-i pass ran with
    # arms = "CSi", so it wrote mrgn.truth.*, mrgn.CSq.* and mrgn.CSa.* as all-NA
    # not-attempted blocks. Comparing those against the master's real values is comparing a
    # column that was computed against one that was deliberately skipped.
    shared <- setdiff(intersect(names(master), names(csi)), csi.cols)
    shared <- shared[!grepl(paste0("^", method, "[.]"), shared)]
    key <- match(csi$dataset, master$dataset)
    if (anyNA(key)) {
        stop(method, ": ", sum(is.na(key)), " CS-i rows have a dataset id that is not in ",
             "the master. The two passes ran over different simulated data.")
    }
    for (cc in shared) {
        a <- master[[cc]][key]; b <- csi[[cc]]
        if (!isTRUE(all.equal(a, b, check.attributes = FALSE))) {
            stop(method, ": column '", cc, "' differs between the master and the CS-i ",
                 "pass. Both were built from the same selection cache, so this means one ",
                 "of them is stale -- do not merge until it is resolved.")
        }
    }
    cat("  ", length(shared), " shared columns agree | matching ", nrow(csi),
        " of ", nrow(master), " master rows on dataset\n", sep = "")

    # NA for master rows the CS-i pass did not cover (MRPC's n = 670 / n = 1000). That is
    # "not run", and it reads the same way the disabled-arm columns already do.
    for (cc in csi.cols) {
        col <- csi[[cc]]
        filled <- rep(if (is.character(col)) NA_character_ else
                      if (is.logical(col)) NA else NA_real_, nrow(master))
        filled[key] <- col
        master[[cc]] <- filled
    }

    covered <- sort(unique(csi$sample.size))
    missing <- setdiff(sort(unique(master$sample.size)), covered)
    cat("  added ", length(csi.cols), " CS-i columns | sizes covered: ",
        paste(covered, collapse = ", "),
        if (length(missing)) paste0(" | NOT covered: ", paste(missing, collapse = ", ")) else "",
        "\n", sep = "")

    assign(paste0(method, ".results"), master)
    base::save(list = paste0(method, ".results"),
               file = file.path(out.dir, sprintf("inference_%s.RData", method)))
    write.csv(master, file.path(out.dir, sprintf("inference_%s.csv", method)),
              row.names = FALSE)
    cat("  wrote inference_", method, ".RData/.csv | ", nrow(master), " rows, ",
        ncol(master), " columns\n", sep = "")
}

cat("\ndone.\n")
