# Run MRGN over the three confounder-structure simulations.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/run_structure_sims.R
#
# Companion to simulation/confounder_structure_simulation.R, which generates the data. That
# script must have been run first.
#
# THREE CASES, MRGN ONLY, n = 670. The fourth structure -- confounders + intermediate +
# common child -- is the main simulation and is not repeated here; its results are already
# in data/mrgn_group_n670.RData.
#
#   u_only   confounders only              filter_int_child = FALSE
#   u_w      confounders + 1 intermediate  filter_int_child = TRUE
#   u_z      confounders + 1 common child  filter_int_child = TRUE
#
# WHY filter_int_child IS OFF FOR u_only AND ON FOR THE OTHER TWO. The filter's job is to
# remove covariates that are intermediates or common children of the trio rather than
# confounders of it. In u_only there are none to remove -- no trio in the group contributes
# a W or Z column, so the pool contains nothing the filter could act on -- and
# get.conf.trios() does not merely no-op in that situation, it stops with "No common child
# or intermediate variables detected". select.confounders() already catches that and falls
# back to FALSE; setting it explicitly here makes the intent visible rather than relying on
# an error handler to reach the right answer. In u_w and u_z the filter is the thing under
# test: those are exactly the covariates it exists to catch.
#
# EACH CASE GETS ITS OWN out.dir. The selection cache is named selection_group_n<size>.RData
# with no other key, so two cases sharing a directory would fight over one cache file --
# selection.cache.mismatch() would catch it on cov.names and silently recompute every time,
# turning the cache into a per-run tax instead of a saving.
#
# Cases run SEQUENTIALLY, each with the full core budget, rather than three at once. MRGN's
# bootstrap already parallelises across the machine, so three concurrent processes would
# oversubscribe it threefold for no gain -- the same reasoning as run_all_inference.R's
# core split, applied across cases instead of across methods.

library(MRGN)
source("bioinfo_revision/simulation_results/inference_config.R")

SIM.DIR  <- "bioinfo_revision/simulation/simulated_data"
DATA.DIR <- "bioinfo_revision/simulation_results/data_structures"

CASES <- list(
    u_only = list(sim = "simulated_trios_u_only.RData", filter.int.child = 0,
                  label = "confounders only"),
    u_w    = list(sim = "simulated_trios_u_w.RData",    filter.int.child = 1,
                  label = "confounders + 1 intermediate"),
    u_z    = list(sim = "simulated_trios_u_z.RData",    filter.int.child = 1,
                  label = "confounders + 1 common child"))

cores <- max(1L, parallel::detectCores() - 2L)
dir.create(log.dir, recursive = TRUE, showWarnings = FALSE)

missing <- Filter(function(nm) !file.exists(file.path(SIM.DIR, CASES[[nm]]$sim)),
                  names(CASES))
if (length(missing) > 0) {
    stop("no simulated data for: ", paste(missing, collapse = ", "),
         "\nRun bioinfo_revision/simulation/confounder_structure_simulation.R first.")
}

started.all <- Sys.time()
for (nm in names(CASES)) {
    case <- CASES[[nm]]
    out <- file.path(DATA.DIR, nm)
    logfile <- file.path(log.dir, paste0("structure_", nm, ".log"))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)

    cat("\n=== ", nm, ": ", case$label, " ===\n", sep = "")
    cat("  sim   : ", file.path(SIM.DIR, case$sim), "\n", sep = "")
    cat("  out   : ", out, "\n", sep = "")
    cat("  filter_int_child: ", case$filter.int.child == 1, " | cores: ", cores, "\n", sep = "")
    cat("  log   : ", logfile, "  (tail it to watch)\n", sep = "")

    args <- c("bioinfo_revision/simulation_results/apply_mrgn.R",
              "--cores", cores,
              "--sizes", "670",
              "--sim-file", file.path(SIM.DIR, case$sim),
              "--out-dir", out,
              "--filter-int-child", case$filter.int.child)

    started <- Sys.time()
    # stdout is redirected through cmd on Windows for the same reason run_all_inference.R
    # does it: system2(stdout = ) opens the log with exclusive access there, so a running
    # job's log cannot be read until it exits.
    status <- if (.Platform$OS.type == "windows") {
        system2("cmd", c("/c", "Rscript", shQuote(args[1]), args[-1],
                         ">", shQuote(logfile), "2>&1"), wait = TRUE)
    } else {
        system2("Rscript", shQuote(args), stdout = logfile, stderr = logfile, wait = TRUE)
    }
    mins <- round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1)

    if (!identical(status, 0L)) {
        warning("case ", nm, " exited with status ", status, "; see ", logfile)
        cat("  FAILED after ", mins, " min -- see the log\n", sep = "")
    } else {
        cat("  done in ", mins, " min -> ",
            file.path(out, "mrgn_group_n670.RData"), "\n", sep = "")
    }
}

cat("\nall cases finished in",
    round(as.numeric(difftime(Sys.time(), started.all, units = "mins")), 1), "min\n")
cat("Tables: bioinfo_revision/simulation_results/results_scripts/confusion_structures.R\n")
