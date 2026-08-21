# Run MRGN, MRPC and GMAC concurrently, then combine.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/run_all_inference.R
#
# Starts apply_mrgn.R, apply_mrpc.R and apply_gmac.R as separate OS processes, waits for
# all three, and writes the joined master. Each apply script also runs perfectly well on
# its own -- this only removes the waiting.
#
# WHY SEPARATE PROCESSES rather than a three-worker cluster: MRGN's bootstrap and GMAC
# both parallelise internally through their own cluster. Putting the methods themselves on
# cluster workers would force those internals to run serially, costing far more than
# method-level concurrency gains.
#
# TWO THINGS THIS DOES that launching three terminals by hand does not:
#
#   1. Confounder selection runs FIRST, to completion. MRGN and MRPC both need it, and
#      starting them together against a cold cache would have each spend ~30 min per group
#      computing the identical thing before racing to write the same file. Done once here,
#      the method processes only ever read it.
#   2. The core budget is divided rather than tripled. Each process otherwise sizes its
#      cluster from detectCores() and the three together oversubscribe the machine.

library(MRGN)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

METHODS      <- ALL.METHODS   # subset to skip a method entirely
POLL.SECONDS <- 30

apply.script <- function(method) {
    file.path("bioinfo_revision", "simulation_results", paste0("apply_", method, ".R"))
}

# Start one method in the background, logging to `logfile`.
#
# On Windows this deliberately does NOT use system2(stdout = logfile). That opens the log
# with exclusive access, so nothing else can read it until the process exits -- the file
# is visible but every read fails with "Permission denied" / "Device or resource busy",
# including from the editor. Redirecting through cmd instead opens it shared, so a running
# job's log can be tailed live:
#
#   PowerShell:  Get-Content bioinfo_revision/simulation_results/logs/apply_mrgn.log -Wait
#
# system2(stdout = ) is fine on Unix, which does not lock like this.
launch.background <- function(script, args, logfile) {
    if (.Platform$OS.type == "windows") {
        system2("cmd", c("/c", "Rscript", shQuote(script), args,
                         ">", shQuote(logfile), "2>&1"), wait = FALSE)
    } else {
        system2("Rscript", c(shQuote(script), args),
                stdout = logfile, stderr = logfile, wait = FALSE)
    }
}

# Core split. MRPC never builds a cluster -- MRPC() is single threaded -- so it gets one
# core and the rest goes to the two methods that can use it. MRGN takes the larger share:
# its bootstrap is n.bootstrap replicates x 3 confounder settings x every trio.
core.budget <- function(methods) {
    total <- max(1L, parallel::detectCores() - 2L)
    share <- setNames(rep(1L, length(methods)), methods)
    rest <- setdiff(methods, "mrpc")
    if (length(rest) == 0) return(share)
    avail <- if ("mrpc" %in% methods) max(1L, total - 1L) else total
    w <- ifelse(rest == "mrgn", 0.6, 0.4)
    share[rest] <- pmax(1L, as.integer(round(avail * w / sum(w))))
    share
}


# ---------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------

sim_data <- loadRData(file = sim.data.file)
group.of <- sapply(sim_data, function(x) x$params$sample.size)
sizes <- if (is.null(sample.sizes)) sort(unique(group.of)) else sample.sizes
cat("loaded", length(sim_data), "datasets |", length(sizes), "sample-size groups\n")

# ---- 1. selection first, once ----
cat("\n=== confounder selection (shared by MRGN and MRPC) ===\n")
if (any(c("mrgn", "mrpc") %in% METHODS)) {
    for (size in sizes) {
        idx <- which(group.of == size)
        # must match what the method processes will use, or their cache lookup misses
        if (!is.null(max.per.group)) idx <- head(idx, max.per.group)
        cat("  n =", size, "|", length(idx), "trios\n")
        invisible(get.selection(sim_data[idx], out.dir = out.dir, sim.file = sim.data.file,
                                rerun = rerun.selection, save = is.null(max.per.group),
                                verbose = TRUE,
                                selection_fdr = selection_fdr, filter_fdr = filter_fdr,
                                alpha = alpha, filter_int_child = filter_int_child))
    }
} else {
    cat("  neither MRGN nor MRPC requested; skipping\n")
}
rm(sim_data); invisible(gc())

# ---- 2. what still needs running ----
# The checkpoints are both the work list and the completion signal: a (method, group) with
# a file already there is not launched and not waited for, and one written after launch is
# done. Nothing else has to be tracked, and a stale checkpoint cannot be mistaken for work
# this run did.
work <- expand.grid(method = METHODS, size = sizes, stringsAsFactors = FALSE)
work$file <- mapply(method.checkpoint, out.dir, work$method, work$size)
work$pending <- rerun.inference | !file.exists(work$file)
todo <- unique(work$method[work$pending])

if (length(todo) == 0) {
    cat("\nevery requested method already has checkpoints for every group.",
        "Set rerun.inference <- TRUE in inference_config.R to redo them.\n")
} else {
    cores <- core.budget(todo)
    dir.create(log.dir, showWarnings = FALSE, recursive = TRUE)
    started <- Sys.time()

    cat("\n=== launching", length(todo), "processes ===\n")
    for (m in todo) {
        logfile <- file.path(log.dir, paste0("apply_", m, ".log"))
        args <- c("--cores", cores[[m]])
        if (!is.null(sample.sizes)) args <- c(args, "--sizes", paste(sizes, collapse = ","))
        if (!is.null(max.per.group)) args <- c(args, "--max-per-group", max.per.group)
        launch.background(apply.script(m), args, logfile)
        cat(sprintf("  %-5s %2d cores -> %s\n", m, cores[[m]], logfile))
    }

    # ---- 3. wait ----
    want <- work[work$pending & work$method %in% todo, ]
    cat("\nwaiting for", nrow(want), "group checkpoints. Logs are live in", log.dir, "\n")
    per.method <- table(want$method)
    repeat {
        done <- file.exists(want$file)
        done[done] <- file.info(want$file[done])$mtime >= started
        done[is.na(done)] <- FALSE
        by.m <- tapply(done, want$method, sum)
        cat(sprintf("\r  %s  (%d/%d) ",
                    paste(sprintf("%s %d/%d", names(by.m), by.m,
                                  per.method[names(by.m)]), collapse = "   "),
                    sum(done), nrow(want)))
        flush.console()
        if (all(done)) break
        Sys.sleep(POLL.SECONDS)
    }
    cat("\n\nall processes finished in",
        round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), "min\n")
}

# ---- 4. combine once ----
# Reads every method checkpoint on disk, not just the ones launched here, so a method run
# separately earlier still makes it into the master.
cat("\n=== combining ===\n")
inference.results <- combine.all(sample.sizes = sample.sizes)
