# create_combined_T1T2_edge_table.R
# ============================================================================
# Creates the combined T1-T2 edge results table (1800 trios) from both the
# <=15 confounder (1500 trios) and 15-50 confounder (300 trios) scenarios.
# Includes Type I and Type II error rates for each method.
#
# Saves to: Manuscript/supplementary_tables/T1.T2.Edge.Results.Combined.csv
#
# NOTE: Set working directory to repository root before running.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
# ============================================================================

library("MRGN")
source('Manuscript/scripts/MRGN_write_up_helper_functions.R')

path_supptabs = "Manuscript/supplementary_tables/"

# ============================================================================
# Load data: <=15 confounder scenario (1500 trios)
# ============================================================================
params = loadRData(file = "Simulation/data/mrpc_v_mrgn_v_gmac_5k_params_all_mods_conf_types.RData")
gmac.cis = loadRData(file = "Simulation/data/gmac_5k_cis_results_all_mods_all_conf_types_proproc.RData")
gmac.trans = loadRData(file = "Simulation/data/gmac_5k_trans_results_all_all_mods_conf_types_preproc.RData")
mrgn.inf = loadRData("Simulation/data/int_and_child_filtered_data/mrgn_5k_inf_results_all_confs_filtered_all_mods.RData")
mrgn.inf.alpha01 = loadRData('Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_inference_results.RData')

# ============================================================================
# Load data: 15-50 confounder scenario (300 trios)
# ============================================================================
many.conf.params = loadRData(file = "Simulation/data/many_conf_data/mrpc_v_mrgn_v_gmac_300_params_all_mods_conf_types.RData")
gmac.cis.many.conf = loadRData("Simulation/data/many_conf_data/gmac_300_cis_results_all_mods_all_conf_types_proproc.RData")
gmac.trans.many.conf = loadRData("Simulation/data/many_conf_data/gmac_300_trans_results_all_all_mods_conf_types_preproc.RData")
mrgn.many.conf.inf = loadRData(file = "Simulation/data/many_conf_data/mrgn_300_inf_results_all_confs_all_mods.RData")
mrgn.mc.inf.alpha01 = loadRData(file = "Simulation/data/many_conf_data/alpha_01_selected_confs_results/mrgn_many_conf_liberal_alpha01_inf_results.RData")

# ============================================================================
# Load MRPC results (combined across both scenarios)
# ============================================================================
mrpc.inf.all = loadRData("Manuscript/other/results_all_MRPC-ADDIS.RData")

# ============================================================================
# Combine both scenarios
# ============================================================================
gt.combined = c(params$model, many.conf.params$model)

mrgn.inf.combined = c(mrgn.inf, mrgn.many.conf.inf)
mrgn.inf.alpha01.combined = c(mrgn.inf.alpha01, mrgn.mc.inf.alpha01)

gmac.05.combined = cbind(c(gmac.cis$output.table$Cis_at_05_cutoff,
                           gmac.cis.many.conf$output.table$Cis_at_05_cutoff),
                         c(gmac.trans$output.table$Trans_at_05_cutoff,
                           gmac.trans.many.conf$output.table$Trans_at_05_cutoff))

gmac.01.combined = cbind(c(gmac.cis$output.table$Cis_at_01_cutoff,
                           gmac.cis.many.conf$output.table$Cis_at_01_cutoff),
                         c(gmac.trans$output.table$Trans_at_01_cutoff,
                           gmac.trans.many.conf$output.table$Trans_at_01_cutoff))

# ============================================================================
# Compute adjacency / edge indicators for all 1800 trios
# ============================================================================

# MRGN
mrgn.adj = lapply(mrgn.inf.combined, get.adj.from.class)
mrgn.adj.alpha01 = lapply(mrgn.inf.alpha01.combined, get.adj.from.class)

# Truth
true.adj = lapply(convert.truth(gt.combined), get.adj.from.class)
true.score = unlist(lapply(true.adj, ind.med.edge))

# MRGN edge indicators
mrgn.edge.ind = unlist(lapply(mrgn.adj, ind.med.edge))
mrgn.edge.ind.alpha01 = unlist(lapply(mrgn.adj.alpha01, ind.med.edge))

# MRPC edge indicators (excluding trios that did not finish)
mrpc.inf2 = unlist(lapply(mrpc.inf.all, function(x) ifelse(is.null(x$model), 'did not finish', x$model)))
idx.dnf = which(mrpc.inf2 == 'did not finish')
mrpc.adj = lapply(mrpc.inf.all, function(x) x$Adj)[-idx.dnf]
mrpc.edge.ind = unlist(lapply(mrpc.adj, ind.med.edge))

# GMAC edge indicators
gmac.edge.ind.at05 = apply(gmac.05.combined, 1, ind.gmac)
gmac.edge.ind.at01 = apply(gmac.01.combined, 1, ind.gmac)

# ============================================================================
# Helper: build confusion matrix with Precision, Recall, Type I & Type II Error
# ============================================================================
build_table <- function(pred, truth, add_zero = FALSE) {
  t <- if (add_zero) rbind(c(0, 0), table(pred, truth)) else table(pred, truth)
  colnames(t) <- c("T1-T2 Absent", "T1-T2 Present")
  rownames(t) <- c("T1-T2 Pred. Absent", "T1-T2 Pred. Present")
  t2 <- rbind(t, Total = colSums(t), Recall = round(diag(t)/colSums(t), 4))
  t3 <- cbind(t2,
              Total = c(rowSums(t2)[1:3], NA),
              Precision = c(round(diag(t)/rowSums(t), 4), rep("", 2)))
  # Type I Error  = FP / (FP + TN) = predicted present when truly absent / total truly absent
  # Type II Error = FN / (FN + TP) = predicted absent when truly present / total truly present
  t3 <- rbind(t3,
              `Type I Error`  = c(round(t[2,1]/colSums(t)[1], 4), "", "", ""),
              `Type II Error` = c("", round(t[1,2]/colSums(t)[2], 4), "", ""))
  return(t3)
}

# ============================================================================
# Build confusion matrices for each method
# ============================================================================
t1  <- build_table(mrgn.edge.ind, true.score)                    # MRGN + CS FDR
t12 <- build_table(mrgn.edge.ind.alpha01, true.score)            # MRGN + CS noFDR
t2  <- build_table(mrpc.edge.ind, true.score[-idx.dnf])          # MRPC-ADDIS
t31 <- build_table(gmac.edge.ind.at05, true.score)               # GMAC alpha<0.05
t32 <- build_table(gmac.edge.ind.at01, true.score)               # GMAC alpha<0.01

# ============================================================================
# Assemble combined display table
# ============================================================================
nr <- nrow(t1)  # 6 rows per block (Absent, Present, Total, Recall, Type I, Type II)

# Label appears at row 3 (Total) of each block.
# Between labels: (nr - 3) remaining rows + 1 separator + 2 pre-label rows = nr
# After last label: (nr - 3) remaining rows
final.t123 <- cbind(
  `Inference Method` = c(rep(NA, 2), "MRGN",
                         rep(NA, nr), "MRGN",
                         rep(NA, nr), "MRPC",
                         rep(NA, nr), "GMAC",
                         rep(NA, nr), "GMAC", rep(NA, nr - 3)),

  `Inference Correction` = c(rep(NA, 2), "None: alpha < 0.01",
                             rep(NA, nr), "None: alpha < 0.01",
                             rep(NA, nr), "ADDIS",
                             rep(NA, nr), "None: alpha < 0.05",
                             rep(NA, nr), "None: alpha < 0.01", rep(NA, nr - 3)),

  `Confounder Selection Correction` = c(rep(NA, 2), "FDR < 0.05",
                                        rep(NA, nr), "None: alpha < 0.01",
                                        rep(NA, nr), "FDR < 0.05",
                                        rep(NA, nr), "FDR < 0.05",
                                        rep(NA, nr), "FDR < 0.05", rep(NA, nr - 3)),

  Description = c(rownames(t1), NA,
                  rownames(t12), NA,
                  rownames(t2), NA,
                  rownames(t31), NA,
                  rownames(t32)),

  rbind(t1, rep(NA, 4),
        t12, rep(NA, 4),
        t2, rep(NA, 4),
        t31, rep(NA, 4),
        t32)
)

# ============================================================================
# Save
# ============================================================================
write.csv(final.t123, file = paste0(path_supptabs, "T1.T2.Edge.Results.Combined.csv"), row.names = FALSE)

cat("Saved:", paste0(path_supptabs, "T1.T2.Edge.Results.Combined.csv"), "\n")
cat("Total trios (MRGN):", length(mrgn.edge.ind), "\n")
cat("Total trios (MRPC):", length(mrpc.edge.ind), "(", length(idx.dnf), "did not finish)\n")
cat("Total trios (GMAC):", length(gmac.edge.ind.at05), "\n")
