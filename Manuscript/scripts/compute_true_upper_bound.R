# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

path = "Manuscript/other/"
library(MRGN)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(MRPC)
library(R.utils)

#load helpers
source('Manuscript/scripts/MRGN_write_up_helper_functions.R')

source("adapted_GMAC_func/gmac_one_trio.R")

#load simulated data from files
params = loadRData(file = "Simulation/data/mrpc_v_mrgn_v_gmac_5k_params_all_mods_conf_types.RData")
small.datasets.all.fields=loadRData(file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types.RData")
small.datasets = lapply(small.datasets.all.fields, function(x) x$data[,c(1:(length(x$data)-2))])
many.conf.data = loadRData(file =
                             "Simulation/data/many_conf_data/mrgn_v_gmac_v_mrpc_300_datasets_all_mods_conf_types.RData")

many.conf.params=loadRData(file="Simulation/data/many_conf_data/mrpc_v_mrgn_v_gmac_300_params_all_mods_conf_types.RData")

# NOTE: This file originates from the MRPC_support repository. Place it at the path below.
output.WB=loadRData(fileName='GTEx/data/all_trios_output_cis.Rdata')
kc = t(as.matrix(output.WB$input.list$known.conf))
many.conf.data2 = lapply(many.conf.data, function(x) x$data[,c(1:(length(x$data)-2))])

otherpath = "Manuscript/other/"

#all params and trios
all.trios = c(small.datasets, many.conf.data2)
trios.with.pcs=lapply(all.trios, function(x,y) cbind.data.frame(x,y), y = kc)
all.params = rbind.data.frame(params, many.conf.params)

#parse ground truth 
true.mod = convert.truth(all.params$model)
true.adj = lapply(convert.truth(all.params$model), get.adj.from.class)
true.ind.score = unlist(lapply(true.adj, ind.med.edge))

#apply mrgn trio to simulated results using ground truth confounders
gt.results = lapply(trios.with.pcs, infer.trio)
save(gt.results, file=paste0(otherpath, 'TUB-mrgn-all-trios.RData'))
gt.results = loadRData(paste0(otherpath, 'TUB-mrgn-all-trios.RData'))
# parse mrgn results 
mrgn.inf = unlist(lapply(gt.results, function(x) x[14]))
mrgn.adj <- lapply(mrgn.inf, get.adj.from.class)
mrgn.ind <- unlist(lapply(mrgn.adj, ind.med.edge))


#compute upper bound of class based perf
x1.1 <- table(MRGN = convert.cats(mrgn.inf), TRUTH = true.mod)
x1.2 <- rbind(x1.1, Total = colSums(x1.1), Recall = round(diag(x1.1)/colSums(x1.1), 4))
x1.3 <- cbind(x1.2, Total = c(rowSums(x1.2[1:7,]), NA), 
              Precision = c(round(diag(x1.2)/rowSums(x1.2[1:5,]), 4), rep(NA, 3)))

write.csv(x1.3, file = paste0(otherpath,'TUB-class-based.csv'))

n = length(all.params$model)
edge.metrics.mrgn <- as.data.frame(matrix(0, nrow = n, ncol = 4))
colnames(edge.metrics.mrgn) <- c("prec_edge", "recall", "EW_prec_edge", "EW_recall")
#compute upper bound of all edge based perf
for(i in 1:n){
  truth_i <- convert.truth(all.params$model)[i]
  
  edge.metrics.mrgn[i, ] <- c(
    get.metrics(Truth = truth_i, Inferred = mrgn.adj[[i]], get.adj.truth = TRUE),
    get.metrics(Truth = truth_i, Inferred = mrgn.adj[[i]], get.adj.truth = TRUE,
                weight.edge.directed.present = 1, weight.edge.present.only = 1)
  )
}    


# Summarize by model class
mod.class <- unique(convert.truth(all.params$model))
mod.class <- sort(mod.class)

summarize_by_class <- function(metrics, model.labels) {
  tab <- sapply(mod.class, function(m) {
    colMeans(metrics[which(convert.truth(model.labels) == m), , drop = FALSE], na.rm = TRUE)
  })
  rownames(tab) <- c("Precision", "Recall", "EW Precision", "EW Recall")
  return(tab)
}

edge.met.tab.mrgn <- summarize_by_class(edge.metrics.mrgn, all.params$model)
write.csv(edge.met.tab.mrgn, file = paste0(otherpath,'TUB-edge-based.csv'))


#compute upper bound of t1-t2 edge based perf

# Get MRGN and true adjacency
mrgn.adj <- lapply(mrgn.inf, get.adj.from.class)
mrgn.adj.alpha01 <- lapply(mrgn.inf.alpha01, get.adj.from.class)
true.adj <- lapply(convert.truth(gt.combined), get.adj.from.class)

# Extract indicators
true.score <- unlist(lapply(true.adj, ind.med.edge))
mrgn.edge.ind <- unlist(lapply(mrgn.adj, ind.med.edge))
mrgn.edge.ind.alpha01 <- unlist(lapply(mrgn.adj.alpha01, ind.med.edge))
mrpc.edge.ind <- unlist(lapply(mrpc.adj, ind.med.edge))
#gmac.edge.ind <- apply(cbind(gmac.cis$output.table$Cis_Sig, gmac.trans$output.table$Trans_Sig), 1, ind.gmac)
gmac.edge.ind.at05.sc <- apply(gmac.05, 1, ind.gmac)
gmac.edge.ind.at01.sc <- apply(gmac.01, 1, ind.gmac)

# Helper to build confusion matrix with PR
build_table <- function(pred, truth, add_zero = FALSE) {
  t <- if (add_zero) rbind(c(0, 0), table(pred, truth)) else table(pred, truth)
  colnames(t) <- c("T1-T2 Absent", "T1-T2 Present")
  rownames(t) <- c("T1-T2 Pred. Absent", "T1-T2 Pred. Present")
  t2 <- rbind(t, Total = colSums(t), Recall = round(diag(t)/colSums(t), 4))
  t3 <- cbind(t2, Total = c(rowSums(t2)[1:3], NA), Precision = c(round(diag(t)/rowSums(t), 4), rep("", 2)))
  return(list(raw = t, formatted = t3))
}

# Tables
t1 <- build_table(mrgn.ind, true.ind.score)


write.csv(t1$formatted, file = paste0(path, 'TUB-t1-t2-edge.csv'))



































































