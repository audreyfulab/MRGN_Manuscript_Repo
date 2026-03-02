# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

## this scipt runs the GMAC algorithm on the simulated data sets with all types of confounding variables
## created by sim_data_all_mods_all_confs_from_graph.R and saves the inference
rm(list=ls())

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

source("adapted_GMAC_func/GMAC_moded/R/GMAC.R")
#cl <- makeCluster(10)
#source("adapted_GMAC_func/GMACpostproc.R")
#library(MRGN)
#library(qvalue)

sim.datasets.all.fields=loadRData(file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types.RData")
sim.datasets = lapply(sim.datasets.all.fields, function(x) x$data)
output.WB=loadRData(fileName='GTEx/data/all_trios_output_cis.Rdata')
kc = t(as.matrix(output.WB$input.list$known.conf))
#load in pc mat from whole blood
# tissue.name="WholeBlood"
# pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
#                           "_AllPC/PCs.matrix.",tissue.name,".RData", sep = ""))
conf.mat = loadRData(file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types_conf_mat.RData")

#perform GMAC inferences
#GMAC
print("Assembling data for gmac...")
num.trios=length(sim.datasets)
sim.datasets.red=lapply(sim.datasets, function(x)x[,1:3])
sim.datasets.mat=do.call("cbind", sim.datasets.red)
snp.idx=seq(1, dim(sim.datasets.mat)[2]-2, 3)
snp.dat.cis=sim.datasets.mat[,snp.idx]
exp.dat=sim.datasets.mat[,-snp.idx]
trios.idx=cbind(c(1:num.trios), matrix(c(1:dim(exp.dat)[2]), nrow = num.trios, ncol = 2, byrow = T))


dim(t(exp.dat))
dim(t(snp.dat.cis))
dim(t(conf.mat))

input.list=list(kc = t(kc), cov.pool=t(conf.mat), snp.dat.cis=t(snp.dat.cis), exp.dat=t(exp.dat), trio.indexes=trios.idx,
                nperm=1000, use.nominal.p=TRUE, fdr = 0.05, fdr_filter = 0.1)


print("applying gmac with cis mediators")
output.cis.med <- gmac(known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                       exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                       trios.idx = input.list$trio.indexes, nperm = input.list$nperm,
                       nominal.p = input.list$use.nominal.p, fdr = input.list$fdr,
                       fdr_filter = input.list$fdr_filter)

print(str(output.cis.med))

print("applying gmac with trans mediators")
output.trans.med <- gmac(known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                         exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                         trios.idx = input.list$trio.indexes[,c(1,3,2)], nperm = input.list$nperm,
                         nominal.p = input.list$use.nominal.p,fdr = input.list$fdr,
                         fdr_filter = input.list$fdr_filter)
#stopCluster(cl)
print(str(output.trans.med))

#reorganize output for saving
#cis results
out.table.cis=cbind.data.frame(output.cis.med[[1]], output.cis.med[[2]])
colnames(out.table.cis)=c(paste0('pval_', colnames(output.cis.med[[1]])),
                          paste0('effect_change_', colnames(output.cis.med[[2]])))

out.list.cis=list(out.table.cis, input.list, output.cis.med[[3]])
names(out.list.cis)=c("output.table", "input.list", "cov.indicator.list")

#trans results
out.table.trans=cbind.data.frame(output.trans.med[[1]], output.trans.med[[2]])
colnames(out.table.cis)=c(paste0('pval_', colnames(output.trans.med[[1]])),
                          paste0('effect_change_', colnames(output.trans.med[[2]])))

out.list.trans=list(out.table.trans, input.list, output.trans.med[[3]])
names(out.list.trans)=c("output.table", "input.list", "cov.indicator.list")



print("Done!...Saving...")

save(out.list.trans, file = "Simulation/data/gmac_5k_trans_results_all_confs_all_mods.RData")
save(out.list.cis, file = "Simulation/data/gmac_5k_cis_results_all_confs_all_mods.RData")
