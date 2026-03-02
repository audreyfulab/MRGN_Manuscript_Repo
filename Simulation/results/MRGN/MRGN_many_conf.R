# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

library(MRGN)
source("adapted_GMAC_func/GMACpostproc.R")


params = loadRData(file = "Simulation/data/many_conf_data/mrpc_v_mrgn_v_gmac_300_params_all_mods_conf_types.RData")
sim.datasets.all.fields=loadRData(file = "Simulation/data/many_conf_data/mrgn_v_gmac_v_mrpc_300_datasets_all_mods_conf_types.RData")
sim.datasets = lapply(sim.datasets.all.fields, function(x) x$data)
just.trios=lapply(sim.datasets, function(x) x[,1:3])
# tissue.name="WholeBlood"
# NOTE: The following path refers to external data not included in this repository.
# pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
#                           "_AllPC/PCs.matrix.",tissue.name,".RData", sep = ""))
conf.mat = loadRData(file = "Simulation/data/many_conf_data/mrgn_v_gmac_v_mrpc_300_datasets_all_mods_conf_types_conf_mat.RData")

confs=get.conf.trios(trios=just.trios, cov.pool=conf.mat, blocksize = 100, selection_fdr = 0.05,
                     filter_fdr = 0.1, filter_int_child = T, adjust_by = 'all')

save(confs, file = "Simulation/data/many_conf_data/confs_mrgn_mrpc_REG.RData")
#confs=loadRData(file = "Simulation/data/many_conf_data/confs_mrgn_mrpc_REG.RData")

conf.list=lapply(confs$sig.asso.covs, function(x,y){y[,x]}, y=conf.mat)

save(conf.list, file = "Simulation/data/many_conf_data/mrgn_300_conf_list_all_confs_all_mods_REG.RData")

print(str(conf.list[[1]]))
trios.with.pcs=mapply(cbind.data.frame, just.trios, conf.list)

#add in known confounders

kc = t(as.matrix(output.WB$input.list$known.conf))

trios.with.pcs2=lapply(trios.with.pcs, function(x,y) cbind.data.frame(x,y), y = kc)
save(trios.with.pcs2, file = "Simulation/data/many_conf_data/mrgn_data_300_with_SELECTED_CONFS_REG.RData")

print(lapply(trios.with.pcs2[1:5], head))

#preform regressions and classify model types
#MRGN
print("Running MRGN")
reg.res=list()
inf.mods=NULL
times=NULL
#regression
# reg.res=sapply(trios.with.pcs2, infer.trio, nperm=1000)
# #model class
# inf.mods=as.vector(unlist(reg.res[14,]))

#get estimate of time to compute each trio
for(i in 1:length(trios.with.pcs2)){
  start.time=Sys.time()
  reg.res[[i]]=infer.trio(trios.with.pcs2[[i]])
  inf.mods[i] = unlist(reg.res[[i]][14])
  end.time=Sys.time()
  times[i]=difftime(end.time, start.time, units = 'mins')
}
params$Time.to.compute.mrgn = times
print("Done!...Saving results...")

#save
save(reg.res, file = "Simulation/data/many_conf_data/mrgn_300_regres_results_all_confs_all_mods.RData")
save(inf.mods, file = "Simulation/data/many_conf_data/mrgn_300_inf_results_all_confs_all_mods.RData")
save(params, file = "Simulation/data/diagnostics/mrgn_many_conf_with_comp_times.RData")
