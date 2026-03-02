# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

## this scipt runs the MRGN algorithm on the simulated data sets with all types of confounding variables
## created by sim_data_all_mods_all_confs_from_graph.R and saves the inference

library(MRGN)
source("adapted_GMAC_func/GMACpostproc.R")


params = loadRData(file = "Simulation/data/mrpc_v_mrgn_v_gmac_5k_params_all_mods_conf_types.RData")
sim.datasets.all.fields=loadRData(file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types.RData")
sim.datasets = lapply(sim.datasets.all.fields, function(x) x$data)
just.trios=lapply(sim.datasets, function(x) x[,1:3])
# tissue.name="WholeBlood"
# pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
#                           "_AllPC/PCs.matrix.",tissue.name,".RData", sep = ""))
conf.mat = loadRData(file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types_conf_mat.RData")

confs=get.conf.trios(trios=just.trios, cov.pool=conf.mat, blocksize = 500, filter_int_child = TRUE,
                     selection_fdr = 0.05, filter_fdr = 0.1, adjust_by = 'none', alpha = 0.01)
#confs2=get.conf(trios=just.trios, PCscores=conf.mat, blocksize = 2000, method = "regression")
save(confs, file = "Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_conf_output.RData")
print(confs$sig.asso.covs[[1]])
#save(confs, file = "Simulation/data/int_and_child_filtered_data/confs_mrgn_mrpc_regression.RData")
conf.list=lapply(confs$sig.asso.covs, function(x,y){y[,x]}, y=conf.mat)
#conf.list2 = lapply(confs$sig.asso.pcs, function(x,y){y[,x[[1]]]}, y=conf.mat)
save(conf.list, file = "Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_conf_list.RData")

print(str(conf.list[[1]]))
trios.with.pcs=mapply(cbind.data.frame, just.trios, conf.list)

#add in known confounders

kc = t(as.matrix(output.WB$input.list$known.conf))

trios.with.pcs2=lapply(trios.with.pcs, function(x,y) cbind.data.frame(x,y), y = kc)
save(trios.with.pcs2, file = "Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_DATA_W_SELECTED_CONFS.RData")

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
save(params, file = "Simulation/data/diagnostics/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_params_w_comptimes.RData")
save(reg.res, file = "Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_regression_results.RData")
save(inf.mods, file = "Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_inference_results.RData")
