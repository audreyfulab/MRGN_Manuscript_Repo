# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

## this scipt runs the MRPC algorithm on the simulated data sets with all types of confounding variables
## created by sim_data_all_mods_all_confs_from_graph.R and saves the inference

source("adapted_GMAC_func/GMACpostproc.R")
library(MRGN)
library(qvalue)

params = loadRData(file = "Simulation/data/mrpc_v_mrgn_v_gmac_5k_params_all_mods_conf_types.RData")
small.datasets.all.fields=loadRData(file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types.RData")
small.datasets = lapply(small.datasets.all.fields, function(x) x$data)
just.trios=lapply(small.datasets, function(x) x[,1:3])
# tissue.name="WholeBlood"
# pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
#                           "_AllPC/PCs.matrix.",tissue.name,".RData", sep = ""))
conf.mat = loadRData(file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types_conf_mat.RData")

#confs=get.conf(trios=just.trios, cov.pool=conf.mat, blocksize = 2000, method = "correlation", selection_fdr = 0.05,
#               filter_fdr = 0.1)
confs = loadRData(file = "Simulation/data/confs_mrgn_mrpc_correlation.RData")
conf.list=lapply(confs$sig.asso.pcs, function(x,y){y[,x]}, y=conf.mat)
save(conf.list, file = "Simulation/data/mrpc_5k_conf_list_all_confs_all_mods.RData")
print(str(conf.list[[1]]))

trios.with.pcs=mapply(cbind.data.frame, just.trios, conf.list)


kc = t(as.matrix(output.WB$input.list$known.conf))

trios.with.pcs2=lapply(trios.with.pcs, function(x,y) cbind.data.frame(x,y), y = kc)
save(trios.with.pcs2, file = "Simulation/data/mrpc_data_with_SELECTED_CONFS.RData")

print(lapply(trios.with.pcs2[1:5], head))
#MRPC
print("Running MRPC...")
#----begin analysis-----#
#truth for M0
#V1-->T1
Truth.M0 <- MRPCtruth$M0
Adj.M0<- as(Truth.M0,"matrix")
#V1-->T2
Adj.M01 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M01) <- colnames(Adj.M01) <- colnames(Adj.M0)
Adj.M01[1,3] <- 1

#truth for M1
#V1-->T1-->T2
Truth.M1 <- MRPCtruth$M1
Adj.M1<- as(Truth.M1,"matrix")
#V1-->T2-->T1
Adj.M11 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M11) <- colnames(Adj.M11) <- colnames(Adj.M1)
Adj.M11[1,3] <- 1
Adj.M11[3,2] <- 1

#truth for M2
#V1-->T1<--T2
Truth.M2 <- MRPCtruth$M2
Adj.M2<- as(Truth.M2,"matrix")
#V1-->T2<--T1
Adj.M21 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M21) <-colnames(Adj.M21) <-colnames(Adj.M2)
Adj.M21[1,3] <- 1
Adj.M21[2,3] <- 1
#truth for M3
#V1-->T1, V1-->T2
Truth.M3 <- MRPCtruth$M3
Adj.M3 <- as(Truth.M3,"matrix")

#truth for M4
#V1-->T1, V1-->T2, T1--T2
Truth.M4 <- MRPCtruth$M4
Adj.M4 <- as(Truth.M4,"matrix")

############################################3
#a function to apply mrpc in parallel
apply.mrpc=function(X){


  n <- nrow (X)
  V <- colnames(X)     # Column names
  print(n)
  # Classical correlation
  suffStat <- list(C = cor(X, use = "complete.obs"),
                   n = n)

  MRPC.fit.FDR <- MRPC(X,
                       suffStat,
                       GV = 1,
                       FDR = 0.05,
                       alpha=0.01,
                       indepTest = 'gaussCItest',
                       labels = V,
                       FDRcontrol = "ADDIS",
                       verbose = FALSE)


  #plot(MRPC.fit_FDR)
  Adj.infe1 <- as( MRPC.fit.FDR@graph,"matrix")
  print(Adj.infe1)
  Adj.infe <- Adj.infe1[1:3,1:3] #only consider snp, cis, trans
  colnames(Adj.infe) <- rownames(Adj.infe) <- colnames(Adj.M01)

  #use adj to determine model
  model="Other"
  if(identical(Adj.M0,Adj.infe)){model="M0.1"}
  if(identical(Adj.M01,Adj.infe)){model="M0.2"}
  if(identical(Adj.M1,Adj.infe)){model="M1.1"}
  if(identical(Adj.M11,Adj.infe)){model="M1.2"}
  if(identical(Adj.M2,Adj.infe)){model="M2.1"}
  if(identical(Adj.M21,Adj.infe)){model="M2.2"}
  if(identical(Adj.M3,Adj.infe)){model="M3"}
  if(identical(Adj.M4,Adj.infe)){model="M4"}
  #handle others

  return(list(model=model, Adj=Adj.infe))

}

# #obtain a subset of the confounders used to generate the model
# #represents the case we dont have the exact correct model
# #compare MRPC and MRGN
# smsmall.datasets=lapply(small.datasets, function(x){
#   if((dim(x)[2]-3)>15){
#     y=x[,c(1:3, sample(4:dim(x)[2], round(runif(1, 1, 15))))]
#     return(y)
#   }else{
#     return(x)
#   }
# })

#save(smsmall.datasets, file = "Simulation/data/mrpc.small.datasets.1k.RData")


#run the MRPC inference
#mrpc.infer=lapply(small.datasets, apply.mrpc)

#get estimate of time to infer each trio
mrpc.infer.list=list()
times = NULL
for(i in 1:length(trios.with.pcs2)){
  start.time=Sys.time()
  mrpc.infer.list[[i]]=apply.mrpc(trios.with.pcs2[[i]])
  end.time=Sys.time()
  times[i]=difftime(end.time, start.time, units = 'mins')
}
params$Time.to.compute.mrpc = times

print("Done!...Saving results...")

#save
save(mrpc.infer.list, file = "Simulation/data/mrpc_5k_small_results_all_confs_all_mods.RData")
#save(params, file = "Simulation/data/mrpc_v_mrgn_v_gmac_10k_params_all_confs_all_mods.RData")
save(params, file = "Simulation/data/diagnostics/params_with_mrpc_comp_times.RData")
