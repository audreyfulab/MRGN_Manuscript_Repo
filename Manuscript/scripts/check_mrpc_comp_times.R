
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

library(MRGN)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(MRPC)
library(R.utils)

source("adapted_GMAC_func/gmac_one_trio.R")

path = "Manuscript/other/"
# NOTE: External path from Doctoral_Thesis repo - update if needed
#savepath = 'C:/Users/Bruin/OneDrive/Documents/GitHub/Doctoral_Thesis/Tables/mrgn/mrgntrio_mrpc_comp_times/'
# NOTE: External file from MRPC_support repo - update path as needed
#source("C:/Users/Bruin/OneDrive/Documents/GitHub/MRPC_support/Analysis with GMAC/GMACpostproc.R")
#library("MRGN")
#library("qvalue")

params = loadRData(file = "Simulation/data/mrpc_v_mrgn_v_gmac_5k_params_all_mods_conf_types.RData")
small.datasets.all.fields=loadRData(file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types.RData")
small.datasets = lapply(small.datasets.all.fields, function(x) x$data)
just.trios=lapply(small.datasets, function(x) x[,1:3])
# tissue.name="WholeBlood"
# NOTE: External server path - update if needed
# pc.matrix=loadRData(paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",tissue.name,
#                           "_AllPC/PCs.matrix.",tissue.name,".RData", sep = ""))
# conf.mat = loadRData(file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types_conf_mat.RData")

# confs=get.conf.trios(trios=just.trios, cov.pool=conf.mat, blocksize = 500, filter_int_child = TRUE,
#                      selection_fdr = 0.05, filter_fdr = 0.1, adjust_by = 'all')
# print(confs$sig.asso.covs[[1]])
# 
# conf.list=lapply(confs$sig.asso.covs, function(x,y){y[,x]}, y=conf.mat)
conf.list = loadRData(file = "Simulation/data/int_and_child_filtered_data/mrpc_5k_conf_list_all_confs_filtered_all_mods.RData")
print(str(conf.list[[1]]))

trios.with.pcs=mapply(cbind.data.frame, just.trios, conf.list)

# NOTE: This file originates from the MRPC_support repository. Place it at the path below.
output.WB=loadRData(fileName='GTEx/data/all_trios_output_cis.Rdata')
kc = t(as.matrix(output.WB$input.list$known.conf))

trios.with.pcs2=lapply(trios.with.pcs, function(x,y) cbind.data.frame(x,y), y = kc)

many.conf.data = loadRData(file =
                             "Simulation/data/many_conf_data/mrgn_v_gmac_v_mrpc_300_datasets_all_mods_conf_types.RData")

many.conf.data2 = lapply(many.conf.data, function(x) x$data)


trios.with.pcs3 = c(trios.with.pcs2, many.conf.data2)

#this reads in the computation times for applying MPRC to 1500 trios from simulations with 1 - 15 confounders
#MRPC was applied to each trio with the !!! Selected confounders !!! using the CSFDR method 
mrpc.params = loadRData('Simulation/data/diagnostics/mrpc_comp_times_with_filtering.RData')

mrgn.params = loadRData('Simulation/data/diagnostics/mrgn_comp_times_with_filtering.RData')



sort_ix = sort(mrpc.params$number.of.Uconfs, decreasing = F, index.return = T)$ix

dat = cbind.data.frame(`Number of Confounders` = mrpc.params$number.of.Uconfs[sort_ix],
                       `Time To Compute Trio (min)` = mrpc.params$Time.to.compute.mrpc[sort_ix],
                       Model = mrpc.params$model[sort_ix])


# Plot shows one trio (#228) took an extremely long time to compute (over 4 hours to compute)
# Additional trios also took more than an hour
P = ggplot(data = dat, aes(x = `Number of Confounders`, y = `Time To Compute Trio (min)`))+
  geom_point(size = 2.5)+
  geom_smooth()+
  theme_hc()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'top', legend.text = element_text(size = 16),axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 16),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5))+
  ggtitle("Time to compute simulated trios with MRPC-ADDIS", subtitle = 'Scenario with 1 - 15 confounders')

pdf(paste0(path, 'mrpc_time_to_compute_trios.pdf'))
plot(P)
dev.off()




bad.trio.idx = which(mrpc.params$Time.to.compute.mrpc == max(mrpc.params$Time.to.compute.mrpc))

bad.trio = trios.with.pcs2[[bad.trio.idx]]
#write bad.trio to file
write.csv(bad.trio, paste0(path, 'Bad_simulated_trio.csv'))





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
apply.mrpc=function(X, FDRmethod = 'ADDIS'){
  
  
  n <- nrow (X)
  V <- colnames(X)     # Column names
  print(n)
  # Classical correlation
  suffStat <- list(C = cor(X, use = "complete.obs"),
                   n = n)
  start.time=Sys.time()
  MRPC.fit.FDR <- MRPC(X,
                       suffStat,
                       GV = 1,
                       FDR = 0.05,
                       alpha=0.01,
                       indepTest = 'gaussCItest',
                       labels = V,
                       FDRcontrol = FDRmethod,
                       verbose = FALSE)
  end.time=Sys.time()
  ttc=difftime(end.time, start.time, units = 'mins')
  
  #plot(MRPC.fit_FDR)
  Adj.infe1 <- as( MRPC.fit.FDR@graph,"matrix")
  #print(Adj.infe1)
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
  
  return(list(model=model, Adj=Adj.infe, time = ttc))
  
}



times = NULL
times.gmac = NULL
results.all.mrpc = list()
results.gmac = list()

for(i in seq_along(trios.with.pcs3)){
  tryCatch({
    result <- withTimeout({
      apply.mrpc(trios.with.pcs3[[i]], FDRmethod = 'ADDIS')
    }, timeout = 120, onTimeout = "error")
    print(paste0('finished results for trio ', i, ': inferred model ', result$model, ' time: ', result$time))
    results.all.mrpc[[i]] <- result
    times[i] <- result$time
  }, error = function(e) {
    message(paste0("Timeout or error in trio ", i, ": ", e$message))
    results.all.mrpc[[i]] <- list(model=NA, Adj=NA, time = NA)
    times[i] <- NA
  })
  
  save(results.all.mrpc, file = paste0(path, 'results_all_MRPC-ADDIS.RData'))
  print(paste0('finished trio ', i))
}

save(results.all.mrpc, file = paste0(path, 'results_all_MRPC-ADDIS.RData'))


results.gmac = list()
times.gmac = NULL
for(i in 1:length(trios.with.pcs3)){
  X = as.matrix(trios.with.pcs3[[i]][,1:3])
  U = as.matrix(trios.with.pcs3[[i]][,-c(1:3)])
  start.time=Sys.time()
  results.gmac[[i]] = gmacOneTrio(trio = X, confounders = U, nominal.p = T, nperm = 1000)
  end.time=Sys.time()
  times.gmac[i] = difftime(end.time, start.time, units = 'secs')
  #print(times.gmac)
}

save(results.gmac, file = paste0(path, 'results_all_gmac.RData'))
save(times.gmac, file = paste0(path, 'results_all_times.RData'))

dat2 = cbind.data.frame(`Number of Confounders` = mrpc.params$number.of.Uconfs[sort_ix],
                       `Time To Compute Trio (min)` = times[sort_ix],
                       Model = mrpc.params$model[sort_ix])


times.mrgn = NULL
for(i in 1:length(trios.with.pcs3)){
  start.time = Sys.time()
  infer.trio(trios.with.pcs3[[i]])
  end.time = Sys.time()
  times.mrgn[i] = difftime(end.time, start.time, units = 'secs')
}
save(times.mrgn, file = paste0(path, 'results_all_times_mrgn.RData'))

# Plot shows one trio (#228) took an extremely long time to compute (over 4 hours to compute)
# Additional trios also took more than an hour
P2 = ggplot(data = dat2, aes(x = `Number of Confounders`, y = `Time To Compute Trio (min)`))+
  geom_point(size = 2.5)+
  geom_smooth()+
  theme_hc()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'top', legend.text = element_text(size = 16),axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 16),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5))+
  ylab('Time To Compute Trio (min)')+
  ggtitle("Time to compute simulated trios with MRPC-LOND", subtitle = 'Scenario with 1 - 15 confounders')

pdf(paste0(path, 'mrpc_LOND_time_to_compute_trios.pdf'))
plot(P2)
dev.off()

