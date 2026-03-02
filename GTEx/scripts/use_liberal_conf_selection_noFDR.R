# NOTE: Set your working directory to the repository root before running this script.
# This script is used to apply the more liberal procedure of confounder selection in which no FDR is applied to
# selected confounders. Instead confounders are selected at the alpha <0.01 cutoff:::: last updated : 6/5/2023
rm(list = ls())

#load in protein coding and lncRNA trios
#source("adapted_GMAC_func/GMACpostproc.R")
tissue.names=read.csv("GTEx/data/tissuenames.csv", header = T)
top5=c(48,40,33,6,1)
#get the top five tissues by sample size
#tissues.vec=tissue.names[top5, 2]
tissues.vec=tissue.names[top5, 2]
#library('psych')
#library('tmvtnorm')
#library('ff')
#library('propagate')
library('MRGN')

# Load known confounders
kclist = loadRData('GTEx/data/kclist_top5_tiss.RData')

for(i in 1:length(tissues.vec)){
  #read in subset of trios pertaining to only protein coding and LncRNA genes
  message("Begining confounder selection for tissue: ",tissues.vec[i]," ...")
  all.unq.snps.ldf.pl.only = loadRData(file = paste0("GTEx/data/trios_subset_data/all.data.unqiue.snps.pclrna.only.",tissues.vec[i],".RData"))
  # NOTE: External GTEx V8 data path - update for your environment
  PCs.matrix=loadRData(paste("GTEx/external_data/GTEx_version_8/",tissues.vec[i],
                            "_AllPC/PCs.matrix.",tissues.vec[i],".RData", sep = ""))
  numcols = dim(all.unq.snps.ldf.pl.only)[2]
  idx.mat = cbind(snp.idx = seq(1, numcols-2, 3),
                  gene1.idx = seq(2, numcols-1, 3),
                  gene2.idx = seq(3, numcols, 3))
  #extract trios into list
  trioslist = lapply(c(1:(numcols/3)), function(x,y,z) y[, z[x,]], y = all.unq.snps.ldf.pl.only,
                 z = idx.mat)
  
  #perform confounder selection - no correction, selected at alpha < 0.05, with filtering 
  # alpha at 0.05 because the pool size is only 600 (in our simulations the pool size was over 10,000 e.g much larger)
  # conf.out=get.conf.trios(trios = trioslist, cov.pool = PCs.matrix, blocksize = 200, filter_int_child = T, 
  #                         selection_fdr = 0.05, filter_fdr = 0.1, alpha = 0.05, adjust_by = 'none', save.list = F)
  conf.out=get.conf.trios(trios = trioslist, cov.pool = PCs.matrix, blocksize = 200, filter_int_child = T, 
                          selection_fdr = 0.05, filter_fdr = 0.1, alpha = 0.01, adjust_by = 'none', save.list = F)
  
  print('Finished Confounder selection...print first 5 sets of selected confs...')
  print(lapply(conf.out$sig.asso.covs[1:5], head))
  
  #extract the confounders into a list of dataframes
  confs.in.list = lapply(conf.out$sig.asso.covs, function(x,y) y[,x], y = PCs.matrix)
  
  #combine each trio with its selected confounders into a list of dataframes
  final.trios.with.confs = lapply(c(1:(numcols/3)), function(w,x,y,z) cbind.data.frame(x[[w]], y, z[[w]]),
                                  x = trioslist, y = kclist[[i]], z = confs.in.list)
  
  print(lapply(final.trios.with.confs[1:5], function(x) if (dim(x)[2]>=7) head(x[,1:7]) else head(x)))
  
  message('done!...Now applying infer.trio')
  
  # apply MRGN 
  reg.res = sapply(final.trios.with.confs, infer.trio, use.perm = TRUE, gamma = 0.05, nperms = 1000,
                   alpha = 0.01)
  
  
  message('saving all results...')
  
  save(conf.out, file = paste0("GTEx/trios_data_prlnc_liberal_conf_sel/list_output/alpha01/conf.output.list.",tissues.vec[i],"alpha01.RData"))
  save(final.trios.with.confs, file = paste0("GTEx/trios_data_prlnc_liberal_conf_sel/data_with_confs/alpha01/trios.with.confs.",tissues.vec[i],"alpha01.RData"))
  save(reg.res, file = paste0("GTEx/trios_data_prlnc_liberal_conf_sel/infer_trio_results/alpha01/infer.trio.results.",tissues.vec[i],"alpha01.RData"))
  
  
  
  
}

























