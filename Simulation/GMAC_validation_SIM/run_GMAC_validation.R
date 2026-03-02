## this scipt runs the GMAC algorithm on the simulated data sets with all types of confounding variables
## created by sim_data_all_mods_all_confs_from_graph.R and saves the inference
rm(list=ls())

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

source("./adapted_GMAC_func/GMAC_moded/R/GMAC.R")
d.vec = c(2,4,6,8,10,12,15)
valid.res.tab.at.01 = as.data.frame(matrix(0, nrow = length(d.vec), ncol = 4))
valid.res.tab.at.05 = as.data.frame(matrix(0, nrow = length(d.vec), ncol = 4))
colnames(valid.res.tab.at.01) = colnames(valid.res.tab.at.05) = c("Simulated Confs", "Type I Error", "Type II Error", "Power")
for(i in 1:length(d.vec)){
  sim.datasets.all.fields=loadRData(file = paste0(file = "./GMAC_validation_SIM/data/gmac_valid_sim_simdata_numconfs_",d.vec[i], ".RData"))
  sim.datasets = lapply(sim.datasets.all.fields, function(x) x$data)
  conf.mat = loadRData(file = paste0(file = "./GMAC_validation_SIM/data/gmac_valid_sim__conf_mat_numconfs_",d.vec[i],".RData"))
  params = loadRData(paste0(file = "./GMAC_validation_SIM/data/gmac_valid_sim_parameters_numcomfs_",d.vec[i],".RData"))
  
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
  
  input.list=list(kc = t(cbind(rnorm(350), rnorm(350), rnorm(350))), cov.pool=t(conf.mat), 
                  snp.dat.cis=t(snp.dat.cis), exp.dat=t(exp.dat), trio.indexes=trios.idx,
                  nperm=1000, use.nominal.p=TRUE, fdr = 0.05, fdr_filter = 0.1)
  
  
  print("applying gmac with cis mediators")
  output.cis.med <- gmac(known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                         exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                         trios.idx = input.list$trio.indexes, #nperm = input.list$nperm,
                         nominal.p = input.list$use.nominal.p, fdr = input.list$fdr,
                         fdr_filter = input.list$fdr_filter, nperm = 100)
  
  print(str(output.cis.med))
  
  inf.at.01.cutoff = output.cis.med$pvals[,2]<0.01
  inf.at.05.cutoff = output.cis.med$pvals[,2]<0.05
  
  table.at.01.cutoff = table(GMAC.Mediation = inf.at.01.cutoff, True.Model = params$model)
  table.at.05.cutoff = table(GMAC.Mediation = inf.at.05.cutoff, True.Model = params$model)
  
  #results table for cutoff of 0.01
  valid.res.tab.at.01$`Simulated Confs`[i] = d.vec[i]
  valid.res.tab.at.01$`Type I Error`[i] = table.at.01.cutoff[2,1]/sum(table.at.01.cutoff[,1]) 
  valid.res.tab.at.01$`Type II Error`[i] = table.at.01.cutoff[1,2]/sum(table.at.01.cutoff[,2]) 
  valid.res.tab.at.01$Power[i] = table.at.01.cutoff[2,2]/sum(table.at.01.cutoff[,2]) 
  
  #results table for cutoff of 0.05
  valid.res.tab.at.05$`Simulated Confs`[i] = d.vec[i]
  valid.res.tab.at.05$`Type I Error`[i] = table.at.05.cutoff[2,1]/sum(table.at.05.cutoff[,1]) 
  valid.res.tab.at.05$`Type II Error`[i] = table.at.05.cutoff[1,2]/sum(table.at.05.cutoff[,2]) 
  valid.res.tab.at.05$Power[i] = table.at.05.cutoff[2,2]/sum(table.at.05.cutoff[,2]) 

  out.table.cis=cbind.data.frame(output.cis.med[[1]], output.cis.med[[2]])
  colnames(out.table.cis)=c(paste0('pval_', colnames(output.cis.med[[1]])),
                            paste0('effect_change_', colnames(output.cis.med[[2]])))
  
  out.list.cis=list(out.table.cis, input.list, output.cis.med[[3]])
  names(out.list.cis)=c("output.table", "input.list", "cov.indicator.list")
  
  print("Done!...Saving...")
  save(out.list.cis, file = paste0("./GMAC_validation_SIM/results/gmac_valid_sim_cis_results_numconfs_",d.vec[i],".RData"))
  
}

final.table = cbind.data.frame(Cutoff = c(rep("alpha < 0.05", length(d.vec)), rep("alpha < 0.01", length(d.vec))),
                               rbind.data.frame(valid.res.tab.at.05, valid.res.tab.at.01))

save(final.table, file = "./GMAC_validation_SIM/gmac_valid_summary_results_table.RData")
write.csv(final.table, file = "./GMAC_validation_SIM/gmac_valid_summary_results_table.csv")

