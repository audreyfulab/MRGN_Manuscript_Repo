# NOTE: Set your working directory to the repository root before running this script.
source("adapted_GMAC_func/GMAC_moded/R/GMAC.R")
library('MRPC')
tissue.names=read.csv("GTEx/data/tissuenames.csv", header = T)
library('MRGN')
wholeblood.num = 48
tissue=tissue.names[wholeblood.num, 2]

#function for MRPC
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

#read in the known confounders
kc.list = loadRData(file = "GTEx/data/kclist_top5_tiss.RData")
#preallocate results table:
res.table=as.data.frame(matrix(0, nrow = length(tissues.vec), ncol = 19))
colnames(res.table)=c("Tissue","M0","M1","M2","M3","M4","Other",
                      "%M0","%M1","%M2","%M3","%M4","%Other",
                      "Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")

med.table=as.data.frame(matrix(0, nrow = length(tissues.vec), ncol = 11))
colnames(med.table)=c("Tissue","M1.1", "M1.2", "%M1.1", "%M1.2",
                      "Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")

inf.vec=vector("list", length = length(tissues.vec))
reg.res.L=vector("list", length = length(tissues.vec))

names(inf.vec)=tissues.vec
names(reg.res.L)=tissues.vec

for(t in 1:length(tissues.vec)){

  print(paste0("reading in data for ",tissues.vec[t], " ..."))

  res.table$Tissue[t]=tissues.vec[t]
  med.table$Tissue[t]=tissues.vec[t]
  #read trio matrix of only protein coding and lncRNA trios
  trios.pr.lr.only = loadRData(file = paste0("GTEx/data/trios_subset_data/all.data.unqiue.snps.pclrna.only.",tissues.vec[t],".RData"))
  #get the total number of trios
  number.of.trios = dim(trios.pr.lr.only)[2]/3
  #matrix of trio indexes
  trio.idx.mat = matrix(c(1:dim(trios.pr.lr.only)[2]), nrow = number.of.trios, ncol = 3, byrow = T)
  #extract all trios into a list
  just.trios = apply(trio.idx.mat, 1, function(x,y) as.data.frame(y[,x]), y = trios.pr.lr.only)
  print(paste0("Number of trios to analyze for ", tissues.vec[t], " is ", length(just.trios)))
  #load in tissue specific PC matrix
  # NOTE: External GTEx V8 data path - update for your environment
  pc.matrix=loadRData(paste("GTEx/external_data/GTEx_version_8/",tissues.vec[t],
                            "_AllPC/PCs.matrix.",tissues.vec[t],".RData", sep = ""))
  #read in tissue specific known confounders
  kc = kc.list[[t]]
  #get the conf list
  print("calculating confounders for MRPC and MRGN...")
  
  # conduct confounder selection and save
  confs.out = get.conf.trios(trios = just.trios, cov.pool = pc.matrix, selection_fdr = 0.05,
                             filter_int_child = T, filter_fdr = 0.1, adjust_by = "all")
  save(confs.out, file = paste0("GTEx/results/MRGN/Regression_conf_select_res_",tissues.vec[t], ".RData"))

  # load precomputed confounders for faster analysis
  confs.out = loadRData(file = paste0("GTEx/results/MRGN/Regression_conf_select_res_",tissues.vec[t], ".RData"))
  #extract the confounders:
  conf.list=lapply(confs.out$sig.asso.covs, function(x,y){y[,x]}, y=pc.matrix)
  #bind trios with selected PCs
  trios.with.pcs=mapply(cbind.data.frame, just.trios, conf.list)
  #bind trios with pcs with the known confs:
  trios.with.pcs2=lapply(trios.with.pcs, function(x,y) cbind.data.frame(x,y), y = kc)
  #get the mean number of PCs
  save(trios.with.pcs2, file = paste0("GTEx/results/MRGN/", tissues.vec[t],"pr.lnc.data.with.pcs.RData"))
  num.pcs=unlist(lapply(trios.with.pcs, function(x){ncol(x)-3}))
  print(summary(num.pcs))
  res.table[t,14:19]=summary(num.pcs)

  #updating tables
  save(num.pcs, file = paste0("GTEx/results/MRGN/plots/num.pcs.wo.pseudo.",tissues.vec[t],".RData"))
  pdf(paste0("GTEx/results/MRGN/plots/num.pcs.wo.pseudo.",tissues.vec[t],".pdf"))
  hist(num.pcs, main=paste0("Selected Cov distribution for ",tissues.vec[t]), xlab="Number Selected")
  dev.off()

  #Run MRGN
  print(paste0("running MRGN on remaining trios: ",tissues.vec[t]))
  out=sapply(trios.with.pcs2, infer.trio, nperms=1000)
  print("done!....updating tables:")
  reg.res.L[[t]]=as.data.frame(out[-14,])
  #model class
  inf.mods=unlist(out[14,])
  inf.vec[[t]]=inf.mods
  #get num of trios for specifically M1.1 and M1.2
  idx=which(inf.mods=="M1.1" | inf.mods=="M1.2")
  med.num.pcs=unlist(lapply(trios.with.pcs2[idx], function(x){ncol(x)-3}))
  med.table[t,6:11]=summary(med.num.pcs)
  #check results
  inf.types=convert.cats(inf.mods)
  u=summary(as.factor(inf.types))
  res.table[t,2:13]=c(u, u/sum(u))
  #allocate mediation counts and percent to table
  u2=summary(as.factor(inf.mods[idx]))
  med.table[t,2:5]=c(u2, u2/sum(u2))

  print(paste0("saving results for ",tissues.vec[t]))
  save(reg.res.L, file = paste0("GTEx/results/MRGN/reg.res.wo.pseudo.RData"))
  save(inf.vec, file = paste0("GTEx/results/MRGN/inf.mods.wo.pseudo.RData"))
  save(res.table, file = paste0("GTEx/results/GMAC/results.table.",tissues.vec[t],".RData"))
  write.csv(res.table,  paste0("GTEx/data/TrioTables/MRGN.summary.results.table.csv"),
            row.names = T)
  write.csv(med.table,  paste0("GTEx/data/TrioTables/MRGN.mediation.results.table.csv"),
            row.names = T)

  
  # Apply MRPC 
  print(paste0("running MRPC on remaining trios: ",tissues.vec[t]))
  mrpc.infer.list=list()
  for(i in 1:length(trios.with.pcs2)){
    mrpc.infer.list[[i]]=apply.mrpc(trios.with.pcs2[[i]])
  }
  print("Done!...Saving results...")

  # save
  save(mrpc.infer.list, file = paste0("GTEx/results/MRPC/mrpc_inference_for_tissue_",tissues.vec[t] ,".RData"))

  
  # Apply GMAC Method
  print("applying gmac with cis genes as mediators...")
  output.cis.med <- gmac(known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                         exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                         trios.idx = input.list$trio.indexes, nperm = input.list$nperm,
                         nominal.p = input.list$use.nominal.p, fdr = input.list$fdr,
                         fdr_filter = input.list$fdr_filter)

  print(str(output.cis.med))

  print("applying gmac with trans genes as mediators...")
  output.trans.med <- gmac(known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                           exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                           trios.idx = input.list$trio.indexes[,c(1,3,2)], nperm = input.list$nperm,
                           nominal.p = input.list$use.nominal.p,fdr = input.list$fdr,
                           fdr_filter = input.list$fdr_filter)
  print(str(output.trans.med))

  # reorganize output for saving
  # cis results
  out.table.cis=cbind.data.frame(output.cis.med[[1]], output.cis.med[[2]])
  colnames(out.table.cis)=c(paste0('pval_', colnames(output.cis.med[[1]])),
                            paste0('effect_change_', colnames(output.cis.med[[2]])))

  out.list.cis=list(out.table.cis, input.list, output.cis.med[[3]])
  names(out.list.cis)=c("output.table", "input.list", "cov.indicator.list")

  # trans results
  out.table.trans=cbind.data.frame(output.trans.med[[1]], output.trans.med[[2]])
  colnames(out.table.cis)=c(paste0('pval_', colnames(output.trans.med[[1]])),
                            paste0('effect_change_', colnames(output.trans.med[[2]])))

  out.list.trans=list(out.table.trans, input.list, output.trans.med[[3]])
  names(out.list.trans)=c("output.table", "input.list", "cov.indicator.list")

  print("Done!...Saving...")

  save(out.list.trans, file = paste0("GTEx/results/GMAC/GMAC_trans_inference_for_tissue_",tissues.vec[t] ,".RData"))
  save(out.list.cis, file = paste0("GTEx/results/GMAC/GMAC_cis_inference_for_tissue_",tissues.vec[t] ,".RData"))

}
