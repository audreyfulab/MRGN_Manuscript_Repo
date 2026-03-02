# NOTE: Set your working directory to the repository root before running this script.
# NOTE: Originally run on HPC cluster. Some external data paths may need adjustment.

#rm(list=ls())
source("adapted_GMAC_func/GMAC_moded/R/GMAC.R")
library('MRGN')
library("missMDA")
tissue.name = "WholeBlood"

#--------function to assemble data for GMAC function---------
assemble.tables=function(trio.table, PCmat, kc, seed=222, which.imp='mice'){

  trio.ref=matrix(colnames(trio.table),nrow=length(colnames(trio.table))/3, ncol=3, byrow=T)
  print(trio.ref[1:10,])

  SNPs=trio.table[,c(1:(dim(trio.table)[2]/3))*3-2]

  SNP.unique=as.data.frame(SNPs[,match(unique(colnames(SNPs)), colnames(SNPs))])

  if(length(attr(na.omit(SNP.unique), 'na.action')) > 0 ){

    if(which.imp=='mice'){

      SNP.unique=complete(mice(SNP.unique, defaultMethod="polyreg"))
      print("...imputation..complete")

    }else{

      SNP.factors=apply.cfactor(SNP.unique)
      imp=imputeMCA(SNP.factors, seed=seed)
      SNP.unique=as.data.frame(apply(imp$completeObs,2,as.numeric))
      print("...imputation..complete")
      print(SNP.unique[1:5,1:5])

    }


  }


  expression.mat = trio.table[, -(c(1:(dim(trio.table)[2]/3))*3-2)]
  express.unique = expression.mat[, match(unique(colnames(expression.mat)), colnames(expression.mat))]


  trio.idx=build.triomap(trio.name.list = trio.ref,
                         express = express.unique,
                         genotype = SNP.unique)

  colnames(trio.idx)=c('snp', 'cis', 'trans')
  colnames(trio.ref)=c('snp', 'cis', 'trans')

  #assemble into list and transpose necessary tables so samples are in columns
  gmac.list=list(kc, t(PCmat), t(express.unique), t(SNP.unique), trio.idx, trio.ref)
  names(gmac.list)=c("known.conf", "cov.pool", "exp.dat", "snp.dat.cis", "trios.idx", 'trio.ref')


  return(gmac.list)

}





# helper which converts all columns of a matrix to factors
apply.cfactor=function(dataframe){

  cnames=colnames(dataframe)

  for(i in 1:dim(dataframe)[2]){

    dataframe[,i]=as.factor(dataframe[,i])

  }

  colnames(dataframe)=cnames

  return(dataframe)

}


#a function to construct the trio index map in the expression and genotype dataframes

build.triomap = function(trio.name.list=NA, express, genotype){

  index.mat=as.data.frame(matrix(0, nrow=dim(trio.name.list)[1], ncol = 3))

  for(i in 1:dim(trio.name.list)[1]){


    index.mat[i,1] = match(trio.name.list[i,1], colnames(genotype))
    index.mat[i,2] = match(trio.name.list[i,2], colnames(express))
    index.mat[i,3] = match(trio.name.list[i,3], colnames(express))

  }

  return(index.mat)

}

#read in the known/clinical covariates
kc.list = loadRData(file = "GTEx/data/kclist_top5_tiss.RData")

print(paste0("reading in data for ",tissues.name, " ..."))

#read trio matrix of only protein coding and lncRNA trios
trios.pr.lr.only = loadRData(file = paste0("GTEx/data/all.data.unqiue.snps.pclrna.only.",tissues.name,".RData"))
#get the total number of trios
number.of.trios = dim(trios.pr.lr.only)[2]/3
#matrix of trio indexes
trio.idx.mat = matrix(c(1:dim(trios.pr.lr.only)[2]), nrow = number.of.trios, ncol = 3, byrow = T)
#extract all trios into a list
just.trios = apply(trio.idx.mat, 1, function(x,y) as.data.frame(y[,x]), y = trios.pr.lr.only)
print(paste0("Number of trios to analyze for ", tissues.name, " is ", length(just.trios)))
#load in tissue specific PC matrix
pc.matrix=loadRData(paste("GTEx/data/PCs.matrix.", tissues.name, ".RData", sep = ""))

#read in tissue specific known confounders
kc = kc.list$WholeBlood
print("Prepping data for GMAC inference...")
print("Running imputation and assembling input list...")
input.list=assemble.tables(trio.table=trios.pr.lr.only, PCmat = pc.matrix, kc = kc, which.imp = 'miss')
print("done...")
save(input.list, file = paste0("GTEx/results/GMAC/GMAC_input_list_for_",tissues.name,".RData"))

# show input shapes
print(str(input.list))
print(dim(pc.matrix))
print(dim(trios.pr.lr.only))
print(dim(kc))

print("applying gmac with cis genes as mediators...")
output.cis.med <- gmac(known.conf = t(input.list$known.conf), cov.pool = input.list$cov.pool,
                        exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                        trios.idx = input.list$trios.idx, nperm = 1000,
                        nominal.p = TRUE, fdr = 0.05,
                        fdr_filter = 0.1)

print(str(output.cis.med))

print("applying gmac with trans genes as mediators...")
output.trans.med <- gmac(known.conf = t(input.list$known.conf), cov.pool = input.list$cov.pool,
                          exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                          trios.idx = input.list$trios.idx[,c(1,3,2)], nperm = 1000,
                          nominal.p = TRUE, fdr = 0.05,
                          fdr_filter = 0.1)

print(str(output.trans.med))
#results treating cis gene as mediator
out.table.cis=cbind.data.frame(output.cis.med[[1]], output.cis.med[[2]])
colnames(out.table.cis)=c(paste0('pval_', colnames(output.cis.med[[1]])),
                          paste0('effect_change_', colnames(output.cis.med[[2]])))

out.list.cis=list(out.table.cis, input.list, output.cis.med[[3]])
names(out.list.cis)=c("output.table", "input.list", "cov.indicator.list")

#results treating trans gene as mediator
out.table.trans=cbind.data.frame(output.trans.med[[1]], output.trans.med[[2]])
colnames(out.table.cis)=c(paste0('pval_', colnames(output.trans.med[[1]])),
                          paste0('effect_change_', colnames(output.trans.med[[2]])))

out.list.trans=list(out.table.trans, input.list, output.trans.med[[3]])
names(out.list.trans)=c("output.table", "input.list", "cov.indicator.list")



print("Done!...Saving...")

save(out.list.trans, file = paste0("GTEx/results/GMAC/GMAC_trans_inference_for_tissue_",tissues.name ,".RData"))
save(out.list.cis, file = paste0("GTEx/results/GMAC/GMAC_cis_inference_for_tissue_",tissues.name ,".RData"))
print("All done!")