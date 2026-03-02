
source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")

library('qvalue', lib="/mnt/ceph/jarredk/Rpackages")
library('MRPC', lib="/mnt/ceph/jarredk/Rpackages")
library('permute', lib="/mnt/ceph/jarredk/Rpackages")
library('prodlim', lib="/mnt/ceph/jarredk/Rpackages")
library('ppcor', lib="/mnt/ceph/jarredk/Rpackages")
library('car', lib = "/mnt/ceph/jarredk/Rpackages")

#library('Morpho', lib = "/mnt/ceph/jarredk/Rpackages")
library('FactoMineR', lib="/mnt/ceph/jarredk/Rpackages")


source("/mnt/ceph/jarredk/ADDIS_verify/ADDIS_Post_Analysis_processing.R")

top5=c(1, 6, 33, 40, 48)
path='/mnt/ceph/jarredk/Reg_Net/'
tissues.vec=tissue.names[top5, 2:3]




#function to run post processing on GMAC output

run.postproc=function(output.pvals=NULL, trio.ref=NULL){

  qvals=qvalue(output.pvals, fdr.level = 0.1)

  num.inferred=sum(qvals$significant)
  #print(qvals$significant)
  sig.trios=trio.ref[qvals$significant,]

  result=cbind.data.frame(qvals$qvalue[qvals$significant], output.pvals[qvals$significant])
  colnames(result)=c("Q.values","Out.p")


  return(list(total.inferred=num.inferred, sig.trios=sig.trios, pq.values=result))
}

# #load in output for the cis trios:
#
output.WB=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/WholeBlood/all_trios_output_cis.Rdata')
output.AS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/AdiposeSubcutaneous/all_trios_output_cis.Rdata')
output.ArT=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/ArteryTibial/all_trios_output_cis.Rdata')
output.MS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/MuscleSkeletal/all_trios_output_cis.Rdata')
output.SSE=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/SkinSunExposed/all_trios_output_cis.Rdata')
# #
# #post process GMAC-cis-mediated
# result.WB=run.postproc(output.pvals = output.WB$output.table[,5], trio.ref = output.WB$output.table[,1:3])
# result.AS=run.postproc(output.pvals = output.AS$output.table[,5], trio.ref = output.AS$output.table[,1:3])
# result.ArT=run.postproc(output.pvals = output.ArT$output.table[,5], trio.ref = output.ArT$output.table[,1:3])
# result.MS=run.postproc(output.pvals = output.MS$output.table[,5], trio.ref = output.MS$output.table[,1:3])
# result.SSE=run.postproc(output.pvals = output.SSE$output.table[,5], trio.ref = output.SSE$output.table[,1:3])
#
#
# #read in output for the trans trios:
#
# #load in output for the cis trios:
#
output.t.WB=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/WholeBlood/all_trios_output_trans.Rdata')
output.t.AS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/AdiposeSubcutaneous/all_trios_output_trans.Rdata')
output.t.ArT=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/ArteryTibial/all_trios_output_trans.Rdata')
output.t.MS=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/MuscleSkeletal/all_trios_output_trans.Rdata')
output.t.SSE=loadRData(fileName='/mnt/ceph/jarredk/GMACanalysis/SkinSunExposed/all_trios_output_trans.Rdata')
# #
# result.t.WB=run.postproc(output.pvals = output.t.WB$output.table[,5], trio.ref = output.t.WB$output.table[,1:3])
# result.t.AS=run.postproc(output.pvals = output.t.AS$output.table[,5], trio.ref = output.t.AS$output.table[,1:3])
# result.t.ArT=run.postproc(output.pvals = output.t.ArT$output.table[,5], trio.ref = output.t.ArT$output.table[,1:3])
# result.t.MS=run.postproc(output.pvals = output.t.MS$output.table[,5], trio.ref = output.t.MS$output.table[,1:3])
# result.t.SSE=run.postproc(output.pvals = output.t.SSE$output.table[,5], trio.ref = output.t.SSE$output.table[,1:3])
#
#
# result.table=cbind.data.frame(tissues.vec[,1],
#                               c(dim(output.AS$output.table)[1],
#                                 dim(output.ArT$output.table)[1],
#                                 dim(output.MS$output.table)[1],
#                                 dim(output.SSE$output.table)[1],
#                                 dim(output.WB$output.table)[1]),
#                               c(result.AS$total.inferred,
#                                 result.ArT$total.inferred,
#                                 result.MS$total.inferred,
#                                 result.SSE$total.inferred,
#                                 result.WB$total.inferred),
#                               c(result.t.AS$total.inferred,
#                                 result.t.ArT$total.inferred,
#                                 result.t.MS$total.inferred,
#                                 result.t.SSE$total.inferred,
#                                 result.t.WB$total.inferred),
#                               c(123, 90, 104, 122, 131),
#                               c(130, 113, 103, 123, 113))
#
#
# colnames(result.table)=c("Tissue",
#                          "Total.Num.Samples",
#                          "GMAC.Total.M1T1",
#                          "GMAC.Total.M1T2",
#                          "LOND.total.M1",
#                          "ADDIS.total.M1")
#
#
#
# write.csv(result.table, file="/mnt/ceph/jarredk/GMACanalysis/master_tables/final_result_table10000.csv")



#a function to match which trios are found in the AM1T1/T2 and LM1T1/T2 tables and which have
#inferred model types that differ

match.trios=function(tissues=tissues.vec[,1]){

  library('prodlim', lib="/mnt/ceph/jarredk/Rpackages")

  num.matched.cis=NULL
  num.matched.trans=NULL

  matched.cis=vector("list", length = length(tissues))
  matched.trans=vector("list", length = length(tissues))

  unmatched.cis=vector("list", length = length(tissues))
  unmatched.trans=vector("list", length = length(tissues))

  for(i in 1:length(tissues)){


    #read in GMAC results
    output.cis=loadRData(fileName=paste('/mnt/ceph/jarredk/GMACanalysis/', tissues[i], '/all_trios_output_cis.Rdata', sep = ""))
    output.trans=loadRData(fileName=paste('/mnt/ceph/jarredk/GMACanalysis/', tissues[i], '/all_trios_output_trans.Rdata', sep=""))

    #process
    postproc.cis=run.postproc(output.pvals = output.cis$output.table[,5],
                            trio.ref = output.cis$output.table[,1:3])
    #print(dim(postproc.cis$sig.trios))

    postproc.trans=run.postproc(output.pvals = output.trans$output.table[,5],
                              trio.ref = output.trans$output.table[,1:3])
    #print(dim(postproc.trans$sig.trios))

    m1t1=read.csv(file = paste("/mnt/ceph/jarredk/Reg_Net/AM1T1/", tissues[i], ".csv", sep = ""),
                  header = T, sep = ",")
    m1t2=read.csv(file = paste("/mnt/ceph/jarredk/Reg_Net/AM1T2/", tissues[i], ".csv", sep = ""),
                  header = T, sep = ",")


    num.matched.cis[i]=length(na.omit(row.match(postproc.cis$sig.trios,m1t1[,2:4])))
    #print(num.matched.cis)
    num.matched.trans[i]=length(na.omit(row.match(m1t2[,2:4],postproc.trans$sig.trios)))

    unmatched.cis[[i]]=postproc.cis$sig.trios[-c(na.omit(row.match(m1t1[,2:4],postproc.cis$sig.trios))),]
    unmatched.trans[[i]]=postproc.trans$sig.trios[-c(na.omit(row.match(m1t2[,2:4],postproc.trans$sig.trios))),]

    matched.cis[[i]]=postproc.cis$sig.trios[c(na.omit(row.match(m1t1[,2:4],postproc.cis$sig.trios))),]
    matched.trans[[i]]=postproc.trans$sig.trios[c(na.omit(row.match(m1t2[,2:4],postproc.trans$sig.trios))),]



  }

  names(unmatched.cis)=tissues
  names(unmatched.trans)=tissues

  names(matched.cis)=tissues
  names(matched.trans)=tissues

  table1=cbind.data.frame(tissues, num.matched.cis, num.matched.trans)
  colnames(table1)=c("tissue","cis.common", "trans.common")

  return(list(table=table1,
              unmatched.cis=unmatched.cis,
              unmatched.trans=unmatched.trans,
              matched.cis=matched.cis,
              matched.trans=matched.trans))


}












#match.trios(tissues=tissues.vec[,1], which.mrpc="ADDIS")
#match.trios(tissues=tissues.vec[,1], which.mrpc="LOND")



#a simple function to find the trio index of a trio based on the SNP name

locate.trio=function(input=NULL, trio.mat=NULL){

  index.vec=NULL

  for(j in 1:dim(input)[1]){

    ind1=which(colnames(trio.mat)==input[j,1])
    index.vec[j]=(ind1-1)/3+1

  }

  return(index.vec)

}

# a simple logical function passed to lapply in which.class()
fn=function(x,b){x==b}

#a simple function to find the addis classification of a trio

which.class=function(ind=NULL, classes.list=NULL){

  models=names(classes.list)
  class.vec=NULL

  for(i in 1:length(ind)){

    vec1=unlist(lapply(lapply(classes.list, fn, ind[i]),any))

    class.vec[i]=ifelse(any(vec1), models[vec1], "Other")

  }

  return(class.vec)

}



#cross analyze mediation test results with ADDIS trios
cross.analyze=function(tissues=tissues.vec[,1], save=FALSE, path="/mnt/ceph/jarredk/GMACanalysis/"){


  tables.list=vector("list", length = length(tissues))
  bk.tbl=as.data.frame(matrix(0, nrow=length(tissues), ncol=6))
  colnames(bk.tbl)=c("MO","M1","M2","M3","M4","Other")
  row.names(bk.tbl)=tissues

  unq.table=as.data.frame(matrix(0, nrow = length(tissues), ncol = 3))
  colnames(unq.table)=c("Unq.Cis", "Unq.Trans", "Both")
  row.names(unq.table)=tissues

  for(i in 1:length(tissues)){

    #read in GMAC results
    output.cis=loadRData(fileName=paste(path, tissues[i], '/all_trios_output_cis.Rdata', sep = ""))
    output.trans=loadRData(fileName=paste(path, tissues[i], '/all_trios_output_trans.Rdata', sep=""))

    #process
    postproc.cis=run.postproc(output.pvals = output.cis$output.table[,5],
                              trio.ref = output.cis$output.table[,1:3])
    #print(dim(postproc.cis$sig.trios))

    postproc.trans=run.postproc(output.pvals = output.trans$output.table[,5],
                                trio.ref = output.trans$output.table[,1:3])

    addis.classes=loadRData(fileName=paste("/mnt/ceph/jarredk/AddisReRunFiles/List.models.", tissues[i], ".all.RData" ,sep=""))


    file1=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
                tissues[i],"_AllPC/data.snp.cis.trans.final.",
                tissues[i],".V8.unique.snps.RData", sep = "")

    trios=loadRData(fileName=file1)

    d1=dim(postproc.cis$sig.trios)[1]
    d2=dim(postproc.trans$sig.trios)[1]

    unq.cis=which(is.na(row.match(postproc.cis$sig.trios, postproc.trans$sig.trios)))
    unq.trans=which(is.na(row.match(postproc.trans$sig.trios, postproc.cis$sig.trios)))

    if(isTRUE(d1>d2)){
      unq.common=na.omit(row.match(postproc.trans$sig.trios, postproc.cis$sig.trios))
      common=postproc.cis$sig.trios[unq.common,]
      perm.reg.p.common=postproc.cis$pq.values$Out.p[unq.common]
    }else{
      unq.common=na.omit(row.match(postproc.cis$sig.trios, postproc.trans$sig.trios))
      common=postproc.trans$sig.trios[unq.common,]
      perm.reg.p.common=postproc.trans$pq.values$Out.p[unq.common]
    }

    pr.p.unq.cis=postproc.cis$pq.values$Out.p[unq.cis]

    pr.p.unq.trans=postproc.trans$pq.values$Out.p[unq.trans]

    #categorize direction
    found=c(rep("Both", dim(common)[1]),
            rep("Cis.Med", dim(postproc.cis$sig.trios[unq.cis,])[1]),
            rep("Trans.Med", dim(postproc.trans$sig.trios[unq.trans,])[1]))

    #bind all information into a single table
    table.pp=rbind.data.frame(common, postproc.cis$sig.trios[unq.cis,], postproc.trans$sig.trios[unq.trans,])
    table2.pp=cbind.data.frame(table.pp,
                               Mediation.type=found,
                               Perm.reg.p=c(perm.reg.p.common, pr.p.unq.cis, pr.p.unq.trans))

    #get the trio indicies and Addis M class
    trio.idx=locate.trio(input=table2.pp[,1:3], trio.mat = trios)
    Mtype=which.class(ind = trio.idx, classes.list = addis.classes)

    #get final table
    table.final=cbind.data.frame(table2.pp, Trio.Num=trio.idx, Addis.Class=Mtype)

    tables.list[[i]]=table.final


    bk.tbl[i,]=c(dim(tables.list[[i]][which(tables.list[[i]]$Addis.Class=="MO"),])[1],
                 dim(tables.list[[i]][which(tables.list[[i]]$Addis.Class=="M1"),])[1],
                 dim(tables.list[[i]][which(tables.list[[i]]$Addis.Class=="M2"),])[1],
                 dim(tables.list[[i]][which(tables.list[[i]]$Addis.Class=="M3"),])[1],
                 dim(tables.list[[i]][which(tables.list[[i]]$Addis.Class=="M4"),])[1],
                 dim(tables.list[[i]][which(tables.list[[i]]$Addis.Class=="Other"),])[1])

    unq.table[i,]=c(dim(tables.list[[i]][which(tables.list[[i]]$Mediation.type=="Cis.Med"),])[1],
                    dim(tables.list[[i]][which(tables.list[[i]]$Mediation.type=="Trans.Med"),])[1],
                    dim(tables.list[[i]][which(tables.list[[i]]$Mediation.type=="Both"),])[1])

  }

  if(save==TRUE){

    save(tables.list, file = paste0(path, "GMAC.results.final.table.RData"))

  }

  bk.tbl$`Total GMAC Inferred`=rowSums(bk.tbl)

  return(list(final.tables=tables.list,
              breakdown.tab=bk.tbl,
              unique.tab=unq.table,
              postproc.cis=postproc.cis$sig.trios,
              postproc.trans=postproc.trans$sig.trios))

}


#plots the distribution of M types for non-overlapping trios classified as mediation in GMAC

plot.dist=function(class.vec.cis=NULL, class.vec.trans=NULL){

  D1=summary(as.factor(class.vec.cis))
  print(D1)
  D1.1=cbind.data.frame(names(D1),D1, rep("cis", length(D1)))
  colnames(D1.1)=c("Class", "Number of Trios", "mediator type")

  D2=summary(as.factor(class.vec.trans))
  print(D2)
  D2.1=cbind.data.frame(names(D2),D2, rep("trans", length(D2)))
  colnames(D2.1)=c("Class", "Number of Trios", "mediator type")

  all.data=rbind.data.frame(D1.1,D2.1)
  all.data$Class=as.factor(all.data$Class)
  all.data$`mediator type`=as.factor(all.data$`mediator type`)

  library(ggpubr)


  ggbarplot(all.data, x = "Class", y = "Number of Trios",
            fill = "mediator type", color = "mediator type",
            palette = c("cadetblue1","red4", "yellow","navyblue","salmon1", "chartreuse3"),
            label = FALSE, ylab = FALSE, position=position_dodge(0.9))+
    theme(axis.text.x = element_text(size=10, hjust=0.5,vjust=0.7),
          axis.text.y = element_text(size=9))+
    scale_y_continuous(breaks = seq(0, 500, 50), lim = c(0, 500))+
    coord_flip()+
    ggtitle("ADDIS Model Classification for Non-Overlapping M1 Trios",
            subtitle = "Inferred by GMAC analysis")


}







#a function to look at the regression on the trans gene in GMAC and ADDIS


cross.regress=function(tissue="WholeBlood", trio.ind=NULL, mod.type="cis", addis.pcs=NULL, plot.it=FALSE,
                       mtype="", verbose=TRUE){

  if(mod.type=="cis"){
    out.data=loadRData(fileName=paste0('/mnt/ceph/jarredk/GMACanalysis/', tissue, '/all_trios_output_cis.Rdata'))

  }else{

    out.data=loadRData(fileName=paste0('/mnt/ceph/jarredk/GMACanalysis/', tissue, '/all_trios_output_trans.Rdata'))

  }

  trio.loc=as.matrix(out.data$input.list$trios.idx)[trio.ind,]
  known.conf=t(out.data$input.list$known.conf)
  SNP=t(as.data.frame(out.data$input.list$snp.dat.cis)[trio.loc[1],])
  GE=t(as.data.frame(out.data$input.list$exp.dat)[trio.loc[-1],])
  which.covs=which(out.data$cov.indicator.list[trio.ind,]==1)
  sel.cov.pool=as.data.frame(out.data$input.list$cov.pool)[which.covs,]

  print("-----------GMAC-Selected-PCs--------------")
  print(row.names(sel.cov.pool))

  sel.cov.pool=t(sel.cov.pool)

  data.mat=cbind.data.frame(GE, SNP, sel.cov.pool, known.conf)
  if(length(addis.pcs)>1){
    which.covs.addis=match(addis.pcs, row.names(out.data$input.list$cov.pool))
    addis.data=cbind.data.frame(GE, SNP, t(out.data$input.list$cov.pool[which.covs.addis,]))
  }else if(length(addis.pcs)==1){
    which.covs.addis=which(row.names(out.data$input.list$cov.pool)==addis.pcs)
    addis.data=cbind.data.frame(GE, SNP, out.data$input.list$cov.pool[which.covs.addis,])
  }else{
    addis.data=cbind.data.frame(GE, SNP)
  }

  if(mod.type=="cis"){

    colnames(data.mat)=c("cis.gene", "trans.gene", "SNP", colnames(sel.cov.pool), colnames(known.conf))
    colnames(addis.data)=c("cis.gene", "trans.gene", "SNP", addis.pcs)

  }else{

    colnames(data.mat)=c("trans.gene", "cis.gene", "SNP", colnames(sel.cov.pool), colnames(known.conf))
    colnames(addis.data)=c("trans.gene", "cis.gene", "SNP", addis.pcs)

  }


  data.mat2=data.mat
  addis.data2=addis.data
  #convert.factors
  #data.mat$SNP=as.factor(data.mat$SNP)
  data.mat$pcr=as.factor(data.mat$pcr)
  data.mat$sex=as.factor(data.mat$sex)
  data.mat$platform=as.factor(data.mat$platform)
  #addis.data$SNP=as.factor(addis.data$SNP)
  if(mod.type=="cis"){
    model.GMAC=lm(trans.gene~., data = data.mat)
    model.ADDIS=lm(trans.gene~., data = addis.data)
  }else{
    model.GMAC=lm(cis.gene~., data = data.mat)
    model.ADDIS=lm(cis.gene~., data = addis.data)
  }

  if(plot.it==TRUE){
    png(paste0("/mnt/ceph/jarredk/GMACanalysis/additional_plots/Reg_assump_plots/GMAC/",
               "GMAC",mod.type,mtype,trio.ind,".png"))
    par(mfrow=c(2,2))
    plot(model.GMAC)
    par(mfrow=c(1,1))
    dev.off()
  }

  if(verbose==TRUE){
    print("--------------------------GMAC------------------------------")
    print(summary(model.GMAC))
  }

  if(plot.it==TRUE){
    png(paste0("/mnt/ceph/jarredk/GMACanalysis/additional_plots/Reg_assump_plots/MRPC/",
               "ADDIS",mod.type,mtype,trio.ind,".png"))
    par(mfrow=c(2,2))
    plot(model.ADDIS)
    par(mfrow=c(1,1))
    dev.off()
  }

  if(verbose==TRUE){
    print("--------------------------ADDIS------------------------------")
    print(summary(model.ADDIS))
  }

  #store relevant regression info
  model.matrix.addis=model.matrix(model.ADDIS)
  model.matrix.gmac=model.matrix(model.GMAC)
  b.gmac=summary(model.GMAC)$coefficients[,1]
  b.addis=summary(model.ADDIS)$coefficients[,1]
  sigma.gmac=summary(model.GMAC)$sigma
  sigma.addis=summary(model.ADDIS)$sigma

  if(mod.type=="cis"){

    p.addis=as.data.frame(summary(model.ADDIS)$coefficients)$`Pr(>|t|)`[which(row.names(as.data.frame(summary(model.ADDIS)$coefficients))=="cis.gene")]
    p.gmac=as.data.frame(summary(model.GMAC)$coefficients)$`Pr(>|t|)`[which(row.names(as.data.frame(summary(model.GMAC)$coefficients))=="cis.gene")]

    bgcg=b.gmac[which(names(b.gmac)=="cis.gene")]
    bacg=b.addis[which(names(b.addis)=="cis.gene")]

  }else{

    p.addis=as.data.frame(summary(model.ADDIS)$coefficients)$`Pr(>|t|)`[which(row.names(as.data.frame(summary(model.ADDIS)$coefficients))=="trans.gene")]
    p.gmac=as.data.frame(summary(model.GMAC)$coefficients)$`Pr(>|t|)`[which(row.names(as.data.frame(summary(model.GMAC)$coefficients))=="trans.gene")]
    bgcg=b.gmac[which(names(b.gmac)=="trans.gene")]
    bacg=b.addis[which(names(b.addis)=="trans.gene")]
  }



  return(list(addis=addis.data2,
              GMAC=data.mat2,
              Regress=list(X.addis=model.matrix.addis,
                           X.gmac=model.matrix.gmac,
                           b.gmac=b.gmac,
                           b.addis=b.addis,
                           sigma.gmac=sigma.gmac,
                           sigma.addis=sigma.addis,
                           p.gmac=p.gmac,
                           p.addis=p.addis,
                           b.gmac.cis.gene=bgcg,
                           b.addis.cis.gene=bacg)))



}









                            #PERMUTED REG
#------------------------------------------------------------------------------

#A helper function to run.permuted.reg
#for parallelization of perm.reg

help.fn1=function(perm.map=NULL, Trio=NULL, Alg=NULL, med.type=NULL){

  #conditioning on each genotype
  # trio0=trio[which(trio$SNP==0),]
  # trio1=trio[which(trio$SNP==1),]
  # trio2=trio[which(trio$SNP==2),]
  #print(trio1$cis.gene)
  #print(trio2$cis.gene)
  #print(trio0$cis.gene[shuffle(trio0$cis.gene)])
  #permute



  trio.permuted=Trio


  if(med.type=="cis"){

    trio.permuted$cis.gene=trio.permuted$cis.gene[perm.map]

  }else{

    trio.permuted$trans.gene=trio.permuted$trans.gene[perm.map]

  }


  if(Alg=="GMAC"){
    trio.permuted$pcr=as.factor(trio.permuted$pcr)
    trio.permuted$sex=as.factor(trio.permuted$sex)
    trio.permuted$platform=as.factor(trio.permuted$platform)
  }






  if(med.type=="cis"){

    coef.mat=as.data.frame(summary(lm(trans.gene~., data=trio.permuted))$coefficients)
    wald.stat=coef.mat$`t value`[which(row.names(coef.mat)=="cis.gene")]
    #print(coef.mat)

  }else{

    coef.mat=as.data.frame(summary(lm(cis.gene~., data=trio.permuted))$coefficients)
    wald.stat=coef.mat$`t value`[which(row.names(coef.mat)=="trans.gene")]
    #print(coef.mat)

  }

  return(wald.stat)

}




#-------------------------------------------------------

#a function to recreate the permutation regression from GMAC

run.permuted.reg=function(trio, nperms=1000, plot=FALSE, filename=NULL, Alg="GMAC", med.type=NULL){

  #allocation of space
  wald.stat=NULL
  cmap=rep(0, dim(trio)[1])
  mediator_perm=matrix(c(1:dim(trio)[1]), nrow=dim(trio)[1], ncol = nperms)


  #get observed statisitc
  if(sum(med.type==c("Cis.Med", "Both"))>0){
    med.type="cis"
  }else{
    med.type="trans"
  }

  trio2=trio
  #convert.factors
  #trio2$SNP=as.factor(trio2$SNP)
  if(Alg=="GMAC"){
    trio2$pcr=as.factor(trio2$pcr)
    trio2$sex=as.factor(trio2$sex)
    trio2$platform=as.factor(trio2$platform)
  }


  if(med.type=="cis"){

    coef.mat=as.data.frame(summary(lm(trans.gene~., data=trio2))$coefficients)
    test.wald=coef.mat$`t value`[which(row.names(coef.mat)=="cis.gene")]

  }else{

    coef.mat=as.data.frame(summary(lm(cis.gene~., data=trio2))$coefficients)
    test.wald=coef.mat$`t value`[which(row.names(coef.mat)=="trans.gene")]

  }


  #preallocate all permutations
  for (j in 0:2) {
    ind <- which(trio$SNP == j)
    if (length(ind) > 1) {
      mediator_perm[ind, ] <- apply(mediator_perm[ind, ], 2, sample)
    }
  }

  #print(mediator_perm)

  #apply permutation regression
  wald.stat=apply(mediator_perm, 2, help.fn1, Trio=trio2, Alg=Alg, med.type=med.type)

  if(plot==TRUE){


    #get colormap for plotting
    for (k in 0:2) {
      ind <- which(trio2$SNP == j)
      if (length(ind) > 1) {
        cmap[ind, ] <- k
      }
    }


    png(filename = filename)

    par(mfrow=c(2,2))
    plot(trio.permuted$cis.gene,
         trio.permuted$trans.gene,
         xlab="Permuted Cis Gene",
         ylab="Trans Gene",
         col=as.factor(cmap),
         bg=as.factor(cmap),
         pch=21)

    legend("topleft", c("Gt=0", "Gt=1","Gt=2"), col = c(1,2,3), pch=c(21,21,21), fill=c(1,2,3))

    plot(trio$cis.gene,
         trio$trans.gene,
         xlab="Cis Gene",
         ylab="Trans Gene",
         col=as.factor(trio$SNP),
         bg=as.factor(trio$SNP),
         pch=21)

    hist(wald.stat)

    par(mfrow=c(1,1))

    dev.off()



  }

  p.value=2*min(sum(wald.stat<test.wald)/length(wald.stat),sum(wald.stat>test.wald)/length(wald.stat) )

  nominal.p.value=2 * (1 - pnorm(abs((test.wald - mean(wald.stat))/sd(wald.stat))))


  return(list(p.value=p.value, nominal.p=nominal.p.value, null.wald.stats=wald.stat, obs.wald.stat=test.wald))

}



#simple function used in check.ns.gmac() for coefficient
#simulation matrix to longform for plotting

longform=function(iters=NULL,varvec1=NULL, cutdata=NULL){

  var1=rep(varvec1, iters)

  long=cbind.data.frame(Coef=var1, Value=as.vector(as.matrix(cutdata)))
  return(long)

}

#a function to check the numerical stability of the regression using all adaptively selected pcs
#in GMAC, and using the coefficients and residual standard error obtained from the regression of a
#specific trio (the output of cross.regress). The trans gene is simulated by randomly sampling the noise
#at each iteration and running the regression on the simulated trans gene

#output: indicator matrix for which regressions had significant SNP (col 1) and significant cis.gene (col 2)

check.ns.gmac=function(input.list.data=NULL, iters=1, alpha=0.01, print.it=FALSE, which.mod="GMAC",
                       plot.it=FALSE, save.name=NULL){

  sig.snp=NULL
  sig.cis=NULL
  if(which.mod=="GMAC"){

    X=as.matrix(input.list.data$Regress$X.gmac)
    b=input.list.data$Regress$b.gmac
    sigma=input.list.data$Regress$sigma.gmac

  }else{

    X=as.matrix(input.list.data$Regress$X.addis)
    b=input.list.data$Regress$b.addis
    sigma=input.list.data$Regress$sigma.addis

  }

  simu.coefs=as.data.frame(matrix(0, nrow=length(b), ncol=iters))
  row.names(simu.coefs)=names(b)
  colnames(simu.coefs)=paste0("iter", c(1:iters))

  for(i in 1:iters){

    if(which.mod=="GMAC"){
      new.data=input.list.data$GMAC
    }else{
      new.data=input.list.data$addis
    }

    trans.gene.new=X%*%b+rnorm(dim(X)[1], mean = 0, sd = sigma)
    #print(trans.gene.new)

    new.data$trans.gene=trans.gene.new
    if(which.mod=="GMAC"){
      new.data$pcr=as.factor(new.data$pcr)
      new.data$sex=as.factor(new.data$sex)
      new.data$platform=as.factor(new.data$platform)
      #new.data$SNP=as.factor(new.data$SNP)
    }

    reg=summary(lm(trans.gene~., data = new.data))
    if(print.it==TRUE){print(reg)}
    coef=as.data.frame(reg$coefficient)
    p.snp=coef$`Pr(>|t|)`[3]
    p.cis=coef$`Pr(>|t|)`[2]

    sig.snp[i]=ifelse(p.snp<alpha, 1, 0)
    sig.cis[i]=ifelse(p.cis<alpha, 1, 0)

    simu.coefs[,i]=coef$Estimate

  }


  if(plot.it==TRUE){


    library('ggplot2', lib="/mnt/ceph/jarredk/Rpackages")
    library('ggridges', lib="/mnt/ceph/jarredk/Rpackages")

    if(which.mod=="GMAC"){
      rm.c=match(c("(Intercept)", "pcr1", "platform1","sex2"),row.names(simu.coefs))
    }else{
      rm.c=match(c("(Intercept)"),row.names(simu.coefs))
    }

    long.data=longform(iters = iters,
                       varvec1 = row.names(simu.coefs)[-c(rm.c)],
                       cutdata = simu.coefs[-c(rm.c),])

    act.data=cbind.data.frame(Coef=names(b)[-rm.c], Value=b[-rm.c])

    if(which.mod=="GMAC"){
      fn1=paste0("/mnt/ceph/jarredk/GMACanalysis/additional_plots/GMAC.ridgeline.trio.", save.name, ".png")
      fn2=paste0("/mnt/ceph/jarredk/GMACanalysis/additional_plots/GMAC.boxplot.trio.", save.name, ".png")
    }else{
      fn1=paste0("/mnt/ceph/jarredk/GMACanalysis/additional_plots/ADDIS.ridgeline.trio.", save.name, ".png")
      fn2=paste0("/mnt/ceph/jarredk/GMACanalysis/additional_plots/ADDIS.boxplot.trio.", save.name, ".png")
    }

    png(fn1)

    p=ggplot(long.data, aes(x = Value, y = Coef)) +
      geom_density_ridges() +
      geom_point(data=act.data, mapping=aes(x = Value, y = Coef, color = "red"))+
      theme_ridges() +
      #scale_x_continuous(breaks = seq(-0.2,0.2,0.01), lim = c(-0.2, 0.2))+
      theme(legend.position = "none")

    print(p)

    dev.off()



    png(fn2)

    p=ggplot(long.data, aes(x = Value, y = Coef,)) +
      geom_boxplot() +
      geom_point(data=act.data, mapping=aes(x = Value, y = Coef, color = "red"))+
      #scale_x_continuous(breaks = seq(-0.2,0.2,0.01), lim = c(-0.2, 0.2))+

      theme(legend.position = "none")

    print(p)

    dev.off()



  }

  sig.ind=cbind.data.frame(sig.snp, sig.cis)

  return(list(sig.indicator.mat=sig.ind, simu.coefs=simu.coefs, act.coefs=b))

}


# A helper function to get the adjacency matrix from a trio
get.adj=function(gmac.data=NULL, verbose=FALSE){

  adj=matrix(0, nrow=2, ncol=2)
  p.trans=summary(lm(trans.gene~., data = gmac.data))$coefficients[2:3,4]
  #print(summary(lm(trans.gene~., data = gmac.data)))
  p.cis=summary(lm(cis.gene~., data = gmac.data))$coefficients[2:3,4]
  if(verbose==TRUE){print(summary(lm(cis.gene~., data = gmac.data)));
    print(summary(lm(trans.gene~., data = gmac.data)))}

  p.ind.trans=ifelse(p.trans<0.01, 1, 0)
  #print(p.ind.trans)
  adj[2,2]=p.ind.trans[2]
  adj[1,2]=p.ind.trans[1]
  p.ind.cis=ifelse(p.cis<0.01, 1, 0)
  #print(p.ind.cis)
  adj[1,1]=p.ind.cis[1]

  return(list(adj=adj, p.mediation=p.trans[2]))

}


# A helper function that can classify an adjaceny matrix as one M0, M2, or M4
# (to be used on unmatched gmac trios therefore no need for M1 or M3)
class.adj=function(adj.mat=NULL){
  type=NA
  M4=matrix(c(1,1,0,1), nrow = 2, ncol = 2, byrow = T)
  M0.1=matrix(c(1,0,0,0), nrow = 2, ncol = 2, byrow = T)
  M0.2=matrix(c(0,1,0,0), nrow = 2, ncol = 2, byrow = T)
  M2.1=matrix(c(0,1,0,1), nrow = 2, ncol = 2, byrow = T)
  M2.2=matrix(c(1,0,0,1), nrow = 2, ncol = 2, byrow = T)
  M3=matrix(c(1,1,0,0), nrow = 2, ncol = 2, byrow = T)

  if(isTRUE(all.equal(M4, adj.mat, check.attributes=FALSE))){type="M4"}
  if(isTRUE(all.equal(M0.1, adj.mat, check.attributes=FALSE)) || isTRUE(all.equal(M0.2, adj.mat, check.attributes=FALSE))){type="M0"}
  if(isTRUE(all.equal(M2.1, adj.mat, check.attributes=FALSE)) || isTRUE(all.equal(M2.2, adj.mat, check.attributes=FALSE))){type="M2"}
  if(isTRUE(all.equal(M3, adj.mat, check.attributes=FALSE))){type="M3"}

  return(type)

}


# A function to reclassify all significant GMAC trios according to possible ADDIS model
# classifications M0:M4 based on the mediation regression
#
#                       Tj = b0 + b1Ci + b2Li + AXij + e
#
# Results are a table of GMAC model types and Addis model types side by side


reclass=function(tissue="WholeBlood", trio.ind.vec=NULL, mod.type="cis", verbose=FALSE){

  #read in data
  if(mod.type=="cis"){
    out.data=loadRData(fileName=paste0('/mnt/ceph/jarredk/GMACanalysis/', tissue, '/all_trios_output_cis.Rdata'))

  }else{

    out.data=loadRData(fileName=paste0('/mnt/ceph/jarredk/GMACanalysis/', tissue, '/all_trios_output_trans.Rdata'))

  }

  #allocate space:

  mod.type.gmac=NULL
  all.data=vector("list", length=length(trio.ind.vec))
  pvals=NULL

  #loop for all trios
  for(i in 1:length(trio.ind.vec)){

    trio.loc=as.matrix(out.data$input.list$trios.idx)[trio.ind.vec[i],]
    known.conf=t(out.data$input.list$known.conf)
    SNP=t(as.data.frame(out.data$input.list$snp.dat.cis)[trio.loc[1],])
    GE=t(as.data.frame(out.data$input.list$exp.dat)[trio.loc[-1],])
    which.covs=which(out.data$cov.indicator.list[trio.ind.vec[i],]==1)
    sel.cov.pool=as.data.frame(out.data$input.list$cov.pool)[which.covs,]

    if(verbose==TRUE){
      print("-----------GMAC-Selected-PCs--------------");print(row.names(sel.cov.pool))
    }

    sel.cov.pool=t(sel.cov.pool)

    data.mat=cbind.data.frame(SNP, GE, sel.cov.pool, known.conf)



    if(mod.type=="cis"){

      colnames(data.mat)=c("SNP","cis.gene", "trans.gene", colnames(sel.cov.pool), colnames(known.conf))
      #colnames(addis.data)=c("cis.gene", "trans.gene", "SNP", paste0("PC", sig.asso.pcs))

    }else{

      colnames(data.mat)=c("SNP", "trans.gene", "cis.gene", colnames(sel.cov.pool), colnames(known.conf))
      #colnames(addis.data)=c("trans.gene", "cis.gene", "SNP", addis.pcs)

    }

    #addis.data2=addis.data
    #convert.factors
    #data.mat$SNP=as.factor(data.mat$SNP)
    data.mat$pcr=as.factor(data.mat$pcr)
    data.mat$sex=as.factor(data.mat$sex)
    data.mat$platform=as.factor(data.mat$platform)
    #addis.data$SNP=as.factor(addis.data$SNP)

    ot=get.adj(data.mat, verbose=verbose)
    adj=ot$adj
    pvals[i]=ot$p.mediation
    if(verbose==TRUE){print(paste0("Trio_Num_",trio.ind.vec[i]));print(adj);print(pvals[i])}
    mod.type.gmac[i]=class.adj(adj)
    all.data[[i]]=data.mat

  }

  mod.inf=cbind.data.frame(mod.type.gmac, trio.ind.vec, pvals)
  colnames(mod.inf)=c("Inferred.Model", "Trio", "Pvalue")
  return(list(mod.inf=mod.inf, all.data=all.data))


}













#a utility function combining many of the above for ease of use


runit=function(indata=NULL, med.type=NULL,trio.number=NULL, mtype=""){

  #space allocation
  cors=as.data.frame(matrix(0, nrow=length(trio.number), ncol = 3))
  colnames(cors)=c("cor.SNP.cis", "cor.SNP.trans", "cor.cis.trans")
  partial.cors.addis=as.data.frame(matrix(0, nrow=length(trio.number), ncol = 3))
  colnames(partial.cors.addis)=c("ADDIS.pcor.SNP.cis", "ADDIS.pcor.SNP.trans", "ADDIS.pcor.cis.trans")
  partial.cors.gmac=as.data.frame(matrix(0, nrow=length(trio.number), ncol = 3))
  colnames(partial.cors.gmac)=c("GMAC.pcor.SNP.cis", "GMAC.pcor.SNP.trans", "GMAC.pcor.cis.trans")
  perm.p.addis=as.data.frame(matrix(0, nrow=length(trio.number), ncol = 2))
  colnames(perm.p.addis)=c("med.p.ADDIS", "med.stat.ADDIS")
  perm.p.gmac=as.data.frame(matrix(0, nrow=length(trio.number), ncol = 2))
  colnames(perm.p.gmac)=c("med.p.GMAC", "med.stat.GMAC")
  num.pcs=as.data.frame(matrix(0, nrow=length(trio.number), ncol = 2))
  colnames(num.pcs)=c("ADDIS.num.pcs", "GMAC.num.pcs")

  for( i in 1:length(trio.number)){
    #trio 133 ADDIS:M3
    types=reclass(tissue="WholeBlood",
                  trio.ind.vec=trio.number[i],
                  mod.type="cis",
                  verbose=FALSE)

    L2A=Lond2Addis.lookup(trio.index=trio.number[i], tissue.name="WholeBlood", with.pc=TRUE)

    cors[i,]=c(L2A$correlation[1, 2:3], L2A$correlation[2,3])

    if(length(colnames(L2A$correlation))>3){
      addis.pcs=colnames(L2A$correlation)[-c(1:3)]
      num.pcs[i,1]=length(addis.pcs)
    }else{
      addis.pcs=NULL
      num.pcs[i,1]=0
    }
    print(trio.number[i])
    print("--------------ADDIS-PCs--------------")
    print(addis.pcs)

    list.data=cross.regress(tissue="WholeBlood", trio.ind=trio.number[i], mod.type=med.type[i], addis.pcs=addis.pcs,
                            mtype = mtype, plot.it = TRUE)

    num.pcs[i,2]=length(colnames(list.data$GMAC)[-c(1:2,(length(colnames(list.data$GMAC))-3):length(colnames(list.data$GMAC)))])

    fname1=paste0("/mnt/ceph/jarredk/GMACanalysis/additional_plots/GMAC.nomatch.permutedREG.trio",
                  trio.number[i],".png")
    out.permreg=run.permuted.reg(trio=list.data$GMAC, nperms=10000, plot=TRUE, med.type = med.type[i], filename=fname1)
    perm.p.gmac[i,1]=out.permreg$p.value
    perm.p.gmac[i,2]=out.permreg$obs.wald.stat

    print("GMAC")
    print(out.permreg$p.value)
    print(out.permreg$obs.wald.stat)

    fname2=paste0("/mnt/ceph/jarredk/GMACanalysis/additional_plots/GMAC.nomatch.permutedREG.trio",
                trio.number[i],".addis.png")
    out.permreg=run.permuted.reg(trio=list.data$addis, nperms=10000, plot=TRUE, med.type = med.type[i], Alg="ADDIS", filename=fname2)
    perm.p.addis[i,1]=out.permreg$p.value
    perm.p.addis[i,2]=out.permreg$obs.wald.stat

    print("ADDIS")
    print(out.permreg$p.value)
    print(out.permreg$obs.wald.stat)

    r=check.ns.gmac(input.list.data=list.data,
                    iters=1000,
                    alpha=0.01,
                    print.it=FALSE,
                    which.mod="GMAC",
                    plot.it=TRUE,
                    save.name=paste0(trio.number[i], mtype))

    ra=check.ns.gmac(input.list.data=list.data,
                     iters=1000,
                     alpha=0.01,
                     print.it=FALSE,
                     which.mod="ADDIS",
                     plot.it=TRUE,
                     save.name=paste0(trio.number[i], mtype, "_A"))

    an.error.occured <- FALSE
    tryCatch( { pcor1 <- pcor(list.data$GMAC); print("------GMAC-pcors-----");print(pcor1$estimate[1:3,1:3]) }
              , error = function(e) {an.error.occured <<- TRUE})

    an.error.occured2 <- FALSE
    tryCatch( { pcor2 <- pcor(list.data$addis); print("------ADDIS-pcors-----");print(pcor2$estimate) }
              , error = function(e) {an.error.occured2 <<- TRUE})

    if(an.error.occured==TRUE){
      partial.cors.gmac[i,]=rep(NA, 3)
    }else{
      pcor1=pcor(list.data$GMAC)$estimate[1:3,1:3]
      partial.cors.gmac[i,]=c(pcor1[3,1:2], pcor1[1,2])
    }


    if(an.error.occured2==TRUE){
      partial.cors.addis[i,]=rep(NA,3)
    }else{
      pcor2=pcor(list.data$addis)$estimate[1:3,1:3]
      partial.cors.addis[i,]=c(pcor2[3,1:2], pcor2[1,2])
    }


    temp.tab=c(trio=trio.number[i], SNP=indata$snp[i], Cis.Gene=indata$cis[i], Trans.Gene=indata$trans[i],
               Class=indata$Addis.Class[i], Type=indata$Mediation.type[i], num.pcs[i,], perm.p.gmac[i,],
               perm.p.addis[i,], GMAC.pcor=partial.cors.gmac$GMAC.pcor.cis.trans[i],
               ADDIS.pcor=partial.cors.addis$ADDIS.pcor.cis.trans[i])

    if(i==1){
      write.table(temp.tab,
                  file = paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/trios.pcors.summary.table",mtype,".txt"),
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = TRUE,
                  sep = "\t")
    }else{

      write.table(temp.tab,
                  file = paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/trios.pcors.summary.table",mtype,".txt"),
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = FALSE,
                  append = TRUE,
                  sep = "\t")
    }
  }

  final.table=cbind.data.frame(trio=trio.number, SNP=indata$snp, Cis.Gene=indata$cis, Trans.Gene=indata$trans,
                               Class=indata$Addis.Class, Type=indata$Mediation.type, num.pcs, cors, perm.p.gmac,
                               perm.p.addis, partial.cors.gmac, partial.cors.addis)

  save(final.table, file = paste0("/mnt/ceph/jarredk/GMACanalysis/master_tables/trios.pcors.summary.table",mtype,".Rdata"))

  return(list(datalist=list.data, table=final.table))

}




#a simple function to identify possible suppressor variables in the set of confounding variables
#included in a trio

id.suppress=function(data=NULL, verbose=FALSE, reg.compare=FALSE){


  supp.index=NULL
  names1=NULL
  #data$pcr=as.factor(data$pcr)
  #data$platform=as.factor(data$platform)
  #data$sex=as.factor(data$platform)

  for(i in 3:(dim(data)[2]-3)){

    modelx1=lm(trans.gene~., data=data[, c(1:3)])
    modelx12=lm(trans.gene~., data=data[,c(1:3,i)])
    modelx2=lm(trans.gene~., data=data[,c(2,i)])

    if(verbose==TRUE){

      print(paste0("========================================Models_", colnames(data)[i],"====================================="))
      print(summary(modelx1))
      print(summary(modelx2))
      print(summary(modelx12))

    }

    regx1=anova(modelx1)
    regx12=anova(modelx12)
    regx2=anova(modelx2)

    ssr.x12=sum(regx12$`Sum Sq`[-length(regx12$`Sum Sq`)])
    ssr.x1=sum(regx1$`Sum Sq`[-length(regx1$`Sum Sq`)])
    ssr.x2.x1=ssr.x12-ssr.x1
    ssr.x2=sum(regx2$`Sum Sq`[-length(regx2$`Sum Sq`)])

    if(verbose==TRUE){

      print("========================================SSRs=====================================")
      print(paste0("SSR(Cis.Gene,",colnames(data)[i],")=",ssr.x12))
      print(paste0("SSR(Cis.Gene)=",ssr.x1))
      print(paste0("SSR(",colnames(data)[i],"| Cis.Gene)=",ssr.x2.x1))
      print(paste0("SSR(",colnames(data)[i],")=",ssr.x2))

    }

    supp.index[i-2]=ssr.x2.x1/ssr.x2

    names1[i-2]=colnames(data)[i]


  }

  names(supp.index)=names1
  suppressors=names(supp.index[supp.index>1])

  if(reg.compare==TRUE){

    for(j in 1:2){

      if(j==1){
        print("=============Regression With Only Suppressors=============")
        supp.data=data[,match(suppressors, colnames(data))]
      }else{
        print("=============Regression Without Suppressors===============")
        supp.data=data[,-match(suppressors, colnames(data))]
      }

      supp.data$trans.gene=data$trans.gene
      supp.data$cis.gene=data$cis.gene
      supp.data$SNP=data$SNP
      supp.data$pcr=as.factor(data$pcr)
      supp.data$platform=as.factor(data$platform)
      supp.data$sex=as.factor(data$sex)

      print(summary(lm(trans.gene~., data = supp.data)))
    }

  }


  return(supp.index)


}
























#===========================================================================
#----------------------------Simulations------------------------------------
#===========================================================================


#A function which simulates the trans gene of a given trio by
#                 Tj=a + b1 Ci + b2 Li + V + errors

#where V represents the PCs significant at level alpha from the original regression of the trio using all
#GMAC adaptively selected PCs


simu1=function(data=NULL,alpha=0.001, mod.type=NULL, verbose=TRUE){



  data$pcr=as.factor(data$pcr)
  data$platform=as.factor(data$platform)
  data$sex=as.factor(data$sex)



  if(mod.type=="Both"){
    mod=lm(trans.gene~., data=data)
  }else if(mod.type=="Cis.Med"){
    mod=lm(trans.gene~., data=data)
  }else{
    mod=lm(cis.gene~., data=data)
  }

  if(verbose==TRUE){
    print("========================Results-for-original-Tj=======================")
    print(summary(mod))
  }

  X=model.matrix(mod)
  #print(colnames(X))
  b=as.data.frame(summary(mod)$coefficients)$Estimate
  names(b)=row.names(as.data.frame(summary(mod)$coefficients))
  idx=which(as.data.frame(summary(mod)$coefficients)$`Pr(>|t|)`[-c(1:3)]<alpha)+3
  keep.pcs=row.names(as.data.frame(summary(mod)$coefficients))[idx]
  sigma=summary(mod)$sigma

  #S=id.suppress(data, verbose=FALSE)
  #p=length(na.omit(match(keep.pcs, names(S[which(S>1)]))))/length(keep.pcs)

  # if(verbose==TRUE){
  #   print('Suppressors:')
  #   print(S[which(S>1)])
  #   print(paste("Sig.PCs:", paste(keep.pcs, collapse = "," )))
  #   print(paste0("Proportion of Sig.PCs That Are Suppressors = ", p))
  # }

  b.gen=b[c(1:3, idx)]
  #print(length(b.gen))
  X.gen=X[,c(1:3, match(keep.pcs, colnames(X)[-1])+1)]
  #print(dim(X.gen))
  errors=rnorm(dim(X.gen)[1], mean = 0, sd = sigma)

  Tj=X.gen%*%b.gen+errors

  data.new=data

  #print(head(data.new))

  if(mod.type=="Both"){

    data.new$trans.gene=Tj
    mod2=lm(trans.gene~., data = data.new)

  }else if(mod.type=="Cis.Med"){

    data.new$trans.gene=Tj
    mod2=lm(trans.gene~., data = data.new)

  }else{

    data.new$cis.gene=Tj
    mod2=lm(cis.gene~., data = data.new)
  }


  if(verbose==TRUE){
    print("==============coefficients-used-to-generate-Tj========================")
    print(b.gen)
    print("=====================model-matrix-for-Tj==============================")
    print(head(X.gen))
    print(dim(X.gen))
    print("========================Results-for-simulated-Tj======================")
    print(summary(mod2))
    print(dim(model.matrix(mod2)))
  }


  sim.coef.mat=as.data.frame(summary(mod2)$coefficients)


  if(mod.type=="Both"){
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="cis.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="cis.gene")]
  }else if(mod.type=="Cis.Med"){
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="cis.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="cis.gene")]
  }else{
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="trans.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="trans.gene")]
  }

  #print(cor(data.new[, (length(colnames(data.new))-2):length(colnames(data.new))])[,1:3])
  return(list(sim.data=data.new,
              pvalue=sim.p,
              b.coef=sim.b,
              dim.diff=(dim(model.matrix(mod2))[2]-1)-(dim(X.gen)[2]-1) ))

}















#=========================================================================
#A function which simulates the trans gene of a given trio by
#                 Tj=a + b1 Ci + b2 Li + V + errors

#where V represents the PC's not selected by addis or GMAC. The goal is to better understand the stability of the regression


simu2=function(tissue = "WholeBlood", data=NULL, seed=NULL, mod.type=NULL, n=10, verbose=TRUE){

  if(isFALSE(is.null(seed))){
    set.seed(seed = seed)
  }

  out.vec=rep(0, 6)
  names(out.vec)=c("Num.PCs.GMAC", "Per.Var", "Med.pvalue.MRPC", "Med.coef.MRPC", "Med.pvalue.GMAC", "Med.coef.GMAC")



  data$pcr=as.factor(data$pcr)
  data$platform=as.factor(data$platform)
  data$sex=as.factor(data$sex)

  file2=paste("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/",
              tissue,"_AllPC/PCs.matrix.",
              tissue, ".RData", sep="")

  PCs=loadRData(fileName=file2)
  match.pcs=na.omit(match(colnames(data), colnames(PCs)))

  out.vec[1]=length(match.pcs)

  PCs2=PCs[,-match.pcs]

  pcs.chosen=PCs2[,sample(c(1:dim(PCs2)[2]), n)]

  new.data=cbind.data.frame(data, pcs.chosen)

  #print(head(new.data))


  if(mod.type=="Both"){
    true.model=lm(trans.gene~., data=new.data)
  }else if(mod.type=="Cis.Med"){
    true.model=lm(trans.gene~., data=new.data)
  }else{
    true.model=lm(cis.gene~., data=new.data)
  }


  if(verbose==TRUE){
    print("========================Results-for-generating-model-of-Tj===========================")
    print(summary(true.model))
  }

  X=model.matrix(true.model)
  b=as.data.frame(summary(true.model)$coefficients)$Estimate
  names(b)=row.names(as.data.frame(summary(true.model)$coefficients))
  sigma=summary(true.model)$sigma

  if(verbose==TRUE){
    print("==============coefficients-used-to-generate-Tj========================")
    print(b)
    print("=====================model-matrix-for-Tj==============================")
    print(head(X))
  }
  errors=rnorm(dim(X)[1],mean = 0, sd = sigma)
  #errors=rnorm(dim(X)[1])

  Tj=X%*%b+errors

  data.new=data

  if(mod.type=="Both"){

    data.new$trans.gene=Tj
    model2=lm(trans.gene~., data = data.new)

  }else if(mod.type=="Cis.Med"){

    data.new$trans.gene=Tj
    model2=lm(trans.gene~., data = data.new)
  }else{

    data.new$cis.gene=Tj
    model2=lm(cis.gene~., data = data.new)
  }



  if(verbose==TRUE){
    print("========================Results-for-simulated-Tj======================")
    print(summary(model2))
  }


  sim.coef.mat=as.data.frame(summary(model2)$coefficients)

  if(mod.type=="Both"){
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="cis.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="cis.gene")]
  }else if(mod.type=="Cis.Med"){
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="cis.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="cis.gene")]
  }else{
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="trans.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="trans.gene")]
  }




  return(list(data.sim=data.new,
              pvalue=sim.p,
              b.coef=sim.b,
              dim.diff=n))
}




#Simulation 3 TGM
#A function to simulate the trans gene under the GMAC parametric model and then
#applying the MRPC-ADDIS parametric model under the simulated trans gene to determine
#The power when only MRPC selected confounders are included - to be wrapped by simu12


simu3=function(MRPC.data=NULL, GMAC.data=NULL, mod.type=NULL, verbose=NULL){

  GMAC.data$pcr=as.factor(GMAC.data$pcr)
  GMAC.data$platform=as.factor(GMAC.data$platform)
  GMAC.data$sex=as.factor(GMAC.data$sex)


  if(mod.type=="Both"){
    true.model=lm(trans.gene~., data=GMAC.data)
  }else if(mod.type=="Cis.Med"){
    true.model=lm(trans.gene~., data=GMAC.data)
  }else{
    true.model=lm(cis.gene~., data=GMAC.data)
  }

  if(verbose==TRUE){
    print("========================Results-for-True-Tj======================")
    print(summary(true.model))
  }

  b=as.data.frame(summary(true.model)$coefficients)$Estimate
  X=model.matrix(true.model)
  sigma=summary(true.model)$sigma

  Tj=X%*%b+rnorm(dim(X)[1], mean = 0, sd=sigma)

  new.data=MRPC.data
  #new.data$trans.gene=Tj


  if(mod.type=="Both"){
    new.data$trans.gene=Tj
    analysis.model=lm(trans.gene~., data = new.data)
  }else if(mod.type=="Cis.Med"){
    new.data$trans.gene=Tj
    analysis.model=lm(trans.gene~., data = new.data)
  }else{
    new.data$cis.gene=Tj
    analysis.model=lm(cis.gene~., data = new.data)
  }


  if(verbose==TRUE){
    print("========================Results-for-simulated-Tj======================")
    print(summary(analysis.model))
  }


  sim.coef.mat=as.data.frame(summary(analysis.model)$coefficients)

  if(mod.type=="Both"){
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="cis.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="cis.gene")]
  }else if(mod.type=="Cis.Med"){
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="cis.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="cis.gene")]
  }else{
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="trans.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="trans.gene")]
  }


  return(list(sim.data=new.data, pvalue=sim.p, b.coef=sim.b))


}






#------------------------------------------

#simulation 4 GSS


simu4=function(GMAC.data=NULL, mod.type=NULL, verbose=NULL){

  GMAC.data$pcr=as.factor(GMAC.data$pcr)
  GMAC.data$platform=as.factor(GMAC.data$platform)
  GMAC.data$sex=as.factor(GMAC.data$sex)


  if(mod.type=="Both"){
    true.model=lm(trans.gene~., data=GMAC.data)
  }else if(mod.type=="Cis.Med"){
    true.model=lm(trans.gene~., data=GMAC.data)
  }else{
    true.model=lm(cis.gene~., data=GMAC.data)
  }

  if(verbose==TRUE){
    print("========================Results-for-True-Tj======================")
    print(summary(true.model)$sigma)
  }

  #simulate trans gene
  b=as.data.frame(summary(true.model)$coefficients)$Estimate
  #print(head(b))
  X=model.matrix(true.model)
  #print(head(X))
  sigma=summary(true.model)$sigma
  #sigma=0.0117
  #print(sigma)

  Tj=X%*%b+rnorm(dim(X)[1], mean = 0, sd=sigma)

  #print("simulated gene")
  #print(head(Tj))
  new.data=GMAC.data
  #new.data$trans.gene=Tj



  if(mod.type=="Both"){
    new.data$trans.gene=Tj
    analysis.model=lm(trans.gene~., data = new.data)
  }else if(mod.type=="Cis.Med"){
    new.data$trans.gene=Tj
    analysis.model=lm(trans.gene~., data = new.data)
  }else{
    new.data$cis.gene=Tj
    analysis.model=lm(cis.gene~., data = new.data)
  }


  if(verbose==TRUE){
    print("========================Results-for-simulated-Tj======================")
    print(summary(analysis.model))
  }


  sim.coef.mat=as.data.frame(summary(analysis.model)$coefficients)

  if(mod.type=="Both"){
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="cis.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="cis.gene")]
  }else if(mod.type=="Cis.Med"){
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="cis.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="cis.gene")]
  }else{
    sim.p=sim.coef.mat$`Pr(>|t|)`[which(row.names(sim.coef.mat)=="trans.gene")]
    sim.b=sim.coef.mat$Estimate[which(row.names(sim.coef.mat)=="trans.gene")]
  }


  return(list(sim.data=new.data, pvalue=sim.p, b.coef=sim.b))


}







#--------------------------------------------------Simulation Wrapper-------------------------------------------------------

#A function which wraps simu1, simu2 and cross.regress together

run.simu12=function(tissue="WholeBlood", trios=NULL, mod.type.vec=NULL, l1.table=NULL, alpha=0.001, n=10, seed=222, num.perms=1000,
                    verbose=FALSE){
  #
  # sim.sets=vector("list", length = length(trios))
  # orig.sets=vector("list", length = length(trios))
  # names(sim.sets)=paste0("trio",trios)
  # names(orig.sets)=paste0("trio",trios)
  #
  #preallocate table
  n1=c("Trio.Num", "Num.pcs.MRPC", "Per.var.MRPC", "Num.pcs.GMAC", "Per.var.GMAC")
  out.mat=as.data.frame(matrix(0, nrow = length(trios), ncol=length(n1)))
  colnames(out.mat)=n1

  #read in badsha peerdata.v8 file
  Peerdata.V8=loadRData(paste0("/mnt/lfs2/mdbadsha/peer_example/SNP_cis_trans_files/GTEx_version_8/PEER_Files_V8/Peerdata.",tissue,".V8.RData"))

  #preform PCA in same method as badsha and retain eigen values
  Peerdata2 <- t(Peerdata.V8)
  Peerdata3 <- Peerdata2[,apply(Peerdata2, 2, var, na.rm=TRUE) != 0]
  # pca matrix
  PCs<- prcomp(Peerdata3,scale=TRUE)

  evs=PCs$sdev^2
  names(evs)=colnames(PCs$x)


  for(i in 1:length(trios)){
    if(verbose==TRUE){
      print("***************************************************************************")
      print(paste("*************************", "Trios Number ", trios[i],"******************************"))
      print("***************************************************************************")
    }

    L2A=Lond2Addis.lookup(trio.index=trios[i], tissue.name=tissue, with.pc=TRUE)

    if(length(colnames(L2A$correlation))>3){
      addis.pcs=colnames(L2A$correlation)[-c(1:3)]
    }else{
      addis.pcs=NULL
    }


    if(mod.type.vec[i]=="Both"){

      list.data=cross.regress(tissue=tissue,
                              trio.ind=trios[i],
                              mod.type="cis",
                              addis.pcs=addis.pcs,
                              verbose = FALSE)

    }else if(mod.type.vec[i]=="Cis.Med"){

      list.data=cross.regress(tissue=tissue,
                              trio.ind=trios[i],
                              mod.type="cis",
                              addis.pcs=addis.pcs,
                              verbose = FALSE)
    }else{

      list.data=cross.regress(tissue=tissue,
                              trio.ind=trios[i],
                              mod.type="trans",
                              addis.pcs=addis.pcs,
                              verbose = FALSE)

    }


    total.var=sum(evs)
    addis.var.act=evs[match(addis.pcs, names(evs))]
    gmac.var.act=evs[na.omit(match(colnames(list.data$GMAC), names(evs)))]

    if(verbose==TRUE){
      print("-------------------variance-accounting---------------------")
      print(paste0("GMAC = ", round(sum(gmac.var.act)/total.var, 6)*100, "%"))
      print(paste0("ADDIS = ", round(sum(addis.var.act)/total.var, 6)*100, "%"))
    }


    #orig.sets[[i]]=list.data$GMAC
    #STM simulation
    print("running Simulation 1...")
    out=simu1(list.data$GMAC, alpha = alpha, mod.type=mod.type.vec[i], verbose=verbose)
    print("...done")

    #sim.sets[[i]]=out

    #LTM simulation
    print("running Simulation 2...")
    if(n=="random"){
      nn=floor(runif(1, 1, 20))
      print(nn)
      out2=simu2(tissue = tissue, data=list.data$GMAC, mod.type=mod.type.vec[i], seed=seed, n=nn, verbose=verbose)
    }else{
      out2=simu2(tissue = tissue, data=list.data$GMAC, mod.type=mod.type.vec[i], seed=seed, n=n, verbose=verbose)
    }

    #run TGM simulation
    print("running Simulation 3...")
    out3=simu3(MRPC.data=list.data$addis, GMAC.data=list.data$GMAC, mod.type=mod.type.vec[i], verbose=F)
    #run GSS simulation
    print("running Simulation 4...")
    out4=simu4(GMAC.data=list.data$GMAC, mod.type=mod.type.vec[i], verbose=T)

    print("...done")

    print("############################################################################")


    match.pcs=na.omit(match(colnames(list.data$GMAC), colnames(PCs$x)))

    #run the permuted regressions

    print("Running Permutation Regressions...")
    #get MRPC nominal and perm p
    perm1=run.permuted.reg(list.data$addis, nperms=num.perms, Alg="ADDIS", med.type = mod.type.vec[i])
    #STM, LTM, and TGM simulations
    perm2=run.permuted.reg(out$sim.data, nperms=num.perms, Alg="GMAC", med.type = mod.type.vec[i])
    perm3=run.permuted.reg(out2$data.sim, nperms=num.perms, Alg="GMAC", med.type = mod.type.vec[i])
    perm4=run.permuted.reg(out3$sim.data, nperms=num.perms, Alg="ADDIS", med.type = mod.type.vec[i])
    #Get GMAC nominal and perm p
    perm5=run.permuted.reg(list.data$GMAC, nperms=num.perms, Alg="GMAC", med.type = mod.type.vec[i])
    #GSS simulation
    perm6=run.permuted.reg(out4$sim.data, nperms=num.perms, Alg="GMAC", med.type = mod.type.vec[i])
    print("-----------checking values to ensure accuracy--------")
    print(perm6$nominal.p)
    print(perm6$p.value)
    print(out4$pvalue)
    print(out4$b.coef)
    print(list.data$GMAC[1:5,1:5])
    print(out4$sim.data[1:5,1:5])

    print("...done")

    #get the permutation and nominal pvalues
    out.mat$Perm.p.MRPC[i]=perm1$p.value
    out.mat$Nominal.p.MRPC[i]=perm1$nominal.p
    #STM
    out.mat$STM.perm.p[i]=perm2$p.value
    out.mat$Nominal.p.STM[i]=perm2$nominal.p
    #LTM
    out.mat$LTM.perm.p[i]=perm3$p.value
    out.mat$Nominal.p.LTM[i]=perm3$nominal.p
    #TGM
    out.mat$TGM.perm.p[i]=perm4$p.value
    out.mat$Nominal.p.TGM[i]=perm4$nominal.p
    #GMAC
    out.mat$Perm.p.GMAC[i]=perm5$p.value
    #GSS
    out.mat$GSS.perm.p[i]=perm6$p.value
    out.mat$Nominal.p.GSS[i]=perm6$nominal.p

    #store
    out.mat$Num.pcs.MRPC[i]=length(addis.pcs)
    out.mat$Per.var.MRPC[i]=round(sum(addis.var.act)/total.var, 6)*100
    out.mat$Num.pcs.GMAC[i]=length(match.pcs)
    out.mat$Per.var.GMAC[i]=round(sum(gmac.var.act)/total.var, 6)*100
    #regression values
    out.mat$Med.pvalue.GMAC[i]=list.data$Regress$p.gmac
    out.mat$Med.pvalue.MRPC[i]=list.data$Regress$p.addis
    out.mat$Med.coef.GMAC[i]=list.data$Regress$b.addis.cis.gene
    out.mat$Med.coef.MRPC[i]=list.data$Regress$b.gmac.cis.gene
    #Small truth simulation
    out.mat$STM.med.p[i]=out$pvalue
    out.mat$STM.med.coef[i]=out$b.coef
    out.mat$STM.dim.diff[i]=out$dim.diff
    #Large truth simulation
    out.mat$LTM.med.p[i]=out2$pvalue
    out.mat$LTM.med.coef[i]=out2$b.coef
    out.mat$LTM.dim.diff[i]=out2$dim.diff
    #True GMAC simulation
    out.mat$TGM.med.p[i]=out3$pvalue
    out.mat$TGM.med.coef[i]=out3$b.coef
    #out.mat$TGM.dim.diff[i]=out3$dim.diff
    #GMAC Self Simulation
    out.mat$GSS.med.p[i]=out4$pvalue
    out.mat$GSS.med.coef[i]=out4$b.coef

  }

  #meta info
  out.mat$Trio.Num=trios
  out.mat$Med.type=mod.type.vec
  out.mat$ADDIS.inf.Class=l1.table$Addis.Class[match(out.mat$Trio.Num,  l1.table$Trio.Num)]
  out.mat$GMAC.Nominal.p=l1.table$Perm.reg.p[match(out.mat$Trio.Num,  l1.table$Trio.Num)]
  out.mat$cis.gene=l1.table$cis[match(out.mat$Trio.Num,  l1.table$Trio.Num)]
  out.mat$trans.gene=l1.table$trans[match(out.mat$Trio.Num,  l1.table$Trio.Num)]
  out.mat$snp=l1.table$snp[match(out.mat$Trio.Num,  l1.table$Trio.Num)]
  #print(out.mat$Perm.p)

  #return(list(data.sim=sim.sets, data.orig=orig.sets, out.mat=out.mat))
  return(out.mat=out.mat)
}



















