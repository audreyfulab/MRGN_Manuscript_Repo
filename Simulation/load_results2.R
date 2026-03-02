# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

library(MRGN)
library(qvalue)
source("adapted_GMAC_func/GMACpostproc.R")

params = loadRData(file="Simulation/data/mrpc_v_mrgn_v_gmac_5k_params_all_mods_conf_types.RData")
gmac.cis = loadRData(file="Simulation/data/gmac_5k_cis_results_all_confs_all_mods.RData")
gmac.trans = loadRData(file="Simulation/data/gmac_5k_trans_results_all_confs_all_mods.RData")
mrgn.inf = loadRData("Simulation/data/int_and_child_filtered_data/mrgn_5k_inf_results_all_confs_filtered_all_mods.RData")
mrgn.confs = loadRData(file = "Simulation/data/int_and_child_filtered_data/mrgn_5k_conf_list_all_confs_filtered_all_mods.RData")
mrgn.inf.alpha01 = loadRData(file = "Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_inference_results.RData")
mrgn.confs.alpha01 = loadRData(file = "Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_conf_list.RData")
mrgn.reg.res = loadRData(file="Simulation/data/int_and_child_filtered_data/mrgn_5k_regres_results_all_confs_filtered_all_mods.RData")
mrgn.reg.res.alpha01 = loadRData(file = "Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_regression_results.RData")
mrgn.resmat.alpha01 = do.call('rbind', mrgn.reg.res.alpha01)
mrpc.inf = loadRData(file="Simulation/data/int_and_child_filtered_data/mrpc_5k_small_results_all_confs_filtered_all_mods.RData")
mrpc.confs = loadRData(file = "Simulation/data/int_and_child_filtered_data/mrpc_5k_conf_list_all_confs_filtered_all_mods.RData")
all.mods = loadRData(file="Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types.RData")


colnames(gmac.trans$output.table) = c(paste0("pval_", colnames(gmac.trans$output.table)[1:2]),
                                      paste0("effect_change_", colnames(gmac.trans$output.table)[1:2]))

mrgn.reg.res = list2DF(mrgn.reg.res)
mrpc.inf2=unlist(lapply(mrpc.inf, function(x) x$model))
all.mods = lapply(all.mods, function(x) x$data)

xx=matrix(make.unique(unlist(lapply(all.mods, function(x) colnames(x)[1:3]))), nrow=length(all.mods), ncol = 3, byrow=T)
colnames(xx)=c("V", "Ti", "Tj")
xx=cbind.data.frame(xx, trio.idx = c(1:length(all.mods)))

#=================================================
replace.tf=function(tf.cols=NULL){

  type.of.med=NULL
  for(i in 1:dim(tf.cols)[1]){

    #if both true = bi-directed edge
    if(isTRUE(all.equal(tf.cols[i,], c(TRUE,TRUE), check.attributes=FALSE))){
      type.of.med[i]="Both"
      #if first true = cis mediation
    }else if(isTRUE(all.equal(tf.cols[i,], c(TRUE,FALSE), check.attributes=FALSE))){
      type.of.med[i]="Cis"
      #if second true = trans mediation
    }else if(isTRUE(all.equal(tf.cols[i,], c(FALSE,TRUE), check.attributes=FALSE))){
      type.of.med[i]="Trans"
      #if no true = no mediation
    }else if(isTRUE(all.equal(tf.cols[i,], c(FALSE,FALSE), check.attributes=FALSE))){
      type.of.med[i]="No.Med"
    }
  }

  #return
  return(type.of.med)


}
#=================================================


#pp.cis = run.postproc(gmac.cis$output.table$pval_Known_sel_pool, trio.ref = xx)
#pp.trans = run.postproc(gmac.trans$output.table$pval_Known_sel_pool, trio.ref = xx)

# apply qvalue correction
qval.cis.obj = qvalue(p = gmac.cis$output.table$pval_Known_sel_pool,
                       fdr.level = 0.1,
                       lambda = seq(0.05, max(gmac.cis$output.table$pval_Known_sel_pool), 0.05))

qval.trans.obj = qvalue(p = gmac.trans$output.table$pval_Known_sel_pool,
                         fdr.level = 0.1,
                         lambda = seq(0.05, max(gmac.trans$output.table$pval_Known_sel_pool), 0.05))

# determine which are significant at the 0.01 level (similar to MRGN)
# cis.sig = gmac.cis$output.table$pval_Known_sel_pool<0.01
# trans.sig = gmac.trans$output.table$pval_Known_sel_pool<0.01

gmac.cis$output.table$Sig = qval.cis.obj$significant
gmac.cis$output.table$qvalue = qval.cis.obj$qvalue
gmac.cis$output.table$at_05_cutoff = gmac.cis$output.table$pval_Known_sel_pool<0.05
gmac.cis$output.table$at_01_cutoff = gmac.cis$output.table$pval_Known_sel_pool<0.01
colnames(gmac.cis$output.table)= paste0("Cis_", colnames(gmac.cis$output.table))

gmac.trans$output.table$Sig = qval.trans.obj$significant
gmac.trans$output.table$qvalue = qval.trans.obj$qvalue
gmac.trans$output.table$at_05_cutoff = gmac.trans$output.table$pval_Known_sel_pool<0.05
gmac.trans$output.table$at_01_cutoff = gmac.trans$output.table$pval_Known_sel_pool<0.01
colnames(gmac.trans$output.table)= paste0("Trans_", colnames(gmac.trans$output.table))

gmac.tf=cbind.data.frame(cis = gmac.cis$output.table$Cis_Sig, trans = gmac.trans$output.table$Trans_Sig)
gmac.tf.at.05=cbind.data.frame(cis = gmac.cis$output.table$Cis_at_05_cutoff, trans = gmac.trans$output.table$Trans_at_05_cutoff)
gmac.tf.at.01=cbind.data.frame(cis = gmac.cis$output.table$Cis_at_01_cutoff, trans = gmac.trans$output.table$Trans_at_01_cutoff)

gmac.inf = replace.tf(tf.cols = gmac.tf)
gmac.inf.at.05 = replace.tf(tf.cols = gmac.tf.at.05)
gmac.inf.at.01= replace.tf(tf.cols = gmac.tf.at.01)


mrgn.num.confs = unlist(lapply(mrgn.confs, function(x) if(is.null(ncol(x))) 0 else ncol(as.matrix(x))))
mrgn.num.confs.alpha01 = unlist(lapply(mrgn.confs.alpha01, function(x) if(is.null(ncol(x))) 0 else ncol(as.matrix(x))))
mrpc.num.confs = unlist(lapply(mrpc.confs, function(x) if(is.null(ncol(x))) 0 else ncol(as.matrix(x))))
gmac.num.confs =  apply(gmac.cis$cov.indicator.list, 1, sum)

master.table=cbind.data.frame(params,
                              gmac.cis$output.table,
                              gmac.trans$output.table,
                              MRGN.pval.b21 = unlist(mrgn.reg.res[8,]),
                              MRGN.pval.b21.CSnoFDR.alpha01 = mrgn.resmat.alpha01[,8],
                              MRGN.pval.b22 = unlist(mrgn.reg.res[10,]),
                              MRGN.pval.b22.CSnoFDR.alpha01 = mrgn.resmat.alpha01[,10],
                              MRPC.inference.CSplusFDR = mrpc.inf2,
                              MRGN.inference.CSplusFDR = mrgn.inf,
                              MRGN.inference.CSnoFDR.alpha01 = mrgn.inf.alpha01,
                              GMAC.inference.with.qval = gmac.inf,
                              GMAC.inference.at.05.cutoff = gmac.inf.at.05,
                              GMAC.inference.at.01.cutoff = gmac.inf.at.01,
                              Number.Confs.Selected.MRPC = mrpc.num.confs,
                              Number.Confs.Selected.MRGN.CSplusFDR = mrgn.num.confs,
                              Number.Confs.Selected.MRGN.CSnoFDR.alpha01 = mrgn.num.confs.alpha01,
                              Number.Confs.Selected.GMAC.CSplusFDR = gmac.num.confs)

save(master.table, file="Simulation/data/master_table_all_confs_all_mods.RData")
save(gmac.trans, file="Simulation/data/gmac_5k_trans_results_all_all_mods_conf_types_preproc.RData")
save(gmac.cis, file="Simulation/data/gmac_5k_cis_results_all_mods_all_conf_types_preproc.RData")
