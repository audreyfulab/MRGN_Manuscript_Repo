
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

#loads all relevant data and merges 15 and > 15 conf results

library(MRGN)

#read in helper functions
source('Manuscript/scripts/MRGN_write_up_helper_functions.R')
source('Manuscript/other/helpers.R')
path_supptabs = "Manuscript/supplementary_tables/"
path_suppfigs = "Manuscript/supplementary_figures/"
path_figs = "Manuscript/figures/"
path_tabs = "Manuscript/tables/"
path_tablescraps = "Manuscript/other/tablescraps/"

#read in simulation data
params=loadRData(file="Simulation/data/mrpc_v_mrgn_v_gmac_5k_params_all_mods_conf_types.RData")

master.table = loadRData(file="Simulation/data/master_table_all_confs_all_mods.RData")


gmac.trans = loadRData(file="Simulation/data/gmac_5k_trans_results_all_all_mods_conf_types_preproc.RData")

gmac.cis=loadRData(file="Simulation/data/gmac_5k_cis_results_all_mods_all_conf_types_proproc.RData")

mrgn.inf=loadRData("Simulation/data/int_and_child_filtered_data/mrgn_5k_inf_results_all_confs_filtered_all_mods.RData")

#mrgn.reg.res=loadRData(file="Simulation/data/int_and_child_filtered_data/mrgn_5k_regres_results_all_confs_filtered_all_mods.RData")

mrpc.inf=loadRData(file="Simulation/data/int_and_child_filtered_data/mrpc_5k_small_results_all_confs_filtered_all_mods.RData")

datasets=loadRData(file="Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types.RData")

conf.list.mrgn=loadRData(file="Simulation/data/int_and_child_filtered_data/mrgn_5k_conf_list_all_confs_filtered_all_mods.RData")

#conf.list.mrpc = loadRData(file="Simulation/data/int_and_child_filtered_data/mrpc_5k_conf_list_all_confs_filtered_all_mods.RData")

conf.mat = loadRData(file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types_conf_mat.RData")
#======================================================================================
#read simulation results for conf selection with no FDR (selected at alpha < 0.01)

conf.list.alpha01 = loadRData('Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_conf_list.RData')
mrgn.inf.alpha01 = loadRData('Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_inference_results.RData')



#======================================================================================
#comp time results
#comp time table::
mrpc.comp.times = loadRData(file = "Simulation/data/diagnostics/mrpc_comp_times_with_filtering.RData")
mrgn.comp.times = loadRData(file = "Simulation/data/diagnostics/mrgn_comp_times_with_filtering.RData")


many.conf.params=loadRData(file="Simulation/data/many_conf_data/mrpc_v_mrgn_v_gmac_300_params_all_mods_conf_types.RData")

many.conf.data = loadRData(file =
                             "Simulation/data/many_conf_data/mrgn_v_gmac_v_mrpc_300_datasets_all_mods_conf_types.RData")

many.conf.data2 = lapply(many.conf.data, function(x) x$data)

many.conf.conf.mat = loadRData("Simulation/data/many_conf_data/mrgn_v_gmac_v_mrpc_300_datasets_all_mods_conf_types_conf_mat.RData")

many.conf.reg.res = loadRData("Simulation/data/many_conf_data/mrgn_300_regres_results_all_confs_all_mods.RData")

many.conf.mrgn.conf.list = loadRData("Simulation/data/many_conf_data/mrgn_300_conf_list_all_confs_all_mods_REG.RData")

mrgn.many.conf.inf = loadRData(file= "Simulation/data/many_conf_data/mrgn_300_inf_results_all_confs_all_mods.RData")

master.table.many.conf = loadRData(file="Simulation/data/many_conf_data/master_table_MANY_CONF_all_confs_all_mods.RData")

params.with.mrgn.comp.time = loadRData("Simulation/data/diagnostics/mrgn_many_conf_with_comp_times.RData")

gmac.trans.many.conf = loadRData("Simulation/data/many_conf_data/gmac_300_trans_results_all_all_mods_conf_types_preproc.RData")

gmac.cis.many.conf = loadRData("Simulation/data/many_conf_data/gmac_300_cis_results_all_mods_all_conf_types_proproc.RData")

gmac.inf.times = loadRData(file = 'Manuscript/other/results_all_times.RData')
#==================================================================
#read in many conf simulation results using MRGN with Confounder Selection and No FDR (alpha <0.01)
mrgn.mc.inf.alpha01 = loadRData(file = "Simulation/data/many_conf_data/alpha_01_selected_confs_results/mrgn_many_conf_liberal_alpha01_inf_results.RData")
mrgn.mc.confs.alpha01 = loadRData(file = "Simulation/data/many_conf_data/alpha_01_selected_confs_results/MRGN_many_conf_liberal_sel_alpha01_conf_list.RData")


#==================================================================
#read in simulation results using MRGN with Confounder Selection and No FDR (alpha <0.01)
conf.list.alpha01 = loadRData('Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_conf_list.RData')
mrgn.inf.alpha01 = loadRData('Simulation/data/int_and_child_filtered_data/alpha_01_selected_confs_results/MRGN_15confs_liberal_alpha01_inference_results.RData')


mrgn.inf.combined = c(mrgn.inf, mrgn.many.conf.inf)
mrgn.inf.alpha01.combined = c(mrgn.inf.alpha01, mrgn.mc.inf.alpha01)
gt.combined = c(params$model, many.conf.params$model)
gmac.inf.combined = c(master.table$GMAC.inference.at.05.cutoff,
                      master.table.many.conf$GMAC.inference.at.05.cutoff)

gmac.raw.combined = rbind.data.frame(cbind(gmac.cis$output.table$Cis_at_01_cutoff,
                                           gmac.trans$output.table$Trans_at_01_cutoff),
                                     cbind(gmac.cis.many.conf$output.table$Cis_at_01_cutoff,
                                           gmac.trans.many.conf$output.table$Trans_at_01_cutoff))

pm1 = params[,1:6]
pm2 = many.conf.params[,1:6]
colnames(pm1)=colnames(pm2) = colnames(params)[1:6]
params.combined = rbind(pm1,pm2)

mrpc.inf.all = loadRData("Manuscript/other/results_all_MRPC-ADDIS.RData")
