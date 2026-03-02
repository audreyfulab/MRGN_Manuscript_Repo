# NOTE: Set your working directory to the repository root before running this script.
# NOTE: Originally run on HPC cluster. Some external data paths may need adjustment.

#library('ggpubr')
library("MRGN")

#===============================================investigate-Permutation=========================================
#read in trio data for wholeblood
#data.with.pcs = loadRData(fileName=paste0("GTEx/data/final_trios_with_PCs/data.with.PCs.WholeBlood.RData"))
data.with.pcs.alpha01 = loadRData(file ="GTEx/trios_data_prlnc_liberal_conf_sel/data_with_confs/trios.with.confs.WholeBlood.RData")
#perform permutation on all trios with liberal confs
#perform permutation on all trios
perm.all = sapply(data.with.pcs.alpha01, infer.trio, gamma = 0.5, use.perm = T, nperms = 1000)
save(perm.all, file = "Permutation_test_analysis/results/perm.all.trios.WB.liberal.confs.alpha01.RData")
#no permutation
no.perm = sapply(data.with.pcs.alpha01, infer.trio, gamma = 0.5, use.perm = F, nperms = 1000)
save(no.perm, file = "Permutation_test_analysis/results/no.perm.all.trios.WB.liberal.confs.alpha01.RData")
