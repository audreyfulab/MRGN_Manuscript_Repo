
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

#library('ggpubr')
#library("MRGN")

library(ggpubr)
library(MRGN)
library(latex2exp)
library(ggthemes)

path_supptabs = "Manuscript/supplementary_tables/"
path_suppfigs = "Manuscript/supplementary_figures/"
path_figs = "Manuscript/figures/"
path_tabs = "Manuscript/tables/"

#===============================================investigate-Permutation=========================================
#read in trio data for wholeblood
#data.with.pcs = loadRData(file = paste0("GTEx/results/MRGN/", "WholeBloodpr.lnc.data.with.pcs.RData"))
#perform permutation on all trios
#perm.all = sapply(data.with.pcs, infer.trio, gamma = 0.5, use.perm = T, nperms = 1000)
#no permutation
#no.perm = sapply(data.with.pcs, infer.trio, gamma = 0.5, use.perm = F, nperms = 1000)

# Read in result for MRGN with permutation test applied to all trios (perm.all)
# --------------------------------------
# Alpha < 0.05
perm.all = loadRData(file = "Permutation_test_analysis/results/perm.all.trios.WB.liberal.confs.alpha05.RData")
write.table(perm.all, file = 'Permutation_test_analysis/results/GTEx_results_Permutation_all_trios_liberal_confs_and_alpha05.txt',
            sep = '\t', col.names = T, row.names = F)
# Alpha < 0.01
perm.all2 = loadRData(file = "Permutation_test_analysis/results/perm.all.trios.WB.RData")
write.table(perm.all2, file = 'Permutation_test_analysis/results/GTEx_results_Permutation_all_trios_Manuscript_settings.txt',
            sep = '\t', col.names = T, row.names = F)

# Read in results for MRGN without permutation test "no perm"
# # --------------------------------------
# Alpha < 0.05
no.perm = loadRData(file = "Permutation_test_analysis/results/no.perm.all.trios.WB.liberal.confs.alpha05.RData")
write.table(no.perm, file = 'Permutation_test_analysis/results/GTEx_results_No_Permutation_all_trios_Manuscript_settings.txt',
            sep = '\t', col.names = T, row.names = F)
# Alpha < 0.01
no.perm2 = loadRData(file = "Permutation_test_analysis/results/no.perm.all.trios.WB.RData")
write.table(no.perm2, file = 'Permutation_test_analysis/results/GTEx_results_No_Permutation_all_trios_liberal_confs_alpha05.txt',
            sep = '\t', col.names = T, row.names = F)
mrgn.tab = loadRData("GTEx/data/TrioTables/all_master_trio_tables.RData")$WholeBlood
write.table(mrgn.tab, file = 'GTEx/results/GTEx_Analysis_All_Results_Table_WholeBlood.txt',
            sep = '\t', col.names = T, row.names = F)


#no.perm  = loadRData(file = paste0("GTEx/results/MRGN/reg.res.wo.pseudo.RData"))$WholeBlood
#inf.no.perm = mrgn.tab$MRGN.Inferred.Model


rare.idx = which(unlist(perm.all[13,])<0.1)

write.table(mrgn.tab[rare.idx,], file = 'GTEx/results/GTEx_Analysis_All_Results_Table_WholeBlood_Rare_Trios_Only.txt',
            sep = '\t', col.names = T, row.names = F)

rare.var = rep("Minor Freq > 0.01", dim(perm.all)[2])
rare.var[rare.idx] = "Minor Freq < 0.01"
inf.perm = unlist(perm.all[14,])
inf.no.perm = unlist(no.perm[14,])

idx.mod.change = which(inf.perm != inf.no.perm)
mod.change = rep("No Change", dim(perm.all)[2])



after.perm = inf.perm[idx.mod.change]
before.perm = inf.no.perm[idx.mod.change]
xx = paste(before.perm, after.perm, sep = " --> ")

classify.change.types = function(x){
  
  if (x[1] == "M0.1" | x[1] == "M0.2"){
    if (x[2] == "M1.1" | x[2] == "M1.2"){
      y = "liberal"
    }else{
      y = "conservative"
    }
  }else if(x[1] == "M1.1" | x[1] == "M1.2"){
    y = "conservative"
    
  }else if(x[1] == "M2.2"){
    y = "conservative"
    
  }else if(x[1] == "M3"){
    
    if(x[2] == "M2.2" | x[2] == "M4"){
      y = "liberal"
    }else{
      y = "conservative"
    }
    
  }else if (x[1] == "M4"){
    y = "conservative"
    
  }else if (x[1] == "Other"){
    y = "liberal"
  }
  
  
}


change.type = apply(cbind.data.frame(before.perm, after.perm), 1, classify.change.types)

mod.change[idx.mod.change] = xx
mod.change2 = mod.change
mod.change2[which(mod.change != "No Change")] = "Model Changed"

counts = summary(as.factor(xx))
nm = strsplit(names(counts), ' --> ')
chty = unlist(lapply(nm, classify.change.types))

mod.change.counts = cbind.data.frame(counts, `type of model change` = names(counts), 
                                     `Category of Change` = chty)

# ct = summary(as.factor(change.type))
# mod.change.types = cbind.data.frame(count = ct, `category of change` = names(ct))


#combine all permutation pvalues into a matrix
pvalues.mat = cbind.data.frame(pb12.with.perm = unlist(perm.all[8,]),
                               pb12.no.perm = unlist(no.perm[8,]),
                               pb22.with.perm = unlist(perm.all[10,]),
                               pb22.no.perm = unlist(no.perm[10,]),
                               Allele.freq = unlist(perm.all[13,]),
                               `Types of Model Changes` = as.factor(mod.change),
                               Inference = as.factor(mod.change2))

pvalues.mat$diff.perm.p = abs(pvalues.mat$pb12.with.perm - pvalues.mat$pb22.with.perm)



h=0.01
A=ggplot(data = pvalues.mat)+
  geom_point(aes(x=log10(pb12.no.perm), y=log10(pb12.with.perm), color=Inference, fill = Inference),
             shape = 19, alpha = 0.6)+
  color_palette(palette = "jco")+
  theme_hc()+
  theme(legend.position = "right")+
  theme(legend.text = element_text(size = 16),axis.text = element_text(size = 18),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
  labs(y = TeX('$\\log_{10}$ pvalue $\\beta_{12}$ with permutation'),
       x = TeX('$\\log_{10}$ pvalue $\\beta_{12}$ without permutation'))+
  scale_x_continuous(limits = c(-15,0))+
  scale_y_continuous(limits = c(-15,0))+
  geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  geom_vline(xintercept = log10(h), color="red", linetype="dashed")+
  geom_hline(yintercept = log10(h), color="red", linetype="dashed")

h=0.01
B=ggplot(data = pvalues.mat)+
  geom_point(aes(x=log10(pb22.no.perm), y=log10(pb22.with.perm), color=Inference, fill = Inference), 
             shape = 19, alpha = 0.6)+
  color_palette(palette = "jco")+
  theme_hc()+
  theme(legend.position = "right")+
  theme(legend.text = element_text(size = 16),axis.text = element_text(size = 18),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
  labs(y = TeX('$\\log_{10}$ pvalue $\\beta_{22}$ with permutation'),
       x = TeX('$\\log_{10}$ pvalue $\\beta_{22}$ without permutation'))+
  scale_x_continuous(limits = c(-15,0))+
  scale_y_continuous(limits = c(-15,0))+
  geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  geom_vline(xintercept = log10(h), color="red", linetype="dashed")+
  geom_hline(yintercept = log10(h), color="red", linetype="dashed")

#pdf("Permutation_test_analysis/results/Perm_check_plot.pdf")
C = ggarrange(A, B, labels = c("A", "B"), nrow = 1, ncol = 2, common.legend = T, legend = "top")
plot(C)

n = dim(pvalues.mat)[1]
pvals.longform = cbind.data.frame(`Permutation Test` = c(rep("b12", n), rep("b22", n)),
                                  pvalue = c(pvalues.mat$pb12.with.perm, pvalues.mat$pb22.with.perm),
                                  `Allele Frequency` = rep(pvalues.mat$Allele.freq, 2),
                                  `Inferred Model Status` = rep(pvalues.mat$Inference, 2))


D=ggplot(data = pvals.longform, aes(x=`Inferred Model Status`, y=`Allele Frequency`))+
  xlab("")+
  ylab(TeX('Minor Allele Frequency'))+
  theme(legend.position = "bottom")+
  stat_compare_means(method = "t.test", label.x = 1.3, label.y = 0.47, color = 'red')+
  stat_compare_means(label.x = 1.3, label.y = 0.42, color = 'red')+
  #scale_y_continuous(limits = c(0, 0.6))+
  theme_hc()+
  theme(legend.text = element_text(size = 16),axis.text = element_text(size = 18),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
  #scale_fill_discrete(labels=c(TeX("$\\H_0: \ \\beta_{12} = 0$"), TeX("$\\H_0: \ \\beta_{22} = 0$")))+
  geom_boxplot(fill = "#7AA6DCFF")
plot(D)


color.codes = c("#0073C2FF","#868686FF")
E = ggplot(data = mod.change.counts, aes(x = `type of model change`, y = counts, fill = `Category of Change`))+
  theme_hc()+
  scale_fill_manual(values = color.codes)+
  scale_color_manual(values = color.codes)+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  geom_bar(stat = "identity", color = "black")+
  theme(legend.position = 'top', legend.text = element_text(size = 16),axis.text = element_text(size = 16),
       legend.title = element_text(size = 16),
       axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
       axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
  xlab("Type of Model Change After Permutation")+
  ylab("Count")
#geom_bar(aes(x = `category of change`, y = count), data = mod.change.types, stat = "identity", fill = c("red", "blue"))
plot(E)

final.plot = ggarrange(C, ggarrange(E, D, labels = c("C", "D"), nrow = 1, ncol = 2), nrow = 2, ncol = 1, vjust = -1)+theme(plot.margin = unit(c(0, 0, 0.45, 0), "cm"))

pdf(paste0(path_suppfigs, "SF3_GTEx_Permutation_and_Model_Changes.pdf"), 
    height = 12, width = 12)
plot(final.plot)
dev.off()

# NOTE: External path from Doctoral_Thesis repo - update if needed
# jpeg("Manuscript/supplementary_figures/SF3_GTEx_Permutation_and_Model_Changes.jpeg", 
#     height = 10, width = 12, units = 'in', res = 1200, quality = 90)
# plot(final.plot)
# dev.off()

