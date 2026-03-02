
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

#sp = 'Manuscript/supplementary_figures/Rare_trio_distributions/'
#sp = 'Manuscript/supplementary_figures/Rare_trio_distributions/'
#sp2 = 'Manuscript/supplementary_tables/'
path_supptabs = "Manuscript/supplementary_tables/"
path_suppfigs = "Manuscript/supplementary_figures/"
path_figs = "Manuscript/figures/"
path_tabs = "Manuscript/tables/"

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(MRGN)
library(cowplot)
#library('ggplot2')
#library('ggpubr')
#library('ggthemes')
#tissue.names=read.csv("GTEx/data/tissuenames.csv", header = T)
tissue.names=read.csv("GTEx/data/tissuenames.csv", header = T)
#library('MRGN')
#library('cowplot')

top5=c(48,40,33,6,1)
tissues.vec=tissue.names[top5, 2]
#perm.all = loadRData(file = "Permutation_test_analysis/results/perm.all.trios.WB.liberal.confs.alpha05.RData")
#trios.pr.lr.only = loadRData(file = paste0("GTEx/data/trios_subset_data/all.data.unqiue.snps.pclrna.only.",tissues.vec[1],".RData"))
perm.all = loadRData(file = "Permutation_test_analysis/results/perm.all.trios.WB.liberal.confs.alpha05.RData")
trios.pr.lr.only = loadRData(file = paste0("GTEx/data/trios_subset_data/all.data.unqiue.snps.pclrna.only.",tissues.vec[1],".RData"))

number.of.trios = dim(trios.pr.lr.only)[2]/3
trio.idx.mat = matrix(c(1:dim(trios.pr.lr.only)[2]), nrow = number.of.trios, ncol = 3, byrow = T)

#mrgn.tab = loadRData("GTEx/data/TrioTables/all_master_trio_tables.RData")$WholeBlood
mrgn.tab = loadRData("GTEx/data/TrioTables/all_master_trio_tables.RData")$WholeBlood
rare.idx = which(unlist(perm.all[13,])<0.1)
no.match = NULL
trio.idx = mrgn.tab$Trio.sub.index[rare.idx]

#find trios whose model was different before and after permutation
model.change.idx = which(mrgn.tab$MRGN.Inferred.Model.libconf.alpha05.all.perm != mrgn.tab$MRGN.Inferred.Model.libconf.alpha05.no.perm)
#how many of these also have rare alleles?
rare.changed = na.omit(match(model.change.idx, rare.idx))

rare.skew = cbind.data.frame(case = rep('Minor Allele Freq < 0.1', length(rare.idx)),
                             index = trio.idx,
                             as.data.frame(matrix(0, nrow = length(rare.idx), ncol = 4)))
colnames(rare.skew) = c('minor.freq', 'index', 'cis.skew', 'cis.kurtosis', 'trans.skew', 'trans.kurtosis')

for(i in 1:length(rare.idx)){
  
  trio = as.data.frame(na.omit(trios.pr.lr.only[, trio.idx.mat[rare.idx[i],]]))
  cn = colnames(trio)
  colnames(trio) = c('SNP Genotype', 'Cis Gene Expression', 'Trans Gene Expression')
  
  rare.skew$cis.skew[i] = propagate::skewness(trio[,2])
  rare.skew$trans.skew[i] = propagate::skewness(trio[,3])
  rare.skew$cis.kurtosis[i] = propagate::kurtosis(trio[,2])
  rare.skew$trans.kurtosis[i] = propagate::kurtosis(trio[,3])
  
  print('====================================================================')
  print(paste0('Cis Gene Skew = ', round(rare.skew$cis.skew[i],4), 
               ' : Cis Gene Kurtosis = ', round(rare.skew$cis.kurtosis[i],4)))
  print('--------------------------------------------------------------------')
  print(paste0('Trans Gene Skew = ', round(rare.skew$trans.skew[i],4), 
               ' : Trans Gene Kurtosis = ', round(rare.skew$trans.kurtosis[i],4)))
  
  newlabs = paste(c('SNP ID', 'Cis Gene', 'Trans Gene'), cn, sep = ': ')
  if(get.freq(trio[,1]) == perm.all[13,rare.idx[i]]){
    
    if(trio.idx[i] > 0){
      color.codes = c("#0073C2FF","#868686FF","#EFC000FF")
      k = ceiling(sqrt(dim(trio)[1]))
      trio$`SNP Genotype` = as.factor(trio$`SNP Genotype`)
      print(paste0('Plotting trio #', trio.idx[i]))
      A = ggplot(data = as.data.frame(trio), aes(x = `Cis Gene Expression`, y = `Trans Gene Expression`,
                                                 color = `SNP Genotype`, fill = `SNP Genotype`))+
        geom_point(size = 3, shape = 21, alpha = 0.7, color = 'black')+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 12))+
        scale_fill_manual(values = color.codes)+
        scale_color_manual(values = color.codes)+
        xlab(newlabs[2])+
        ylab(newlabs[3])+
        ggtitle(newlabs[1])+
        theme_classic2()+
        theme(title = element_text(size = 16),
              legend.text = element_text(size = 16),
              axis.text = element_text(size = 18),
              legend.title = element_text(size = 16),
              axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
              axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))
      
      B = ggplot(data = as.data.frame(trio))+
        geom_histogram(aes(x = `Cis Gene Expression`),
                       color = 'black', fill = '#0073C2FF', bins = k)+
        theme(legend.text = element_text(size = 16),
              axis.text = element_text(size = 18),
              legend.title = element_text(size = 16),
              axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
              axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
        ylab('Frequency')+
        theme_hc()
      
      C = ggplot(data = as.data.frame(trio))+
        geom_histogram(aes(x = `Trans Gene Expression`),
                       color = 'black', fill = '#0073C2FF', bins = k)+
        theme(legend.text = element_text(size = 16),axis.text = element_text(size = 18),
              legend.title = element_text(size = 16),
              axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
              axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
        ylab('Frequency')+
        theme_hc()
      
      D = ggarrange(B,C, nrow = 1)
      pdf(paste0(path_suppfigs,'rare_trio_distributions/trio_number_', trio.idx[i], '.pdf'), height = 10, width = 10)
      final = ggarrange(A,D, nrow = 2, ncol = 1, common.legend = T)+theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm"))
      plot(final)
      dev.off()
    }
  }else{
    warning('Frequencies did NOT match!!! ... identifying trio ...')
    no.match = append(no.match, rare.idx[i])
  }
}




common.idx = which(unlist(perm.all[13,])>=0.1)
trio.idx = mrgn.tab$Trio.sub.index[common.idx]
common.skew = cbind.data.frame(case = rep('Minor Allele Freq > 0.1', length(common.idx)),
                               index = trio.idx,
                             as.data.frame(matrix(0, nrow = length(common.idx), ncol = 4)))
colnames(common.skew) = c('minor.freq', 'index', 'cis.skew', 'cis.kurtosis', 'trans.skew', 'trans.kurtosis')

for(i in 1:length(common.idx)){
  
  
  trio = as.data.frame(na.omit(trios.pr.lr.only[, trio.idx.mat[common.idx[i],]]))
  cn = colnames(trio)
  colnames(trio) = c('SNP Genotype', 'Cis Gene Expression', 'Trans Gene Expression')
  
  
  common.skew$cis.skew[i] = propagate::skewness(trio[,2])
  common.skew$trans.skew[i] = propagate::skewness(trio[,3])
  common.skew$cis.kurtosis[i] = propagate::kurtosis(trio[,2])
  common.skew$trans.kurtosis[i] = propagate::kurtosis(trio[,3])
  
  print('====================================================================')
  print(paste0('Cis Gene Skew = ', round(common.skew$cis.skew[i],4), 
               ' : Cis Gene Kurtosis = ', round(common.skew$cis.kurtosis[i],4)))
  print('--------------------------------------------------------------------')
  print(paste0('Trans Gene Skew = ', round(common.skew$trans.skew[i],4), 
               ' : Trans Gene Kurtosis = ', round(common.skew$trans.kurtosis[i],4)))
  
  newlabs = paste(c('SNP ID', 'Cis Gene', 'Trans Gene'), cn, sep = ': ')
  if(get.freq(trio[,1]) == perm.all[13,common.idx[i]]){
    color.codes = c("#0073C2FF","#868686FF","#EFC000FF")
    k = ceiling(sqrt(dim(trio)[1]))
    trio$`SNP Genotype` = as.factor(trio$`SNP Genotype`)
    print(paste0('Plotting trio #', common.idx[i]))
    # A = ggplot(data = as.data.frame(trio), aes(x = `Cis Gene Expression`, y = `Trans Gene Expression`,
    #                                            color = `SNP Genotype`, fill = `SNP Genotype`))+
    #   geom_point(size = 3, shape = 21, alpha = 0.7, color = 'black')+
    #   theme(axis.text = element_text(size = 12),
    #         axis.title = element_text(size = 12))+
    #   scale_fill_manual(values = color.codes)+
    #   scale_color_manual(values = color.codes)+
    #   xlab(newlabs[2])+
    #   ylab(newlabs[3])+
    #   ggtitle(newlabs[1])+
    #   theme_classic2()+
    #   theme(legend.text = element_text(size = 16),axis.text = element_text(size = 18),
    #         axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
    #         axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))
    # 
    # B = ggplot(data = as.data.frame(trio))+
    #   geom_histogram(aes(x = `Cis Gene Expression`),
    #                  color = 'black', fill = '#0073C2FF', bins = k)+
    #   theme(legend.text = element_text(size = 16),axis.text = element_text(size = 18),
    #         axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
    #         axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
    #   ylab('Frequency')+
    #   theme_hc()
    # 
    # C = ggplot(data = as.data.frame(trio))+
    #   geom_histogram(aes(x = `Trans Gene Expression`),
    #                  color = 'black', fill = '#0073C2FF', bins = k)+
    #   theme(legend.text = element_text(size = 16),axis.text = element_text(size = 18),
    #         axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
    #         axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
    #   ylab('Frequency')+
    #   theme_hc()
    # 
    # D = ggarrange(B,C, nrow = 1)
    # pdf(paste0(sp,'/rare_trio_distributions/common_allele_trios/','trio_number_', trio.idx[i], '.pdf'), height = 6, width = 10)
    # final = ggarrange(A,D, nrow = 2, ncol = 1, common.legend = T)
    # plot(final)
    # dev.off()
  }else{
    warning('Frequencies did NOT match!!! ... identifying trio ...')
    no.match = append(no.match, common.idx[i])
  }
}




cis.skew.all = cbind.data.frame(gene.position = rep('Cis', number.of.trios),
                            rbind.data.frame(common.skew[,c(1,3)], rare.skew[,c(1,3)]))
trans.skew.all = cbind.data.frame(gene.position = rep('Trans', number.of.trios),
                                  rbind.data.frame(common.skew[, c(1,5)], rare.skew[,c(1,5)]))
colnames(cis.skew.all) = colnames(trans.skew.all) = c('Gene Type','Variant Frequency', 'skew.stat')
final.skew = rbind(cis.skew.all, trans.skew.all)

color.codes = c("#0073C2FF","#EFC000FF")
A = ggplot(data = final.skew, aes(x = skew.stat,# y = after_stat(density), 
                                  fill = `Gene Type`, color = `Gene Type`))+
  geom_histogram(color = 'black', bins = 20, alpha = 0.8, position = 'stack')+
  #geom_density(alpha = 0, linewidth = 0.9)+
  scale_fill_manual(values = color.codes)+
  scale_color_manual(values = color.codes)+
  facet_wrap(~`Variant Frequency`)+
  theme(legend.text = element_text(size = 16),axis.text = element_text(size = 18),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
  xlab('Sample Skewness')+
  ylab('Frequency')+
  theme_hc()
plot(A)
  


cis.kurt.all = cbind.data.frame(gene.position = rep('Cis', number.of.trios),
                                rbind.data.frame(common.skew[,c(1,4)], rare.skew[,c(1,4)]))
trans.kurt.all = cbind.data.frame(gene.position = rep('Trans', number.of.trios),
                                  rbind.data.frame(common.skew[, c(1,6)], rare.skew[,c(1,6)]))
colnames(cis.kurt.all) = colnames(trans.kurt.all) = c('Gene Type','Variant Frequency', 'kurtosis')
final.kurtosis = rbind(cis.kurt.all, trans.kurt.all)


summary_hp_tests = cbind.data.frame(`Gene` = c('Cis', 'Trans', 'Cis', 'Trans'),
                                    `Summary Statistic` = c('Sample Skewness', 'Sample Skewness', 
                                                            'Sample Excess Kurtosis', 'Sample Excess Kurtosis'),
                                    Test = rep('Minor Allele < 0.1 (Rare) - Minor Allele > 0.1 (Common)', 4),
                                    `Mean Rare` = c(mean(subset(cis.skew.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$skew.stat),
                                                    mean(subset(trans.skew.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$skew.stat),
                                                    mean(subset(cis.kurt.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$kurtosis),
                                                    mean(subset(trans.kurt.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$kurtosis)),
                                    `SD Rare` = c(sd(subset(cis.skew.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$skew.stat),
                                                  sd(subset(trans.skew.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$skew.stat),
                                                  sd(subset(cis.kurt.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$kurtosis),
                                                  sd(subset(trans.kurt.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$kurtosis)),
                                    `Mean Common` = c(mean(subset(cis.skew.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$skew.stat),
                                                      mean(subset(trans.skew.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$skew.stat),
                                                      mean(subset(cis.kurt.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$kurtosis),
                                                      mean(subset(trans.kurt.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$kurtosis)),
                                    `SD Common` = c(sd(subset(cis.skew.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$skew.stat),
                                                    sd(subset(trans.skew.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$skew.stat),
                                                    sd(subset(cis.kurt.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$kurtosis),
                                                    sd(subset(trans.kurt.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$kurtosis)),
                                    `Wilcoxen Sum-Rank Statistic` = c(wilcox.test(x = subset(cis.skew.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$skew.stat,
                                                                                  y = subset(cis.skew.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$skew.stat, 
                                                                                  alternative = 'two.sided')$statistic,
                                                                      wilcox.test(x = subset(trans.skew.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$skew.stat,
                                                                                  y = subset(trans.skew.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$skew.stat, 
                                                                                  alternative = 'two.sided')$statistic,
                                                                      wilcox.test(x = subset(cis.kurt.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$kurtosis,
                                                                                  y = subset(cis.kurt.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$kurtosis, 
                                                                                  alternative = 'two.sided')$statistic,
                                                                      wilcox.test(x = subset(trans.kurt.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$kurtosis,
                                                                                  y = subset(trans.kurt.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$kurtosis, 
                                                                                  alternative = 'two.sided')$statistic),
                                                                      
                                    `Pvalue` = c(wilcox.test(x = subset(cis.skew.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$skew.stat,
                                                             y = subset(cis.skew.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$skew.stat, 
                                                             alternative = 'two.sided')$p.value,
                                                 wilcox.test(x = subset(trans.skew.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$skew.stat,
                                                             y = subset(trans.skew.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$skew.stat, 
                                                             alternative = 'two.sided')$p.value,
                                                 wilcox.test(x = subset(cis.kurt.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$kurtosis,
                                                             y = subset(cis.kurt.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$kurtosis, 
                                                             alternative = 'two.sided')$p.value,
                                                 wilcox.test(x = subset(trans.kurt.all, `Variant Frequency` == 'Minor Allele Freq < 0.1')$kurtosis,
                                                             y = subset(trans.kurt.all, `Variant Frequency` == 'Minor Allele Freq > 0.1')$kurtosis, 
                                                             alternative = 'two.sided')$p.value))
                                    
summary_hp_tests$`Bonferroni Corrected Pvalue` = p.adjust(summary_hp_tests$Pvalue, method = 'bonferroni')

A2 = ggplot(data = final.kurtosis, aes(x = kurtosis,# y = after_stat(density), 
                                       fill = `Gene Type`, color = `Gene Type`))+
  geom_histogram(color = 'black', bins = 20, alpha = 0.8, position = 'stack')+
  #geom_density(alpha = 0, linewidth = 0.9)+
  scale_fill_manual(values = color.codes)+
  scale_color_manual(values = color.codes)+
  facet_wrap(~`Variant Frequency`)+
  theme(legend.text = element_text(size = 16),axis.text = element_text(size = 18),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
  xlab('Sample (Excess) Kurtosis')+
  ylab('Frequency')+
  theme_hc()
plot(A2)

pdf(paste0(path_suppfigs, 'SF8_Trio_Skew_Stat_Distribtuions.pdf'), height = 8, width = 10)
plot(ggarrange(A, A2, nrow = 2, common.legend = T)+theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")))
dev.off()


write.csv(summary_hp_tests, file = paste0(path_supptabs, 'ST_sign_rank_tests_skewness_in_trans_genes.csv'))






final.skew.mat = rbind.data.frame(common.skew, rare.skew)

















































