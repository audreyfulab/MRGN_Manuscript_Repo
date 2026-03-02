# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

# model mispecification simulation scenario pvalue plots
library(prodlim)
library(ggpubr)
library(MRGN)
library(ggthemes)
library(gridExtra)
library(latex2exp)


path_supptabs = "Manuscript/supplementary_tables/"
path_suppfigs = "Manuscript/supplementary_figures/"
path_figs = "Manuscript/figures/"
path_tabs = "Manuscript/tables/"

source('Manuscript/scripts/MRGN_write_up_helper_functions.R')
mrgn.tab = loadRData("GTEx/data/TrioTables/all_master_trio_tables.RData")$WholeBlood
# NOTE: This file originates from the MRPC_support repository. Place it at the path below.
sim=read.csv(file="Manuscript/other/TRIOS_imbalanced_genotypes_WB_Updated_12_10_2021.csv", header = T)

A = mrgn.tab[,c(1,2,8)]
B = cbind.data.frame(sim$snp, sim$cis.gene, sim$trans.gene)
colnames(A) = colnames(B) = c("SNP", "Cis.Gene", "Trans.Gene")
ix = na.omit(apply(A, 1, row.match, table = B))

sim2 = sim[ix,]
sim2$ha[which(is.na(sim2$ha))] = 0
sim2$alternative=(sim2$het+2*sim2$ha)/(2*rowSums(sim2[,37:39]))
sim2$reference=1-sim2$alternative
cv1=rep(0, dim(sim2)[1])
cv2=replace(cv1, which(sim2$alternative<=0.1 | sim2$reference<=0.1), 1)
sim2$`Minor Allele Frequency` = cv2
sim2$`Minor Allele Frequency`[which(cv2 == 1)] = "Less Than 10%"
sim2$`Minor Allele Frequency`[which(cv2 == 0)] = "Greater Than 10%"


color.codes2 = c("#0073C2FF", "#EFC000FF", "#868686FF")
require(grid)
h=0.05
A3=ggplot(data = sim2)+
  geom_point(aes(x=log10(Nominal.p.STM), y=log10(STM.median.p), fill=`Minor Allele Frequency`),
             shape = 21, alpha = 0.8, size = 2.5, color = 'black')+
  theme(legend.position = "none")+
  scale_fill_manual(values = color.codes2)+
  scale_color_manual(values = color.codes2)+
  theme_hc()+
  labs(y = expression(log10('Median Sim. Para P-value')),
       x = expression(log10('True Nominal P-value')),
       subtitle = "Scenario: STMS")+
  scale_x_continuous(limits = c(-16,0))+
  scale_y_continuous(limits = c(-16,0))+
  geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  geom_vline(xintercept = log10(h), color="red", linetype="dashed")+
  geom_hline(yintercept = log10(h), color="red", linetype="dashed")+
  theme(legend.text = element_text(size = 16),axis.text = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))



B3=ggplot(data = sim2)+
  geom_point(aes(x=log10(Nominal.p.LTM), y=log10(LTM.median.p), fill=`Minor Allele Frequency`),
             shape = 21, alpha = 0.8, size = 2.5, color = 'black')+
  theme(legend.position = "none")+
  scale_fill_manual(values = color.codes2)+
  scale_color_manual(values = color.codes2)+
  theme_hc()+
  labs(y = expression(log10('Median Sim. Para P-value')),
       x = expression(log10('True Nominal P-value')),
       subtitle = "Scenario: LTMS")+
  scale_x_continuous(limits = c(-16,0))+
  scale_y_continuous(limits = c(-16,0))+
  geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  geom_vline(xintercept = log10(h), color="red", linetype="dashed")+
  geom_hline(yintercept = log10(h), color="red", linetype="dashed")+
  theme(legend.text = element_text(size = 16), axis.text = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))

C3=ggplot(data = sim2)+
  geom_point(aes(x=log10(Nominal.p.TGM), y=log10(TGM.median.p), fill=`Minor Allele Frequency`),
             shape = 21, alpha = 0.8, size = 2.5, color = 'black')+
  theme(legend.position = "none")+
  scale_fill_manual(values = color.codes2)+
  scale_color_manual(values = color.codes2)+
  theme_hc()+
  labs(y = expression(log10('Median Sim. Para P-value')),
       x = expression(log10('True Nominal P-value')),
       subtitle = "Scenario: SIMS")+
  scale_x_continuous(limits = c(-16,0))+
  scale_y_continuous(limits = c(-16,0))+
  geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  geom_vline(xintercept = log10(h), color="red", linetype="dashed")+
  geom_hline(yintercept = log10(h), color="red", linetype="dashed")+
  theme(legend.text = element_text(size = 16), axis.text = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))

D3=ggplot(data = sim2)+
  geom_point(aes(x=log10(Nominal.p.GSS), y=log10(GSS.median.p), fill=`Minor Allele Frequency`),
             shape = 21, alpha = 0.8, size = 2.5, color = 'black')+
  theme(legend.position = "none")+
  scale_fill_manual(values = color.codes2)+
  scale_color_manual(values = color.codes2)+
  theme_hc()+
  labs(y = expression(log10('Median pvalue from Simulation')),
       x = expression(log10('True pvalue')),
       subtitle = "Scenario: CMS")+
  scale_x_continuous(limits = c(-16,0))+
  scale_y_continuous(limits = c(-16,0))+
  geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  geom_vline(xintercept = log10(h), color="red", linetype="dashed")+
  geom_hline(yintercept = log10(h), color="red", linetype="dashed")+
  theme(legend.text = element_text(size = 16), axis.text = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))

figure <- ggarrange(A3 + rremove("ylab") + rremove("xlab"), B3 + rremove("ylab") + rremove("xlab"), 
                    D3 + rremove("ylab") + rremove("xlab"), C3 + rremove("ylab") + rremove("xlab"), 
                    # remove axis labels from plots
                    labels = c("A","B","C","D"),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "bottom",
                    align = "hv", 
                    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))+theme(plot.margin = unit(c(0, 0, 0.45, 0), "cm"))

pdf(paste0(path_suppfigs, "SF6_Model_Misspecification.pdf"), 
    height = 10, width = 12)
figure2 = annotate_figure(figure, left = textGrob(TeX('Median Observed $\\log_{10}$ Pvalue for $\\beta_{22}$ From Parametric Test'), 
                                                  rot = 90, vjust = 0.3, gp = gpar(cex = 1.5)),
                          bottom = textGrob(TeX('Observed $\\log_{10}$ Pvalue for $\\beta_{22}$ From Permutation Test'), 
                                            gp = gpar(cex = 1.5)))
plot(figure2)

dev.off()


ttests = dim(sim2)[1]
#create FP TP table
typeII.param = c(sum(sim2$STM.median.p>0.05)/ttests, sum(sim2$LTM.median.p>0.05)/ttests, 
                 sum(sim2$GSS.median.p>0.05)/ttests, sum(sim2$TGM.median.p>0.05)/ttests) 

typeII.perm = c(sum(sim2$Nominal.p.STM>0.05)/ttests, sum(sim2$Nominal.p.LTM>0.05)/ttests, 
                sum(sim2$Nominal.p.GSS>0.05)/ttests, sum(sim2$Nominal.p.TGM>0.05)/ttests)




cn = c("Test Type", "Scenario STMS", "Scenario LTMS", "Scenario CMS", "Scenario SIMS")
rn = c("Permutation Test (Type II Error)", "Parametric Test (Type II Error)")



table1 = cbind.data.frame(rn, rbind.data.frame(typeII.perm, typeII.param))
colnames(table1) = cn

write.csv(table1, file = paste0(path_supptabs,'/ST10_Model_miss_Type_II_error_rates.csv'))

