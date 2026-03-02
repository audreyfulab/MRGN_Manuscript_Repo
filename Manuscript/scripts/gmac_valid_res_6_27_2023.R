# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

library(MRGN)
library(ggpubr)
library(ggthemes)


path_supptabs = "Manuscript/supplementary_tables/"
path_suppfigs = "Manuscript/supplementary_figures/"
path_figs = "Manuscript/figures/"
path_tabs = "Manuscript/tables/"


x = loadRData("Simulation/GMAC_validation_SIM/gmac_valid_summary_results_table.RData")


x$`Significance Cutoff For Mediation` = x$Cutoff

color.codes = c("#0073C2FF","#EFC000FF")

A = ggplot(aes(x = `Simulated Confs`, y = `Type I Error`, color = `Significance Cutoff For Mediation`, 
           linetype = `Significance Cutoff For Mediation`), data = x)+
  geom_line(linewidth = 1.5)+
  geom_point(size = 3)+
  scale_fill_manual(values = color.codes)+
  scale_color_manual(values = color.codes)+
  scale_x_continuous(limits = c(2,16), breaks = seq(2, 16, 2))+
  theme_hc()+
  theme(legend.position = 'top', legend.text = element_text(size = 16),axis.text = element_text(size = 18),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 8, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 8, l = 0), size = 18))+
  xlab('Number of Confounders Simulated Per Trio')+
  ylab('Type I Error Rate For Mediation')


pdf(paste0(path_suppfigs, 'SF9_GMAC_Valid_TypeI_error.pdf'),
    height = 6, width = 8)

plot(A)

dev.off()