# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

library(MRGN)

path_supptabs = "Manuscript/supplementary_tables/"
path_suppfigs = "Manuscript/supplementary_figures/"
path_figs = "Manuscript/figures/"
path_tabs = "Manuscript/tables/"


x = loadRData(file = "Simulation/data/diagnostics/times_to_select_confs.RData")

data.longform = cbind(Implementation = c(rep("GMAC", 6), rep("MRGN", 6)),
                                 rbind(setNames(x[,-4], c(colnames(x)[1:2],"Time.in.Min")),
                                       setNames(x[,-3], c(colnames(x)[1:2],"Time.in.Min"))))


data.longform2 = data.longform
data.longform2$Number_of_covs_in_pool[data.longform2$Number_of_covs_in_pool==100] = "100 candidate covariates"
data.longform2$Number_of_covs_in_pool[data.longform2$Number_of_covs_in_pool==500] = "500 candidate covariates"
data.longform2$Number_of_covs_in_pool[data.longform2$Number_of_covs_in_pool==1000] = "1000 candidate covariates"

data.longform2$Number_of_covs_in_pool_f = as.factor(data.longform2$Number_of_covs_in_pool)

library(ggpubr)
library(grid)
library(ggthemes)

# A=ggplot(data = data.longform, aes(x = as.factor(Number_of_covs_in_pool), y = Time.in.Min, fill = Implementation))+
#   geom_bar(stat = "identity", color = "black", position=position_dodge())+
#   scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
#   theme_minimal()+
#   ylab("Time To Compute Confounders (Mins)")+
#   xlab("Number of Covariates In Candidate Pool")
B=ggplot(data = data.longform2, aes(x = as.factor(Number_of_trios), y = Time.in.Min, fill = Implementation))+
  geom_bar(stat = "identity", color = "black", position=position_dodge())+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_hc()+
  ylab("Time To Compute Confounders (Mins)")+
  xlab("Number of Trios")+
  theme(legend.position = "top", legend.text = element_text(size = 16),
        axis.text = element_text(size = 16), legend.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 18))+
  facet_grid(.~factor(Number_of_covs_in_pool, levels = c("100 candidate covariates","500 candidate covariates","1000 candidate covariates")))+
  #geom_smooth()+
  theme(panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1))

# figure = ggarrange(A+rremove("ylab"), B+rremove("ylab"), labels = c("A", "B"), nrow = 1, ncol = 2,
#           heights = c(8,8), widths = c(10,10), common.legend = T, legend = "top")
#
# annotate_figure(figure, left = textGrob("Time To Compute Confounders (Min)", rot = 90, vjust = 1, gp = gpar(cex = 0.95)))
#





pdf(paste0(path_suppfigs, "SF1_comp_times_for_conf_selection_updated.pdf"), onefile = F)

plot(B)

dev.off()









