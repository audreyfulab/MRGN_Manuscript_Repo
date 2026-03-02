
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

library('ggpubr')
library('ggthemes')
library('gridExtra')
library('MRGN')

path_supptabs = "Manuscript/supplementary_tables/"
path_suppfigs = "Manuscript/supplementary_figures/"
path_figs = "Manuscript/figures/"
path_tabs = "Manuscript/tables/"

reg.res.L = loadRData(file = "GTEx/results/MRGN/reg.res.wo.pseudo.list.all.tissues.RData")
reg.res.libconf = loadRData(file = 'GTEx/trios_data_prlnc_liberal_conf_sel/infer_trio_results/infer.trio.results.WholeBlood.RData')

source('Manuscript/scripts/MRGN_write_up_helper_functions.R')
mrgn.tab = loadRData("GTEx/data/WholeBloodmaster_table.RData")
gmac.tab = loadRData("GTEx/results/GMAC/gmac.results.tables.combined.RData")
gmac.tab = gmac.tab$WholeBlood

#gmac.cis.tab = loadRData("GTEx/results/trio_tables/WholeBlood.gmac.results.tables.cis.RData")
#gmac.trans.tab = loadRData("GTEx/results/trio_tables/WholeBlood.gmac.results.tables.trans.RData")

write.csv(mrgn.tab, file = paste0(path_supptabs, 'ST_GTEx_all_trios_master.csv'))
#for plots showing freqency of 5 models
sum.mrgn = summary(as.factor(convert.cats(mrgn.tab$MRGN.Inferred.Model)))
sum.mrgn.libconf = summary(as.factor(convert.cats(mrgn.tab$ MRGN.Inferred.Model.libconf.alpha05)))
sum.mrpc = summary(as.factor(convert.cats(mrgn.tab$MRPC.Addis.Inferred.Model)))
models = names(sum.mrgn)
#bind to table
model.res = cbind.data.frame(Method = c(rep("MRGN", 6),
                                        rep("GMAC", 6), rep("MRPC-ADDIS", 6)),
                             Count = c(sum.mrgn.libconf,
                                       rep(0, 6), sum.mrpc),
                             `Inferred Model` = rep(models, 3))


#for cis mediation trans mediation breakdown
m1.mrgn = subset(mrgn.tab, MRGN.Inferred.Model.libconf.alpha05 == "M1.1" | MRGN.Inferred.Model.libconf.alpha05 == "M1.2")$MRGN.Inferred.Model.libconf.alpha05
m1.mrpc = subset(mrgn.tab, MRPC.Addis.Inferred.Model == "M1.1" | MRPC.Addis.Inferred.Model == "M1.2")$MRPC.Addis.Inferred.Model
#bind to table
model.res2 = cbind.data.frame(Method = c(rep("MRGN", 2),
                                         rep("GMAC", 2), rep("MRPC-ADDIS", 2)),
                             Count = c(summary(as.factor(m1.mrgn)),
                                       rep(0, 2), summary(as.factor(m1.mrpc))),
                             `Inferred Model` = c("Via The Cis Gene", "Via the Trans Gene"))

model.res$Method = factor(model.res$Method, levels = c("MRGN", "GMAC", "MRPC-ADDIS"))
model.res2$Method = factor(model.res2$Method, levels = c("MRGN", "GMAC", "MRPC-ADDIS"))


#5models mrgn and mrpc
color.codes = c("#0073C2FF", "#EFC000FF", "#868686FF")
A = ggplot(data=subset(model.res, Method != 'GMAC'), aes(x=`Inferred Model`, y=Count, fill=Method)) +
  geom_bar(stat="identity", position = position_dodge(), color = 'black')+
  scale_fill_manual(values = color.codes[-2])+
  scale_color_manual(values = color.codes[-2])+
  theme_hc()+
  theme(legend.position = 'top',legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
  xlab("Model")


B = ggplot(data=subset(model.res2, Method != 'GMAC'), aes(x=`Inferred Model`, y=Count, fill=Method)) +
  geom_bar(stat="identity", position = position_dodge(), color = 'black')+
  scale_fill_manual(values = color.codes[-2])+
  scale_color_manual(values = color.codes[-2])+
  theme_hc()+
  theme(legend.position = 'top',legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
  xlab("Type of Mediation")


mrgn.adj = list()
mrpc.adj=lapply(mrgn.tab$MRPC.Addis.Inferred.Model, get.adj.from.class)
for(i in 1:dim(mrgn.tab)[1]){
  
  #get the adjacency for MRGN and the truth
  mrgn.adj[[i]]=get.adj.from.class(mrgn.tab$MRGN.Inferred.Model.libconf.alpha05[i], 
                                   reg.vec = unlist(reg.res.L$WholeBlood[1:6, i]))
}



mrgn.edge.ind = unlist(lapply(mrgn.adj, ind.med.edge))
mrpc.edge.ind = unlist(lapply(mrpc.adj, ind.med.edge))
gmac.edge.ind = apply(cbind(gmac.tab$cis.sig.at.01.cutoff,
                            gmac.tab$trans.sig.at.01.cutoff),
                      1, ind.gmac)

#gmac.edge.ind.01 = apply(cbind(gmac.tab$cis.sig.at.01.cutoff,
#                            gmac.tab$trans.sig.at.01.cutoff),
#                      1, ind.gmac)

summ.t1t2.mrgn = summary(as.factor(mrgn.edge.ind))
summ.t1t2.mrpc = summary(as.factor(mrpc.edge.ind))
summ.t1t2.gmac = summary(as.factor(gmac.edge.ind))

t1t2.res = cbind.data.frame(Method = c(rep("MRGN",2), rep("GMAC",2), rep("MRPC-ADDIS",2)),
                            Count = c(summ.t1t2.mrgn, summ.t1t2.gmac, summ.t1t2.mrpc),
                            `T1 T2 Edge Prediction` = rep(c("Absent", "Present"), 3))

t1t2.res$Method = factor(t1t2.res$Method, levels = c("MRGN", "GMAC", "MRPC-ADDIS"))

C = ggplot(data=t1t2.res, aes(x=`T1 T2 Edge Prediction`, y=Count, fill=Method)) +
  geom_bar(stat="identity", position = position_dodge(), color = 'black')+
  scale_fill_manual(values = color.codes)+
  scale_color_manual(values = color.codes)+
  theme_hc()+
  theme(legend.position = 'top',legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))+
  xlab("T1 - T2 Edge")


D = ggplot(data=mrgn.tab, aes(x= MRGN.libconf.alpha05.number.of.PCs)) +
  geom_histogram( color = "black", fill = '#0073C2FF')+
  theme_hc()+
  xlab("Number of Selected PCs")+
  ylab("Count")+
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))



get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}




legend_obj <- get_legend(C)

E = ggarrange(A+theme(legend.position = 'none'),
              B+theme(legend.position = 'none'),
              C+theme(legend.position = 'none'),
              D,
              labels = c("A", "B", "C", "D"),
              nrow = 2, ncol = 2,
              common.legend = F,
              legend = 'top',
              legend.grob = legend_obj,
              font.label = list(size = 16, face = "bold", color = "black"))+theme(plot.margin = unit(c(0, 0, 0.4, 0), "cm"))
pdf(paste0(path_figs, "MF3_GTEx.model.and.t1t2.edge.bargraphs.pdf"), height = 10,
    width = 12)
plot(E)
dev.off()


# 
# 
# jpeg(paste0(figsave_path, "MF3_GTEx.model.and.t1t2.edge.bargraphs.jpeg"),
#      height = 8, width = 12, units = 'in', res = 1000, quality = 90)
# plot(E)
# dev.off()










# data.with.pcs = loadRData(file = paste0("GTEx/results/MRGN/WholeBloodpr.lnc.data.with.pcs.RData"))
# source('adapted_GMAC_func/gmac_one_trio.R')
# #========================================================================
# #looking into GMAC inference differences
# 
# new.data = cbind.data.frame(only5.perm = mrgn.tab$MRGN.Inferred.Model, 
#                             all.perm = perm.all.mods,
#                             gmac.mods = gmac.tab$type.of.med,
#                             maf = unlist(reg.res.L$WholeBlood[13,]),
#                             gmac.ind = gmac.edge.ind, 
#                             only5.mrgn.ind = mrgn.edge.ind, 
#                             mrgn.all.perm = mrgn.perm.all.ind)
# 
# 
# T1=as.data.frame(table(GMAC = new.data$gmac.mods, MRGN.with.perm.all = convert.cats(new.data$all.perm)))
# T2=as.data.frame(table(GMAC.Edge.Prediction = new.data$gmac.ind, MRGN.with.perm.all = convert.cats(new.data$all.perm)))
# T3=as.data.frame(table(GMAC = new.data$gmac.mods, MRGN.with.perm.all = new.data$mrgn.all.perm))
# T4=as.data.frame(table(GMAC = new.data$gmac.ind, MRGN.with.perm.all = new.data$mrgn.all.perm))
# 
# tlist = list(T1, T2, T3, T4)
# library(writexl)
# write_xlsx(tlist, path = 'Manuscript/supplementary_tables/S6_MRGN_GMAC_perm_differences.xlsx')
# 
# x1 = rep(0, dim(new.data)[1])
# ix = which(new.data$only5.mrgn.ind != new.data$mrgn.all.perm)
# x1[ix] = "edge prediction change"
# x1[-ix] = "no change"
# 
# new.data$changes.with.perm = x1
# 
# x2 = rep(0, dim(new.data)[1])
# x3 = rep(0, dim(new.data)[1])
# ix2 = which(new.data$gmac.ind != new.data$mrgn.all.perm)
# x2[ix2] = paste(new.data$all.perm[ix2], new.data$gmac.mods[ix2], sep = " :: ")
# x2[-ix2] = paste(new.data$all.perm[ix2], new.data$gmac.mods[ix2], sep = " :: ")
# x3[ix2] = "T1 - T2 edge prediction different"
# x3[-ix2] = "T1 - T2 edge prediction same "
# 
# new.data$mrgn.gmac.diff.pred.inf = x2
# new.data$mrgn.gmac.diff.pred.label = x3

# ch = summary(as.factor(x2[ix2]))
# nch = summary(as.factor(x2[-ix2]))
# plt.table = cbind.data.frame(count = c(ch, nch), Inference = c(names(ch), names(nch)),  

# class.T1T2.edge.change = function(x){
#   y = NULL
#   if(x[1] == 1){
#     y = "T1 - T2 present"
#   }else{
#     y = "T1 - T2 absent"
#   }
#   
#   return(y)
# }
# 
# 
# edge.pred.changes.perm$`type of change` = apply(edge.pred.changes[,5:6], 1, class.T1T2.edge.change)
# 
# 
# 
# 
# # grid.arrange(arrangeGrob(A + theme(legend.position="none"),
# #                          B + theme(legend.position="none"), nrow=1),
# #              legend_obj,
# #              nrow=2,heights=c(10, 1))





# 
# xdat = na.omit(data.with.pcs[[129]])
# tr129 = as.matrix(xdat[,1:3])
# tr129.pcs = as.matrix(xdat[,-c(1:3)])
# gmacOneTrio(trio = tr129, confounders = tr129.pcs, nperm = 1000, nominal.p = T)














