
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

path_supptabs = "Manuscript/supplementary_tables/"
path_suppfigs = "Manuscript/supplementary_figures/"
path_figs = "Manuscript/figures/"
path_tabs = "Manuscript/tables/"

source("Manuscript/scripts/loadedResults.R")
source("Manuscript/scripts/create_main_figs_helpers.R")
reso = 2000



mrpc.inf.w.f2=unlist(lapply(mrpc.inf.all, function(x) x$model))
mrpc.w.f.adj=lapply(mrpc.inf.all, function(x) x$Adj)
mrgn.w.f.adj=lapply(mrgn.inf.combined, get.adj.from.class)
mrgn.w.f.libconf.adj=lapply(mrgn.inf.alpha01.combined, get.adj.from.class)
true.adj=lapply(convert.truth(params.combined$model), get.adj.from.class)




true.score = unlist(lapply(true.adj, ind.med.edge))
mrgn.edge.ind.w.f = unlist(lapply(mrgn.w.f.adj, ind.med.edge))
mrgn.edge.ind.w.f.libconf = unlist(lapply(mrgn.w.f.libconf.adj, ind.med.edge))
mrpc.edge.ind.w.f = unlist(lapply(mrpc.w.f.adj, function(x) ifelse(is.null(x), 0, ind.med.edge(x))))
gmac.edge.ind = apply(gmac.raw.combined, 1, ind.gmac)

mrgn.metrics.t1t2 = cbind.data.frame(Method = rep("MRGN + CSFDR", 1800),
                                     truth = true.score,
                                     inf = mrgn.edge.ind.w.f)
mrgn.metrics.t1t2.libconf = cbind.data.frame(Method = rep("MRGN + CSnoFDR", 1800),
                                             truth = true.score,
                                             inf = mrgn.edge.ind.w.f.libconf)
mrpc.metrics.t1t2 = cbind.data.frame(Method = rep("MRPC-ADDIS + CSFDR", 1800),
                                     truth = true.score,
                                     inf = mrpc.edge.ind.w.f)
gmac.metrics.t1t2 = cbind.data.frame(Method = rep("GMAC", 1800),
                                     truth = true.score,
                                     inf = gmac.edge.ind)
colnames(params.combined)[2:6] =c("Residual SD", "Minor Allele Frequency", "SNP Signal", "T1 - T2 Edge Signal", "Number of Simulated Confounders")
#T1 - T2 edge inf plots first:






library(grid)
upper.lim = 1.1
#Residual SD
At1t2=plot.sim.metrics(metrics.mrgn = mrgn.metrics.t1t2,
                       metrics.mrgn.libconf = mrgn.metrics.t1t2.libconf,
                   metrics.mrpc = mrpc.metrics.t1t2,
                   metrics.gmac = gmac.metrics.t1t2,
                   params = params.combined,
                   which.param = 2,
                   by.class=FALSE,
                   use.metric = "t1t2",
                   which.class=NULL,
                   plot.it=TRUE,
                   spline.it=TRUE,
                   spline.int=10,
                   dodge=0.01,
                   return.means=TRUE,
                   rplot=2,
                   nplot=1,
                   xbrks = seq(0.15, 1.5, 0.25),
                   brks1 = seq(0.2, 1, 0.2),
                   brks2 = seq(0.2, 1, 0.2),
                   lmts1 = c(0.2,upper.lim),
                   lmts2 = c(0.2,upper.lim),
                   include.gmac = T,
                   include.mrpc = T,
                   sp.method = "loess",
                   save.plot=FALSE,
                   remove.ylab = TRUE)
#Minor Allele Freq
Bt1t2=plot.sim.metrics(metrics.mrgn = mrgn.metrics.t1t2,
                       metrics.mrgn.libconf = mrgn.metrics.t1t2.libconf,
                       metrics.mrpc = mrpc.metrics.t1t2,
                       metrics.gmac = gmac.metrics.t1t2,
                       params = params.combined,
                       which.param = 3,
                       by.class=FALSE,
                       use.metric = "t1t2",
                       which.class=NULL,
                       plot.it=TRUE,
                       spline.it=TRUE,
                       spline.int=10,
                       dodge=0.01,
                       return.means=TRUE,
                       rplot=2,
                       nplot=1,
                       xbrks = seq(0, 0.5, 0.1),
                       brks1 = seq(0.4, 1, 0.1),
                       brks2 = seq(0.5, 1, 0.1),
                       lmts1 = c(0.4, upper.lim),
                       lmts2 = c(0.5, upper.lim),
                       include.gmac = T,
                       include.mrpc = T,
                       sp.method = "loess",
                       save.plot=FALSE,
                       remove.ylab = TRUE)

#T1 - T2 Signal
Ct1t2=plot.sim.metrics(metrics.mrgn = mrgn.metrics.t1t2,
                       metrics.mrgn.libconf = mrgn.metrics.t1t2.libconf,
                       metrics.mrpc = mrpc.metrics.t1t2,
                       metrics.gmac = gmac.metrics.t1t2,
                       params = params.combined,
                       which.param = 5,
                       by.class=FALSE,
                       use.metric = "t1t2",
                       which.class=NULL,
                       plot.it=TRUE,
                       spline.it=TRUE,
                       spline.int=20,
                       dodge=0.01,
                       return.means=TRUE,
                       rplot=2,
                       nplot=1,
                       xbrks = seq(0.5, 1, 0.1),
                       brks1 = seq(0.2, 1, 0.2),
                       brks2 = seq(0.2, 1, 0.2),
                       lmts1 = c(0.2,upper.lim),
                       lmts2 = c(0.2,upper.lim),
                       include.gmac = T,
                       include.mrpc = T,
                       sp.method = "glm",
                       save.plot=FALSE)

#Number of Hidden Confounders
Dt1t2=plot.sim.metrics(metrics.mrgn = mrgn.metrics.t1t2,
                       metrics.mrpc = mrpc.metrics.t1t2,
                       metrics.mrgn.libconf = mrgn.metrics.t1t2.libconf,
                       metrics.gmac = gmac.metrics.t1t2,
                       params = params.combined,
                       which.param = 6,
                       by.class=FALSE,
                       use.metric = "t1t2",
                       which.class=NULL,
                       plot.it=TRUE,
                       spline.it=TRUE,
                       spline.int=10,
                       dodge=0.01,
                       return.means=TRUE,
                       rplot=2,
                       nplot=1,
                       xbrks = seq(0, 50, 5),
                       brks1 = seq(0.4, 1, 0.1),
                       brks2 = seq(0.5, 1, 0.1),
                       lmts1 = c(0.4,upper.lim),
                       lmts2 = c(0.5,upper.lim),
                       include.gmac = T,
                       include.mrpc = T,
                       sp.method = "loess",
                       save.plot=FALSE)

ABCD = ggarrange(Dt1t2$plot, Bt1t2$plot, Ct1t2$plot, At1t2$plot,
                 labels = c("A", "B", "C", "D"), align = "hv", legend = "top",
                 common.legend = T, legend.grob = At1t2$legend, 
                 hjust = -3.2,
                 vjust = 1, font.label = list(size = 16, face = "bold", color = "black"))+theme(plot.margin = unit(c(0, 0, 0.45, 0), "cm"))

pdf(paste0(path_figs, "MF1_MRGN.GMAC.MRPC.all.params.t1t2.pdf"), height = 10, width = 12)
plot(ABCD)
dev.off()

# jpeg("Manuscript/figures/MF1_MRGN.GMAC.MRPC.all.params.t1t2.jpeg",
#      height = 8, width = 12, units = 'in', res = reso, quality = 90)
# plot(ABCD)
# dev.off()





#edge - based metrics plot MRGN v MRPC

edge.metrics.mrgn.wf=as.data.frame(matrix(NA, nrow = 1800, ncol = 2))
colnames(edge.metrics.mrgn.wf)=c("prec_edge", "recall")
edge.metrics.mrgn.wf.libconf=as.data.frame(matrix(NA, nrow = 1800, ncol = 2))
colnames(edge.metrics.mrgn.wf.libconf)=c("prec_edge", "recall")
edge.metrics.mrpc.wf=as.data.frame(matrix(NA, nrow = 1800, ncol = 2))
colnames(edge.metrics.mrpc.wf)=c("prec_edge", "recall")

temp.inf.mrgn = NULL
temp.inf.mrgn.libconf = NULL
temp.mods = NULL
temp.mrpc = NULL

for(i in 1:1800){
  if (!is.null(mrpc.inf.all[[i]])){
    edge.metrics.mrgn.wf[i,] = get.metrics(Truth = convert.truth(params.combined$model)[i],
                                           Inferred = mrgn.w.f.adj[[i]],
                                           get.adj.truth = TRUE)
    edge.metrics.mrgn.wf.libconf[i,] = get.metrics(Truth = convert.truth(params.combined$model)[i],
                                                   Inferred = mrgn.w.f.libconf.adj[[i]],
                                                   get.adj.truth = TRUE)
    res = get.metrics(Truth = convert.truth(params.combined$model)[i],
                      Inferred = mrpc.inf.all[[i]]$Adj,
                      get.adj.truth = TRUE)
    #print(res)
    edge.metrics.mrpc.wf[i,] = res
    
    temp.inf.mrgn[i] = mrgn.inf.combined[i]
    temp.inf.mrgn.libconf[i] = mrgn.inf.alpha01.combined[i]
    temp.mods[i] = convert.truth(params.combined$model)[i]
    temp.mrpc[i] = ifelse(any(is.na(res)), NA,  mrpc.inf.all[[i]]$model)
  }
  
}
ln = length(temp.mods[-idx])
idx = attr(na.omit(temp.mrpc), "na.action")
edge.metrics.mrgn.wf = cbind.data.frame(Method = rep("MRGN", ln),
                                        truth = temp.mods[-idx],
                                        inf = temp.inf.mrgn[-idx],
                                        edge.metrics.mrgn.wf[-idx, 1:2])

edge.metrics.mrgn.wf.libconf = cbind.data.frame(Method = rep("MRGN", ln),
                                        truth = temp.mods[-idx],
                                        inf = temp.inf.mrgn.libconf[-idx],
                                        edge.metrics.mrgn.wf.libconf[-idx, 1:2])

edge.metrics.gmac.wf = cbind.data.frame(Method = rep("MRPC", ln),
                                        truth = temp.mods[-idx],
                                        inf = rep(1,ln),
                                        matrix(rep(1, ln*2), nrow = ln, ncol = 2))


edge.metrics.mrpc.wf = cbind.data.frame(Method = rep("MRPC", ln),
                                        truth = temp.mods[-idx],
                                        inf = na.omit(temp.mrpc),
                                        na.omit(edge.metrics.mrpc.wf)[,1:2])


upper.lim2 = 1.01
###########################################################
At1t2.1=plot.sim.metrics(metrics.mrgn = edge.metrics.mrgn.wf,
                       metrics.mrpc = edge.metrics.mrpc.wf,
                       metrics.mrgn.libconf = edge.metrics.mrgn.wf.libconf,
                       metrics.gmac = edge.metrics.gmac.wf,
                       params = params.combined[-idx, ],
                       which.param = 2,
                       by.class=FALSE,
                       use.metric = "edge",
                       which.class=NULL,
                       plot.it=TRUE,
                       spline.it=TRUE,
                       spline.int=10,
                       dodge=0.01,
                       return.means=TRUE,
                       rplot=2,
                       nplot=1,
                       xbrks = seq(0.15, 1.5, 0.25),
                       brks1 = seq(0.75, 1, 0.05),
                       brks2 = seq(0.5, 1, 0.1),
                       lmts1 = c(0.75, upper.lim2),
                       lmts2 = c(0.5, upper.lim2),
                       include.gmac = F,
                       include.mrpc = T,
                       sp.method = "loess",
                       save.plot=FALSE,
                       remove.ylab = TRUE)
#Minor Allele Freq
Bt1t2.1=plot.sim.metrics(metrics.mrgn = edge.metrics.mrgn.wf,
                         metrics.mrgn.libconf = edge.metrics.mrgn.wf.libconf,
                       metrics.mrpc = edge.metrics.mrpc.wf,
                       metrics.gmac = edge.metrics.gmac.wf,
                       params = params.combined[-idx, ],
                       which.param = 3,
                       by.class=FALSE,
                       use.metric = "edge",
                       which.class=NULL,
                       plot.it=TRUE,
                       spline.it=TRUE,
                       spline.int=10,
                       dodge=0.01,
                       return.means=TRUE,
                       rplot=2,
                       nplot=1,
                       xbrks = seq(0, 0.5, 0.1),
                       brks1 = seq(0.75, 1, 0.05),
                       brks2 = seq(0.5, 1, 0.1),
                       lmts1 = c(0.8, upper.lim2),
                       lmts2 = c(0.5, upper.lim2),
                       include.gmac = F,
                       include.mrpc = T,
                       sp.method = "loess",
                       save.plot=FALSE,
                       remove.ylab = TRUE)

#T1 - T2 Signal
Ct1t2.1=plot.sim.metrics(metrics.mrgn = edge.metrics.mrgn.wf,
                         metrics.mrgn.libconf = edge.metrics.mrgn.wf.libconf,
                       metrics.mrpc = edge.metrics.mrpc.wf,
                       metrics.gmac = edge.metrics.gmac.wf,
                       params = params.combined[-idx, ],
                       which.param = 5,
                       by.class=FALSE,
                       use.metric = "edge",
                       which.class=NULL,
                       plot.it=TRUE,
                       spline.it=TRUE,
                       spline.int=20,
                       dodge=0.01,
                       return.means=TRUE,
                       rplot=2,
                       nplot=1,
                       xbrks = seq(0.5, 1, 0.1),
                       brks1 = seq(0.75, 1, 0.05),
                       brks2 = seq(0.5, 1, 0.1),
                       lmts1 = c(0.75, upper.lim2),
                       lmts2 = c(0.5, upper.lim2),
                       include.gmac = F,
                       include.mrpc = T,
                       sp.method = "glm",
                       save.plot=FALSE)

#Number of Hidden Confounders
Dt1t2.1=plot.sim.metrics(metrics.mrgn = edge.metrics.mrgn.wf,
                         metrics.mrgn.libconf = edge.metrics.mrgn.wf.libconf,
                       metrics.mrpc = edge.metrics.mrpc.wf,
                       metrics.gmac = edge.metrics.gmac.wf,
                       params = params.combined[-idx,],
                       which.param = 6,
                       by.class=FALSE,
                       use.metric = "edge",
                       which.class=NULL,
                       plot.it=TRUE,
                       spline.it=TRUE,
                       spline.int=10,
                       dodge=0.01,
                       return.means=TRUE,
                       rplot=2,
                       nplot=1,
                       xbrks = seq(0, 50, 5),
                       brks1 = seq(0.75, 1, 0.05),
                       brks2 = seq(0.5, 1, 0.1),
                       lmts1 = c(0.8, upper.lim2),
                       lmts2 = c(0.5, upper.lim2),
                       include.gmac = F,
                       include.mrpc = T,
                       sp.method = "loess",
                       save.plot=FALSE)

ABCD.1 = ggarrange(Dt1t2.1$plot, Bt1t2.1$plot+rremove('ylab'), Ct1t2.1$plot, At1t2.1$plot+rremove('ylab'), labels = c("A", "B", "C", "D"),
                 align = "hv", legend = "top", common.legend = T, legend.grob = At1t2.1$legend, hjust = -3.2,
                 vjust = 0.55, font.label = list(size = 16, face = "bold", color = "black"))+theme(plot.margin = unit(c(0, 0, 0.45, 0), "cm"))

pdf(paste0(path_suppfigs, "SF_MRGN.MRPC.all.params.edge.based.pdf"), height = 10, width = 12)
plot(ABCD.1)
dev.off()

# jpeg("Manuscript/figures/MF2_MRGN.MRPC.all.params.edge.based.jpeg",
#     height = 8, width = 12, units = 'in', res = reso)
# plot(ABCD.1)
# dev.off()



















































# 
# 
# upper.lim3 = 1
# At1t2.3=plot.sim.metrics(metrics.mrgn = mrgn.many.metrics.t1t2,
#                          metrics.mrgn.libconf = mrgn.many.metrics.t1t2.libconf,
#                          metrics.mrpc = mrpc.many.metrics.t1t2,
#                          metrics.gmac = gmac.many.metrics.t1t2,
#                          params = many.conf.params,
#                          which.param = 2,
#                          by.class=FALSE,
#                          use.metric = "t1t2",
#                          which.class=NULL,
#                          plot.it=TRUE,
#                          spline.it=TRUE,
#                          spline.int=10,
#                          dodge=0.01,
#                          return.means=TRUE,
#                          rplot=2,
#                          nplot=1,
#                          xbrks = seq(0.15, 1.5, 0.25),
#                          brks1 = seq(0, 1, 0.2),
#                          brks2 = seq(0.75, 1, 0.05),
#                          lmts1 = c(0, upper.lim3),
#                          lmts2 = c(0.75, upper.lim3),
#                          include.gmac = T,
#                          include.mrpc = F,
#                          sp.method = "loess",
#                          save.plot=FALSE,
#                          remove.ylab = TRUE)
# #Minor Allele Freq
# Bt1t2.3=plot.sim.metrics(metrics.mrgn = mrgn.many.metrics.t1t2,
#                          metrics.mrgn.libconf = mrgn.many.metrics.t1t2.libconf,
#                          metrics.mrpc = mrpc.many.metrics.t1t2,
#                          metrics.gmac = gmac.many.metrics.t1t2,
#                          params = many.conf.params,
#                          which.param = 3,
#                          by.class=FALSE,
#                          use.metric = "t1t2",
#                          which.class=NULL,
#                          plot.it=TRUE,
#                          spline.it=TRUE,
#                          spline.int=10,
#                          dodge=0.01,
#                          return.means=TRUE,
#                          rplot=2,
#                          nplot=1,
#                          xbrks = seq(0, 0.5, 0.1),
#                          brks1 = seq(0.2, 1, 0.2),
#                          brks2 = seq(0.8, 1, 0.05),
#                          lmts1 = c(0.2, upper.lim2),
#                          lmts2 = c(0.8, upper.lim2),
#                          include.gmac = T,
#                          include.mrpc = F,
#                          sp.method = "loess",
#                          save.plot=FALSE,
#                          remove.ylab = TRUE)
# 
# #T1 - T2 Signal
# Ct1t2.3=plot.sim.metrics(metrics.mrgn = mrgn.many.metrics.t1t2,
#                          metrics.mrgn.libconf = mrgn.many.metrics.t1t2.libconf,
#                          metrics.mrpc = mrpc.many.metrics.t1t2,
#                          metrics.gmac = gmac.many.metrics.t1t2,
#                          params = many.conf.params,
#                          which.param = 5,
#                          by.class=FALSE,
#                          use.metric = "t1t2",
#                          which.class=NULL,
#                          plot.it=TRUE,
#                          spline.it=TRUE,
#                          spline.int=20,
#                          dodge=0.01,
#                          return.means=TRUE,
#                          rplot=2,
#                          nplot=1,
#                          xbrks = seq(0.5, 1, 0.1),
#                          brks1 = seq(0, 1, 0.2),
#                          brks2 = seq(0.75, 1, 0.05),
#                          lmts1 = c(0, upper.lim2),
#                          lmts2 = c(0.75, upper.lim2),
#                          include.gmac = T,
#                          include.mrpc = F,
#                          sp.method = "glm",
#                          save.plot=FALSE)
# 
# #Number of Hidden Confounders
# Dt1t2.3=plot.sim.metrics(metrics.mrgn = mrgn.many.metrics.t1t2,
#                          metrics.mrgn.libconf = mrgn.many.metrics.t1t2.libconf,
#                          metrics.mrpc = mrpc.many.metrics.t1t2,
#                          metrics.gmac = gmac.many.metrics.t1t2,
#                          params = many.conf.params,
#                          which.param = 6,
#                          by.class=FALSE,
#                          use.metric = "t1t2",
#                          which.class=NULL,
#                          plot.it=TRUE,
#                          spline.it=TRUE,
#                          spline.int=10,
#                          dodge=0.01,
#                          return.means=TRUE,
#                          rplot=2,
#                          nplot=1,
#                          xbrks = seq(15, 50, 5),
#                          brks1 = seq(0.2, 1, 0.2),
#                          brks2 = seq(0.8, 1, 0.05),
#                          lmts1 = c(0.2, upper.lim2),
#                          lmts2 = c(0.8, upper.lim2),
#                          include.gmac = T,
#                          include.mrpc = F,
#                          sp.method = "loess",
#                          save.plot=FALSE)
# 
# ABCD.3 = ggarrange(Dt1t2.3$plot, Bt1t2.3$plot, Ct1t2.3$plot, At1t2.3$plot, labels = c("A", "B", "C", "D"),
#                    align = "hv", legend = "top", common.legend = T, legend.grob = At1t2.3$legend, hjust = -3.2,
#                    vjust = 1,hjust - 2,
#                    font.label = list(size = 16, face = "bold", color = "black"))+theme(plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))
# 
# pdf(paste0(path_suppfigs, "SF_MRGN.GMAC.all.params.t1t2.50confs.pdf"), height = 10, width = 12)
# plot(ABCD.3)
# dev.off()

# jpeg("Manuscript/supplementary_figures/SF2_MRGN.GMAC.all.params.t1t2.50confs.jpeg",
#     height = 8, width = 12, units = 'in', res = reso)
# plot(ABCD.3)
# dev.off()


















