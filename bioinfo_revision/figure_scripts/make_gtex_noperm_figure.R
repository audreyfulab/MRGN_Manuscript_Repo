# GTEx Whole Blood, MRGN WITHOUT the permutation test.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/figure_scripts/make_gtex_noperm_figure.R
#
# The manuscript Figure 3 (Manuscript/figures/MF3_GTEx.model.and.t1t2.edge.bargraphs.pdf,
# built by Manuscript/scripts/create_GTEx_figs.R) plots the master table column
# MRGN.Inferred.Model.libconf.alpha05. This script rebuilds that figure panel for panel on
# MRGN.Inferred.Model.no.perm -- the same 3,248 trios and the same confounder selection,
# with the permutation test switched off.
#
# MRPC-ADDIS and GMAC have no perm / no-perm variants, so their bars are identical to the
# ones in Figure 3. They are kept so the two figures can be read side by side.
#
# ONE DELIBERATE DEPARTURE FROM create_GTEx_figs.R. Panel C needs a trio adjacency, which
# get.adj.from.class() refines using the six MRGN indicators (b11, b12, b21, b22, V1:T1,
# V1:T2). Figure 3 pairs its libconf.alpha05 model labels with the indicators in
# reg.res.wo.pseudo.list.all.tissues.RData, which come from the PERMUTATION run -- label
# and indicators from two different runs. Here both come from the same no-permutation run:
# no.perm[14, ] is the model label (verified below to be identical to the master table
# column) and no.perm[1:6, ] are its indicators.
#
# Writes fig_gtex_noperm.{pdf,png} to bioinfo_revision/reports/figures/.

library('ggpubr')
library('ggthemes')
library('gridExtra')
library('MRGN')

source('Manuscript/scripts/MRGN_write_up_helper_functions.R')

FIG.DIR <- "bioinfo_revision/reports/figures"
dir.create(FIG.DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------

mrgn.tab = loadRData("GTEx/data/WholeBloodmaster_table.RData")

# 14 x 3248. Rows 1:6 are the indicators b11, b12, b21, b22, V1:T1, V1:T2; row 14 is the
# inferred model. Written by Permutation_test_analysis/; read here for BOTH.
no.perm = loadRData("Permutation_test_analysis/results/no.perm.all.trios.WB.RData")

gmac.tab = loadRData("GTEx/results/GMAC/gmac.results.tables.combined.RData")
gmac.tab = gmac.tab$WholeBlood

# The master table column and row 14 of the no-perm run are the same call. If this ever
# breaks, one of the two files has been regenerated without the other, and panel C would be
# silently pairing a label with indicators from a different run.
stopifnot(all(unlist(no.perm[14, ]) == mrgn.tab$MRGN.Inferred.Model.no.perm))
stopifnot(nrow(mrgn.tab) == ncol(no.perm))
stopifnot(nrow(mrgn.tab) == nrow(gmac.tab))

cat("=== GTEx no-permutation figure ===\n")
cat("  trios:", nrow(mrgn.tab), "\n")

# ---------------------------------------------------------------------------------------
# panel A -- the six model classes
# ---------------------------------------------------------------------------------------

sum.mrgn = summary(as.factor(convert.cats(mrgn.tab$MRGN.Inferred.Model.no.perm)))
sum.mrpc = summary(as.factor(convert.cats(mrgn.tab$MRPC.Addis.Inferred.Model)))
models = names(sum.mrgn)

# GMAC returns a mediation call, not a model label, so it is padded with zeros here and
# filtered out of panels A and B -- exactly as create_GTEx_figs.R does it.
model.res = cbind.data.frame(Method = c(rep("MRGN", 6),
                                        rep("GMAC", 6), rep("MRPC-ADDIS", 6)),
                             Count = c(sum.mrgn,
                                       rep(0, 6), sum.mrpc),
                             `Inferred Model` = rep(models, 3))

# ---------------------------------------------------------------------------------------
# panel B -- which gene the mediation runs through
# ---------------------------------------------------------------------------------------

m1.mrgn = subset(mrgn.tab, MRGN.Inferred.Model.no.perm == "M1.1" | MRGN.Inferred.Model.no.perm == "M1.2")$MRGN.Inferred.Model.no.perm
m1.mrpc = subset(mrgn.tab, MRPC.Addis.Inferred.Model == "M1.1" | MRPC.Addis.Inferred.Model == "M1.2")$MRPC.Addis.Inferred.Model

# M1.1 is mediation through the cis gene, M1.2 through the trans gene. Named by the path
# rather than by the gene, so the direction is readable off the axis.
model.res2 = cbind.data.frame(Method = c(rep("MRGN", 2),
                                         rep("GMAC", 2), rep("MRPC-ADDIS", 2)),
                              Count = c(summary(as.factor(m1.mrgn)),
                                        rep(0, 2), summary(as.factor(m1.mrpc))),
                              `Inferred Model` = c("V->Cis->Trans", "V->Trans->Cis"))

model.res$Method = factor(model.res$Method, levels = c("MRGN", "GMAC", "MRPC-ADDIS"))
model.res2$Method = factor(model.res2$Method, levels = c("MRGN", "GMAC", "MRPC-ADDIS"))

# ---------------------------------------------------------------------------------------
# panel C -- the T1-T2 edge
# ---------------------------------------------------------------------------------------

mrgn.adj = list()
mrpc.adj = lapply(mrgn.tab$MRPC.Addis.Inferred.Model, get.adj.from.class)
for(i in 1:dim(mrgn.tab)[1]){

  # label and indicators both from the no-permutation run; see the header note
  mrgn.adj[[i]] = get.adj.from.class(mrgn.tab$MRGN.Inferred.Model.no.perm[i],
                                     reg.vec = unlist(no.perm[1:6, i]))
}

mrgn.edge.ind = unlist(lapply(mrgn.adj, ind.med.edge))
mrpc.edge.ind = unlist(lapply(mrpc.adj, ind.med.edge))
gmac.edge.ind = apply(cbind(gmac.tab$cis.sig.at.01.cutoff,
                            gmac.tab$trans.sig.at.01.cutoff),
                      1, ind.gmac)

summ.t1t2.mrgn = summary(as.factor(mrgn.edge.ind))
summ.t1t2.mrpc = summary(as.factor(mrpc.edge.ind))
summ.t1t2.gmac = summary(as.factor(gmac.edge.ind))

t1t2.res = cbind.data.frame(Method = c(rep("MRGN",2), rep("GMAC",2), rep("MRPC-ADDIS",2)),
                            Count = c(summ.t1t2.mrgn, summ.t1t2.gmac, summ.t1t2.mrpc),
                            `T1 T2 Edge Prediction` = rep(c("Absent", "Present"), 3))

t1t2.res$Method = factor(t1t2.res$Method, levels = c("MRGN", "GMAC", "MRPC-ADDIS"))

# ---------------------------------------------------------------------------------------
# what the figure is built from, printed so a changed input is visible in the run log
# ---------------------------------------------------------------------------------------

# lettered as the figure displays them; see the PANEL ORDER note at the ggarrange call
cat("  A | selected PCs       : median=", stats::median(mrgn.tab$MRGN.number.of.PCs),
    " range=", min(mrgn.tab$MRGN.number.of.PCs), "-", max(mrgn.tab$MRGN.number.of.PCs),
    "\n", sep = "")
cat("  B | MRGN no.perm models:",
    paste(sprintf("%s=%d", models, sum.mrgn), collapse = " "), "\n")
cat("  C | MRGN T1-T2 edge    : absent=", summ.t1t2.mrgn[1],
    " present=", summ.t1t2.mrgn[2], "\n", sep = "")
cat("  D | MRGN M1 split      : V->Cis->Trans=", summary(as.factor(m1.mrgn))[1],
    " V->Trans->Cis=", summary(as.factor(m1.mrgn))[2], "\n", sep = "")

stopifnot(sum(sum.mrgn) == nrow(mrgn.tab))
stopifnot(sum(summ.t1t2.mrgn) == nrow(mrgn.tab))

# ---------------------------------------------------------------------------------------
# the panels -- style carried over verbatim from create_GTEx_figs.R
# ---------------------------------------------------------------------------------------

# get_legend() and ggarrange() below both call ggplot_build(), which measures text and so
# needs an open device. Under Rscript there is none, so R opens the default pdf() and
# leaves a stray Rplots.pdf in the working directory. Send that to the null device for the
# duration of the build; closed again just before the real output devices are opened.
grDevices::pdf(NULL)

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


# MRGN.number.of.PCs, not MRGN.libconf.alpha05.number.of.PCs: the no-perm call was made on
# the standard confounder set, and panel D has to describe the set panels A-C were run on.
D = ggplot(data=mrgn.tab, aes(x= MRGN.number.of.PCs)) +
  geom_histogram( color = "black", fill = '#0073C2FF')+
  theme_hc()+
  xlab("Number of Selected PCs")+
  ylab("Count")+
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.title.x = element_text(margin = margin(t = 0, r = 35, b = 0, l = 0), size = 18))


# ggpubr::get_legend broke on ggplot2 3.5; this digs the grob out directly. Verified
# against ggplot2 4.0.3, where the guide grob is still named "guide-box".
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


# taken from C, the only panel carrying all three methods
legend_obj <- get_legend(C)

# PANEL ORDER. The grid reads
#
#     selected PCs  |  model classes
#     T1-T2 edge    |  type of mediation
#
# i.e. the plot objects go in as D, A, C, B. The displayed letters stay in reading order
# A-D, following create_main_figs.R:162, which likewise passes its panels out of order and
# keeps labels = c("A","B","C","D"). So displayed "A" is the PC histogram, "B" the model
# classes, "C" the T1-T2 edge, "D" the mediation type.
E = ggarrange(D,
              A+theme(legend.position = 'none'),
              C+theme(legend.position = 'none'),
              B+theme(legend.position = 'none'),
              labels = c("A", "B", "C", "D"),
              nrow = 2, ncol = 2,
              common.legend = F,
              legend = 'top',
              legend.grob = legend_obj,
              font.label = list(size = 16, face = "bold", color = "black"))+theme(plot.margin = unit(c(0, 0, 0.4, 0), "cm"))

# ---------------------------------------------------------------------------------------
# write
# ---------------------------------------------------------------------------------------

invisible(grDevices::dev.off())   # the null device opened above

pdf.path <- file.path(FIG.DIR, "fig_gtex_noperm.pdf")
pdf(pdf.path, height = 10, width = 12)
plot(E)
invisible(dev.off())
cat(sprintf("  wrote %s (12 x 10 in)\n", pdf.path))

png.path <- file.path(FIG.DIR, "fig_gtex_noperm.png")
png(png.path, height = 10, width = 12, units = 'in', res = 300)
plot(E)
invisible(dev.off())
cat(sprintf("  wrote %s (12 x 10 in, 300 dpi)\n", png.path))

cat("\ndone.\n")
