# GTEx Whole Blood Figure 4, rebuilt on the inference rerun WITHOUT the known covariates
# (pcr, platform, sex).
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/gtex_rerun/make_gtex_fig_no_known_covs.R --arm=CSq       (Figure 4)
#   Rscript bioinfo_revision/gtex_rerun/make_gtex_fig_no_known_covs.R --arm=CSalpha   (Supp. Figure 9)
#
# Run bioinfo_revision/gtex_rerun/rerun_gtex_no_known_covs.R first. Style is carried over
# from Manuscript/scripts/create_GTEx_figs.R; the panel order and mediation labels are those
# of the manuscript's Figure 4 (as in bioinfo_revision/figure_scripts/make_gtex_noperm_figure.R):
#
#     A  selected PCs (--arm)  |  B  inferred models, MRGN vs MRPC
#     C  T1-T2 edge, 3 methods |  D  type of mediation, MRGN vs MRPC
#
# --arm picks the confounder set MRGN was run on (default CSq):
#   CSq      MRGN.Inferred.Model.CSq.noKC, mrgn_noKC_CSq.RData, PCs = MRGN.number.of.PCs
#   CSalpha  MRGN.Inferred.Model.Csalpha.noKC, mrgn_noKC_CSalpha.RData,
#            PCs = MRGN.libconf.alpha05.number.of.PCs
# Panel C's MRGN indicators always come from the SAME run as its model labels.
# GMAC and MRPC are the same in both figures:
#   GMAC  cis / trans mediation p < 0.01 (GMAC.cis.pval.noKC, GMAC.trans.pval.noKC)
#   MRPC  NOT rerun (most fits at n = 670 with the CS-q confounders did not finish in 10 min);
#         the original MRPC.Addis.Inferred.Model column is plotted, as in create_GTEx_figs.R.
#
# Writes fig_gtex_noKC_MRGN_<arm>.{pdf,png} to bioinfo_revision/gtex_rerun/updated_figures/.

suppressMessages({
  library('ggpubr')
  library('ggthemes')
  library('gridExtra')
  library('MRGN')
})

source('Manuscript/scripts/MRGN_write_up_helper_functions.R')

RES.DIR <- "bioinfo_revision/gtex_rerun/updated_results"
FIG.DIR <- "bioinfo_revision/gtex_rerun/updated_figures"
dir.create(FIG.DIR, recursive = TRUE, showWarnings = FALSE)

args = commandArgs(trailingOnly = TRUE)
ARM  = sub("^--arm=", "", c(grep("^--arm=", args, value = TRUE), "--arm=CSq")[1])
ARMS = list(CSq     = list(model.col = "MRGN.Inferred.Model.CSq.noKC",
                           res.file  = "mrgn_noKC_CSq.RData",
                           pcs.col   = "MRGN.number.of.PCs"),
            CSalpha = list(model.col = "MRGN.Inferred.Model.Csalpha.noKC",
                           res.file  = "mrgn_noKC_CSalpha.RData",
                           pcs.col   = "MRGN.libconf.alpha05.number.of.PCs"))
if (!ARM %in% names(ARMS)) stop("--arm must be one of: ", paste(names(ARMS), collapse = ", "))
arm = ARMS[[ARM]]

# ---------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------

mrgn.tab = read.csv("bioinfo_revision/gtex_rerun/TableS3_GTEx_all_trios_master_updatedOct1.csv",
                    check.names = FALSE)
mrgn.res = loadRData(file.path(RES.DIR, arm$res.file))

# label and indicators must come from the same run, and the table from that file
stopifnot(nrow(mrgn.tab) == 3248, ncol(mrgn.res) == 3248)
stopifnot(all(unlist(mrgn.res["Inferred.Model", ]) == mrgn.tab[[arm$model.col]]))

mrgn.mod = mrgn.tab[[arm$model.col]]
mrgn.tab$n.pcs = mrgn.tab[[arm$pcs.col]]
mrpc.mod = mrgn.tab$MRPC.Addis.Inferred.Model
mrpc.ok  = !is.na(mrpc.mod)
stopifnot(all(mrpc.ok))

cat("=== GTEx figure, no known covariates -- MRGN arm:", ARM, "(MRGN, GMAC rerun; MRPC original) ===\n")
cat("  trios:", nrow(mrgn.tab), "\n")

models = c("M0", "M1", "M2", "M3", "M4", "Other")
count.models = function(x) as.vector(table(factor(convert.cats(x), levels = models)))

# ---------------------------------------------------------------------------------------
# panel B -- the six model classes
# ---------------------------------------------------------------------------------------

sum.mrgn = count.models(mrgn.mod)
sum.mrpc = count.models(mrpc.mod[mrpc.ok])

model.res = cbind.data.frame(Method = c(rep("MRGN", 6),
                                        rep("GMAC", 6), rep("MRPC-ADDIS", 6)),
                             Count = c(sum.mrgn,
                                       rep(0, 6), sum.mrpc),
                             `Inferred Model` = rep(models, 3))

# ---------------------------------------------------------------------------------------
# panel D -- which gene the mediation runs through
# ---------------------------------------------------------------------------------------

m1 = c("M1.1", "M1.2")
m1.mrgn = as.vector(table(factor(mrgn.mod[mrgn.mod %in% m1], levels = m1)))
m1.mrpc = as.vector(table(factor(mrpc.mod[mrpc.mod %in% m1], levels = m1)))

model.res2 = cbind.data.frame(Method = c(rep("MRGN", 2),
                                         rep("GMAC", 2), rep("MRPC-ADDIS", 2)),
                              Count = c(m1.mrgn,
                                        rep(0, 2), m1.mrpc),
                              `Inferred Model` = c("V->Cis->Trans", "V->Trans->Cis"))

model.res$Method = factor(model.res$Method, levels = c("MRGN", "GMAC", "MRPC-ADDIS"))
model.res2$Method = factor(model.res2$Method, levels = c("MRGN", "GMAC", "MRPC-ADDIS"))

# ---------------------------------------------------------------------------------------
# panel C -- the T1-T2 edge
# ---------------------------------------------------------------------------------------

mrgn.adj = lapply(seq_len(nrow(mrgn.tab)), function(i)
  get.adj.from.class(mrgn.mod[i], reg.vec = unlist(mrgn.res[1:6, i])))
mrgn.edge.ind = unlist(lapply(mrgn.adj, ind.med.edge))
mrpc.edge.ind = unlist(lapply(lapply(mrpc.mod, get.adj.from.class), ind.med.edge))
gmac.edge.ind = apply(cbind(mrgn.tab$GMAC.cis.pval.noKC < 0.01,
                            mrgn.tab$GMAC.trans.pval.noKC < 0.01),
                      1, ind.gmac)

count.edge = function(x) as.vector(table(factor(x, levels = c(0, 1))))
summ.t1t2.mrgn = count.edge(mrgn.edge.ind)
summ.t1t2.mrpc = count.edge(mrpc.edge.ind)
summ.t1t2.gmac = count.edge(gmac.edge.ind)

t1t2.res = cbind.data.frame(Method = c(rep("MRGN",2), rep("GMAC",2), rep("MRPC-ADDIS",2)),
                            Count = c(summ.t1t2.mrgn, summ.t1t2.gmac, summ.t1t2.mrpc),
                            `T1 T2 Edge Prediction` = rep(c("Absent", "Present"), 3))

t1t2.res$Method = factor(t1t2.res$Method, levels = c("MRGN", "GMAC", "MRPC-ADDIS"))

# ---------------------------------------------------------------------------------------
# the numbers the manuscript text quotes, printed so they can be updated
# ---------------------------------------------------------------------------------------

cat("  A | selected PCs (", ARM, "): mean=", round(mean(mrgn.tab$n.pcs), 2),
    " sd=", round(sd(mrgn.tab$n.pcs), 2), "\n", sep = "")
cat("  B | MRGN models  :", paste(sprintf("%s=%d", models, sum.mrgn), collapse = " "), "\n")
cat("  B | MRPC models  :", paste(sprintf("%s=%d", models, sum.mrpc), collapse = " "), "\n")
cat("  C | T1-T2 absent/present  MRGN=", paste(summ.t1t2.mrgn, collapse = "/"),
    "  GMAC=", paste(summ.t1t2.gmac, collapse = "/"),
    "  MRPC=", paste(summ.t1t2.mrpc, collapse = "/"), "\n", sep = "")
cat("  D | cis/trans mediation  MRGN=", paste(m1.mrgn, collapse = "/"),
    "  MRPC=", paste(m1.mrpc, collapse = "/"), "\n", sep = "")
cat("  D | mediation trios called by both MRGN and MRPC: cis=",
    sum(mrgn.mod == "M1.1" & mrpc.mod %in% "M1.1"), " trans=",
    sum(mrgn.mod == "M1.2" & mrpc.mod %in% "M1.2"), "\n", sep = "")

stopifnot(sum(sum.mrgn) == nrow(mrgn.tab), sum(summ.t1t2.mrgn) == nrow(mrgn.tab))

# ---------------------------------------------------------------------------------------
# the panels -- style carried over verbatim from create_GTEx_figs.R
# ---------------------------------------------------------------------------------------

# ggplot_build() needs an open device; keep Rscript from leaving a stray Rplots.pdf
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


# confounder counts for the arm panels B-D were inferred on (unchanged by the rerun)
D = ggplot(data=mrgn.tab, aes(x= n.pcs)) +
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

# taken from C, the only panel carrying all three methods
legend_obj <- get_legend(C)

# plot objects go in as D, A, C, B so the displayed A-D read as in the manuscript figure
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

invisible(grDevices::dev.off())   # the null device opened above

# ---------------------------------------------------------------------------------------
# write
# ---------------------------------------------------------------------------------------

pdf.path <- file.path(FIG.DIR, paste0("fig_gtex_noKC_MRGN_", ARM, ".pdf"))
pdf(pdf.path, height = 10, width = 12)
plot(E)
invisible(dev.off())
cat(sprintf("  wrote %s (12 x 10 in)\n", pdf.path))

png.path <- file.path(FIG.DIR, paste0("fig_gtex_noKC_MRGN_", ARM, ".png"))
png(png.path, height = 10, width = 12, units = 'in', res = 300)
plot(E)
invisible(dev.off())
cat(sprintf("  wrote %s (12 x 10 in, 300 dpi)\n", png.path))

cat("\ndone.\n")
