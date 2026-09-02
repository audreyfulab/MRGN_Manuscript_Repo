cat("Script started at:", format(Sys.time()), "\n")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript GTEx_wholeblood_baycn_cluster.R <iterations> <thinTo> <trioIdx> <run_number>")
}
iterations <- as.numeric(args[1])
thinTo <- as.numeric(args[2])
trioIdx <- as.numeric(args[3])
run_number <- as.numeric(args[4])

cat("Iterations:", iterations, "\n")
cat("ThinTo:", thinTo, "\n")
cat("Trio index:", trioIdx, "\n")
cat("Run number:", run_number, "\n")

# Define home directory
home_dir <- "/wsu/home/ht/ht26/ht2699"

# Load libraries
library ("baycn", lib=file.path(home_dir, "Rpackages"))

# load data of trios
load(file.path(home_dir, "GTEx_MRGN/all.data.unqiue.snps.pclrna.only.WholeBlood.RData"))
trios.final <- all.unq.snps.ldf.pl.only
rm (all.unq.snps.ldf.pl.only)

# load known covariates
load(file.path(home_dir, "GTEx_MRGN/kclist_top5_tiss.RData"))
# extract known covariates for whole blood
kc.wb <- kc.list[[1]]

# load gene expression PCs
load(file.path(home_dir, "GTEx_MRGN/PCs.matrix.WholeBlood.RData"))
dim (PCs.matrix.WholeBlood)

# load PCs identified to be confounders for each trio
load(file.path(home_dir, "GTEx_MRGN/List.Match.significant.trios.WholeBlood.RData"))
trios.pc <- List.significant.asso1
rm (List.significant.asso1)

# load trio inference results and find mediation trios
mrgn.gtex <- read.delim(file.path(home_dir, "GTEx_MRGN/TableS2_GTEx_all_trios_master.csv"), header = TRUE, row.names = 1, sep = ",")
id.m1 <- which (mrgn.gtex$MRGN.Inferred.Model.no.perm=="M1.1" | mrgn.gtex$MRGN.Inferred.Model.no.perm=="M1.2")

# compile data for a specific trio
#for (i in beginIdx:endIdx) {
#  print (i)
  # calculate the indices for the trio in trios.final
  trio.id <- (trioIdx-1)*3+1:3
  # extract trio data and merge with known and identified confounders
  trio.data <- data.frame(trios.final[,trio.id], kc.wb, PCs.matrix.WholeBlood[,trios.pc[[trioIdx]]])
  
  # generate the adj matrix for baycn
  trio.am <- matrix(0, nrow = ncol(trio.data), ncol = ncol(trio.data))
  trio.am[1,2] <- 1
  trio.am[1,3] <- 1
  trio.am[2,3] <- 1
  trio.am[3,2] <- 1
  trio.am[4:nrow(trio.am), 2] <- 1
  trio.am[4:nrow(trio.am), 3] <- 1
  
  # run baycn
  trio.baycn <- mhEdge(data = as.matrix(trio.data),
                       adjMatrix = trio.am,
                       prior = c(0.05,
                                 0.05,
                                 0.9),
                       nCPh = 0,
                       nGV = 1,
                       pmr = TRUE,
                       burnIn = 0.2,
                       iterations = iterations,
                       thinTo = thinTo,
                       progress = TRUE)
  
  save(trio.baycn, file=paste0(home_dir, "/GTEx_MRGN/GTEx_wholeblood_baycn_trio_", trioIdx, "_",run_number,".RData"))
#}


cat("Script completed at:", format(Sys.time()), "\n")
