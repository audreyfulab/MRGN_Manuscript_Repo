# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

## this script generates data using the simData.from.graph function to simulate one of the 5 model topologies
## with confounder, intermediate and common child variables
## these models are then inferred using MRGN, GMAC, and MRPC with the confounder selection program get.conf()
## the inference scripts are GMAC_sim2.R ; MRGN_sim2.R and MRPC_sim2.R

library("MRGN")
library("mvtnorm")

set.seed(234)
#pre-allocate simulation conditions
model.types=c("model0","model1","model2","model3","model4")
number.of.datasets=1500
#preallocate the sim parameters for later saving::
params=as.data.frame(matrix(0, nrow = number.of.datasets, ncol = 6))
colnames(params)=c("model", "SD", "minor.freq", "b.snp","b.med",
                   "number.of.Uconfs")

#preallocate list for all simulated datasets
sim.datasets=list()

clinical.covs = loadRData("GTEx/data/kclist_top5_tiss.RData")
kc = clinical.covs$WholeBlood


print("Simulating Data...")
#pc.matrix=as.matrix(pc.matrix)
#simulate all datasets
resamp.vec = NULL
pc.name.list = list()
for(j in 1:number.of.datasets){
  init=0
  k=1
  resamples=NULL
  #simulate trio under parameters
  #use while statement to resample until all genotypes are represented
  #simulate parameters
  params$model[j]=sample(model.types, 1)
  #simulate minor allele frequency between 1 and 50%
  params$minor.freq[j]=sample(seq(0.01, 0.5, 0.01),1)
  #simulate snp effects and mediation effects
  params$b.snp[j]=sample(seq(0.5, 1.5, 0.1), 1)
  params$b.med[j]=sample(seq(0.5, 1, 0.1), 1)
  #simulate Std error of residuals (noise) to be between 1/3 and 1 + 1/2 times mediation signal
  params$SD[j]=(sample(seq(0.3, 1.5, 0.1), 1)*params$b.med[j])
  #generate the number of each type of confounding variable so that
  #the max number of confounding variables possible is 33
  #known.confs = round(runif(1, min = 0, max = 3))
  d = round(runif(1, min = 0.5, max = 15.5))
  U = rmvnorm(n = dim(kc)[1], mean = rep(0, d), sigma = diag(d))

  params$number.of.Uconfs[j] = d
  #create the conf mat to be passed to simData.from.graph:
  KU.mat = cbind.data.frame(kc, U)
  colnames(KU.mat) = c(paste0("K", 1:dim(kc)[2]), paste0("U", c(1:d)))
  row.names(KU.mat) = c(1:dim(kc)[1])

  #checkpoint
  print(paste0("printing params for sim trio ", j))
  print(params[j,])

  while(length(init)<=2){




    #simulate data under parameters
    pdf(paste0("Simulation/sim_graphs/simset_",j,".pdf"))
    X = simData.from.graph(model = params$model[j],
                           theta = params$minor.freq[j],
                           b0.1 = 0,
                           b.snp = params$b.snp[j],
                           b.med = params$b.med[j],
                           sd.1 = params$SD[j],
                           conf.num.vec = c(K = 3, U = d,
                                            W = 1,
                                            Z = 1),
                           simulate.confs = FALSE,
                           conf.mat = KU.mat,
                           plot.graph = TRUE,
                           conf.coef.ranges = list(K = c(0, 0),
                                                   U = c(0.15, 0.5),
                                                   W = c(0.15, 0.5),
                                                   Z = c(1, 1.5)))
    dev.off()
    init=unique(X$data$V1)
    resamples[k]=k
    k = k+1
  }

  print(paste0("resampled ", k, " times before 3 genotypes were represented"))

  colnames(X$data)[c(1:3)+3] = colnames(kc)
  colnames(X$data)[7:(dim(X$data)[2]-2)] = paste0(colnames(KU.mat)[-c(1:3)],".",j)
  colnames(X$data)[(dim(X$data)[2]-1):dim(X$data)[2]] = paste0(colnames(X$data)[(dim(X$data)[2]-1):dim(X$data)[2]],".", j)
  sim.datasets[[j]]=X
}

####post-processing####
sim.datasets2 = lapply(sim.datasets, function(x) x$data)

confs.only.list = lapply(sim.datasets2, function(x) x[,-c(1:6)])
conf.mat = do.call("cbind", confs.only.list)
print("Done!...Saving...")

save(sim.datasets, file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types.RData")
save(conf.mat, file = "Simulation/data/mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types_conf_mat.RData")
save(params, file = "Simulation/data/mrpc_v_mrgn_v_gmac_5k_params_all_mods_conf_types.RData")

print("Done!")
