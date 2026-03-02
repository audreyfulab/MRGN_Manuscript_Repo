# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

#Simulate All Data and Parameters:
#source("adapted_GMAC_func/GMACpostproc.R")
library("MRGN")
library("mvtnorm")

model.types=c("model0","model1")
number.of.datasets=2000
#preallocate the sim parameters for later saving::
params=as.data.frame(matrix(0, nrow = number.of.datasets, ncol = 7))
colnames(params)=c("model", "SD", "b.intercept", "minor.freq", "b.snp","b.med",
                   "number.of.Uconfs")
sim.datasets=list()



print("Simulating Data...")
resamp.vec = NULL
pc.name.list = list()
d.vec = c(2,4,6,8,10,12,15)
for(i in 1:length(d.vec)){
  U.mat = rmvnorm(n = 350, mean = rep(0, 350), sigma = diag(350))
  colnames(U.mat) = paste("U", c(1:350), sep = ".")
  prop.sampled = runif(1, 0, 0.2)
  confs.for.trios =
  for(j in 1:number.of.datasets){
    init=0
    k=1
    resamples=NULL
    #simulate trio under parameters
    #use while statement to resample until all genotypes are represented
    #simulate parameters
    params$model[j]=sample(model.types, 1)
    #simulate minor allele frequency between 1 and 50%
    #params$minor.freq[j]=sample(seq(0.01, 0.5, 0.01),1)
    params$minor.freq[j]=0.1
    #simulate snp effects and mediation effects
    params$b.snp[j]=sample(seq(0.5, 1.5, 0.1), 1)
    params$b.med[j]=sample(seq(0.5, 1, 0.1), 1)
    params$b.intercept[j]=sample(seq(0.5, 1, 0.1), 1)
    #simulate Std error of residuals (noise) to be between 1/3 and 1 + 1/2 times mediation signal
    params$SD[j]=sample(seq(0.5, 1, 0.1), 1)
    params$number.of.Uconfs[j] = d.vec[i]

    #checkpoint
    print(paste0("printing params for sim trio ", j))
    print(params[j,])

    while(length(init)<=2){
      #simulate data under parameters
      X = simData.from.graph(model = params$model[j],
                             theta = params$minor.freq[j],
                             b0.1 = params$b.intercept[j],
                             b.snp = params$b.snp[j],
                             b.med = params$b.med[j],
                             sd.1 = params$SD[j],
                             conf.num.vec = c(K = 0, U = params$number.of.Uconfs[j],
                                              W = 0,
                                              Z = 0),
                             simulate.confs = T,
                             sample.size = 350,
                             #conf.mat = KU.mat,
                             plot.graph = FALSE,
                             conf.coef.ranges = list(K = c(0, 0),
                                                     U = c(0.15, 0.5),
                                                     W = c(0.15, 0.5),
                                                     Z = c(1, 1.5)))
      init=unique(X$data$V1)
      resamples[k]=k
      k = k+1
    }


    print(paste0("resampled ", k, " times before 3 genotypes were represented"))
    colnames(X$data)[-c(1:3)] = paste0(colnames(X$data)[-c(1:3)],".",j)
    sim.datasets[[j]]=X
  }


  ####post-processing####
  sim.datasets2 = lapply(sim.datasets, function(x) x$data)

  confs.only.list = lapply(sim.datasets2, function(x) x[,-c(1:3)])
  #conf.mat = cbind.data.frame(pc.matrix, do.call("cbind", confs.only.list))
  conf.mat = do.call("cbind", confs.only.list)
  print("Done!...Saving...")

  save(sim.datasets, file = paste0("Simulation/GMAC_validation_SIM/gmac_valid_sim_simdata_numconfs_",d.vec[i], ".RData"))
  save(conf.mat, file = paste0("Simulation/GMAC_validation_SIM/gmac_valid_sim__conf_mat_numconfs_",d.vec[i],".RData"))
  save(params, file =  paste0("Simulation/GMAC_validation_SIM/gmac_valid_sim_parameters_numcomfs_",d.vec[i],".RData"))

  print("Done!")
}
