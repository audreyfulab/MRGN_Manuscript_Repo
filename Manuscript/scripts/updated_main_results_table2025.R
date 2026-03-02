# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")

library(kableExtra)
library(knitr)
library(latex2exp)
library(gridExtra)

library("MRGN")
library("qvalue")
library("caret")


source("Manuscript/scripts/loadedResults.R")

mrgn.inf.combined = c(mrgn.inf, mrgn.many.conf.inf)
mrgn.inf.alpha01.combined = c(mrgn.inf.alpha01, mrgn.mc.inf.alpha01)
gt.combined = c(params$model, many.conf.params$model)
gmac.inf.combined = c(master.table$GMAC.inference.at.05.cutoff, 
                      master.table.many.conf$GMAC.inference.at.05.cutoff)

mrpc.inf.all = loadRData("Manuscript/other/results_all_MRPC-ADDIS.RData")



#==========================================
#-----------Class-Based-Metrics------------
#==========================================
#debug(generate_class_based_metrics)
generate_class_based_metrics(mrpc.inf = mrpc.inf.all,
                             mrgn.inf = mrgn.inf.combined,
                             mrgn.inf.alpha01 = mrgn.inf.alpha01.combined,
                             params.model = params$model,
                             params.model.mc = many.conf.params$model,
                             path_outfile = paste0(path_supptabs, 'Results-class-based.csv'))


#==========================================
#-----------Edge-Based-Metrics------------
#==========================================
#debug(generate_edge_based_metrics)
output = generate_edge_based_metrics(mrpc.inf = mrpc.inf.all,
                             mrgn.inf = mrgn.inf.combined,
                             mrgn.inf.alpha01 = mrgn.inf.alpha01.combined,
                             params.model = params$model,
                             params.model.mc = many.conf.params$model,
                             path_outfile = paste0(path_supptabs, 'Results-edge-based.csv'))

edge.metrics.mrgn = output$edge.metrics.mrgn
edge.metrics.mrgn.alpha01 = output$edge.metrics.mrgn.alpha01
edge.metrics.mrpc = output$edge.metrics.mrpc

#==========================================
#-----T1-T2-Edge-Only-Based-Metrics--------
#==========================================
gmac.05.combined = cbind(c(gmac.cis$output.table$Cis_at_05_cutoff,
                             gmac.cis.many.conf$output.table$Cis_at_05_cutoff),
                           c(gmac.trans$output.table$Trans_at_05_cutoff,
                             gmac.trans.many.conf$output.table$Trans_at_05_cutoff))

gmac.01.combined = cbind(c(gmac.cis$output.table$Cis_at_01_cutoff,
                            gmac.cis.many.conf$output.table$Cis_at_01_cutoff),
                          c(gmac.trans$output.table$Trans_at_01_cutoff,
                            gmac.trans.many.conf$output.table$Trans_at_01_cutoff))


get_inference_type <- function(df) {
  # Check that df has exactly two columns
  if (ncol(df) != 2) stop("Input must be a two-column data frame")
  inf = NULL
  for(i in 1:nrow(df)){
    inf[i] <- ifelse(df[i,1] & df[i,2], 'both',
                  ifelse(df[i,1] & !df[i, 2], 'cis.med',
                         ifelse(!df[i, 1] & df[i,2], 'trans.med',
                                'no.med')))
    
  }
  
  return(inf)
  
}

gmac.05.cats = get_inference_type(gmac.05.combined)
gmac.01.cats = get_inference_type(gmac.01.combined)

mrpc.inf.model = unlist(lapply(mrpc.inf.all, function(x) ifelse(is.null(x), 'did not finish', x$model)))
mrpc.inf.time = unlist(lapply(mrpc.inf.all, function(x) ifelse(is.null(x), 'did not finish', round(as.numeric(x$time)*60, 5))))

mrgn.times = loadRData('Manuscript/other/results_all_times_mrgn.RData')
sim.master = cbind.data.frame(rbind(params, many.conf.params), 
                              MRGN.model.CSFDR = mrgn.inf.combined,
                              MRGN.model.CSnoFDR = mrgn.inf.alpha01.combined,
                              MRGN.time.to.compute.sec = mrgn.times,
                              gmac.alpha05.inf = gmac.05.cats,
                              gmac.alpha01.inf = gmac.01.cats,
                              gmac.time.to.compute.sec = gmac.inf.times,
                              MPRC.model = mrpc.inf.model, 
                              MPRC.time.to.compute.sec = mrpc.inf.time)

write.csv(sim.master, paste0(path_supptabs, "ST_all_results_simulation.csv"), row.names = F)

times.table = cbind.data.frame(Method = c('MRGN', 'MRPC-ADDIS', 'GMAC'),
                               `Median Time To Compute A Trio (sec)` = c(median(mrgn.times),
                                                                         median(as.numeric(mrpc.inf.time[mrpc.inf.time != 'did not finish'])),
                                                                         median(gmac.inf.times)),
                               `IQR Time To Compute A Trio (sec)` = c(IQR(mrgn.times),
                                                                      IQR(as.numeric(mrpc.inf.time[mrpc.inf.time != 'did not finish'])),
                                                                      IQR(gmac.inf.times)))

write.csv(times.table, paste0(path_supptabs, "ST_compute_times.csv"))

#debug(generate_t1_t2_results)
generate_t1_t2_results(params = params,
                       params.mc = many.conf.params,
                       mrpc.inf = mrpc.inf.all,
                       mrgn.inf = mrgn.inf.combined,
                       mrgn.inf.alpha01 = mrgn.inf.alpha01.combined,
                       gmac.05 = gmac.05.combined,
                       gmac.01 = gmac.01.combined,
                       edge.metrics.mrgn,
                       edge.metrics.mrgn.alpha01, 
                       edge.metrics.mrpc,
                       path_supptabs,
                       path_tabs)























































































































