

# Convert GMAC inference to binary
#gmac.binary <- gmac.inf
#gmac.binary[gmac.binary != "No.Med"] <- "Mediation"





generate_class_based_metrics <- function(
    mrpc.inf, 
    mrgn.inf, 
    mrgn.inf.alpha01,
    params.model,
    params.model.mc, 
    path_outfile = NULL
) {
  
  #combine gt
  gt.combined = c(params.model, params.model.mc)
  
  # Decompose MRPC inference
  mrpc.inf2 <- unlist(lapply(mrpc.inf, function(x) ifelse(is.null(x$model), 'did not finish', x$model)))
  mrpc.adj <- lapply(mrpc.inf, function(x) x$Adj)
  
  # Get MRGN and true adjacency
  mrgn.adj <- lapply(mrgn.inf, get.adj.from.class)
  mrgn.adj.alpha01 <- lapply(mrgn.inf.alpha01, get.adj.from.class)
  true.adj <- lapply(convert.truth(params.model), get.adj.from.class)
  
  # MRGN (no alpha correction)
  x1.1 <- table(MRGN = convert.cats(mrgn.inf), TRUTH = convert.truth(gt.combined))
  x1.2 <- rbind(x1.1, Total = colSums(x1.1), Recall = round(diag(x1.1)/colSums(x1.1), 4))
  x1.3 <- cbind(x1.2, Total = c(rowSums(x1.2[1:7,]), NA), 
                Precision = c(round(diag(x1.2)/rowSums(x1.2[1:5,]), 4), rep(NA, 3)))
  
  # MRGN alpha = 0.01
  x12.1 <- table(MRGN = convert.cats(mrgn.inf.alpha01), TRUTH = convert.truth(gt.combined))
  x12.2 <- rbind(x12.1, Total = colSums(x12.1), Recall = round(diag(x12.1)/colSums(x12.1), 4))
  x12.3 <- cbind(x12.2, Total = c(rowSums(x12.2[1:7,]), NA), 
                 Precision = c(round(diag(x12.2)/rowSums(x12.2[1:5,]), 4), rep(NA, 3)))
  
  # MRPC
  x2.1 <- table(MRPC = convert.cats(mrpc.inf2), TRUTH = convert.truth(gt.combined))
  x2.2 <- rbind(x2.1, Total = colSums(x2.1), Recall = round(diag(x2.1[2:7,])/colSums(x2.1[2:7,]), 4))
  x2.3 <- cbind(x2.2, Total = c(rowSums(x2.2[1:7,]), NA, NA), 
                Precision = c(NA, round(diag(x2.2[2:7,])/rowSums(x2.2[2:6,]), 4), rep(NA, 3)))
  
  # Combine into a single table
  final.x12 <- cbind(
    `Inference Method` = c(rep(NA, 3), "MRGN", 
                           rep(NA, 8), "MRGN", 
                           rep(NA, 9), "MRPC", rep(NA, 4)),
    `Inference Correction` = c(rep(NA, 3), "None: $\\alpha <0.01$", 
                               rep(NA, 8), "None: $\\alpha <0.01$", 
                               rep(NA, 9), "ADDIS", rep(NA, 4)),
    `Confounder Selection Correction` = c(rep(NA, 3), "FDR $ < 0.05$", 
                                          rep(NA, 8), "None: $\\alpha <0.01$", 
                                          rep(NA, 9), "FDR $< 0.05$", rep(NA, 4)),
    Description = c(rownames(x1.3), NA, rownames(x12.3), NA, rownames(x2.3)),
    rbind.data.frame(x1.3, rep(NA, 7), x12.3, rep(NA, 7), x2.3)
  )
  
  # Optional write
  if (!is.null(path_outfile)) {
    write.csv(final.x12, path_outfile, row.names = FALSE)
  }
  
  return(final.x12)
}


generate_edge_based_metrics <- function(
    mrpc.inf, 
    mrgn.inf, 
    mrgn.inf.alpha01,
    params.model,
    params.model.mc, 
    path_outfile = NULL
) {
  
  
  #combine gt
  m = length(params.model)
  gt.combined = c(params.model, params.model.mc)
  print(paste('this is the length of gt combined', length(gt.combined)))
  
  # Decompose MRPC inference
  mrpc.inf2 <- unlist(lapply(mrpc.inf, function(x) ifelse(is.null(x$model), 'did not finish', x$model)))
  mrpc.adj <- lapply(mrpc.inf, function(x) {
    if (is.null(x$Adj)) {
      'did not finish'
    } else {
      x$Adj
    }
  })
  idx = which(mrpc.inf2 == 'did not finish')
  # Get MRGN and true adjacency
  mrgn.adj <- lapply(mrgn.inf, get.adj.from.class)
  mrgn.adj.alpha01 <- lapply(mrgn.inf.alpha01, get.adj.from.class)
  true.adj <- lapply(convert.truth(gt.combined), get.adj.from.class)
  
  
  n <- length(gt.combined)
  edge.metrics.mrgn <- as.data.frame(matrix(0, nrow = n, ncol = 4))
  colnames(edge.metrics.mrgn) <- c("prec_edge", "recall", "EW_prec_edge", "EW_recall")
  
  edge.metrics.mrgn.alpha01 <- edge.metrics.mrgn
  edge.metrics.mrpc <- edge.metrics.mrgn
  
  # Compute metrics
  for(i in 1:n){
    truth_i <- convert.truth(gt.combined)[i]
    
    edge.metrics.mrgn[i, ] <- c(
      get.metrics(Truth = truth_i, Inferred = mrgn.adj[[i]], get.adj.truth = TRUE),
      get.metrics(Truth = truth_i, Inferred = mrgn.adj[[i]], get.adj.truth = TRUE,
                  weight.edge.directed.present = 1, weight.edge.present.only = 1)
    )
    
    edge.metrics.mrgn.alpha01[i, ] <- c(
      get.metrics(Truth = truth_i, Inferred = mrgn.adj.alpha01[[i]], get.adj.truth = TRUE),
      get.metrics(Truth = truth_i, Inferred = mrgn.adj.alpha01[[i]], get.adj.truth = TRUE,
                  weight.edge.directed.present = 1, weight.edge.present.only = 1)
    )
  }    
  for(i in 1:length(gt.combined[-idx])){
      truth_i = convert.truth(gt.combined[-idx])[i]
      edge.metrics.mrpc[i, ] <- c(
        get.metrics(Truth = truth_i, Inferred = mrpc.adj[-idx][[i]], get.adj.truth = TRUE),
        get.metrics(Truth = truth_i, Inferred = mrpc.adj[-idx][[i]], get.adj.truth = TRUE,
                    weight.edge.directed.present = 1, weight.edge.present.only = 1)
      )
    }
    
  
  
  # Summarize by model class
  mod.class <- unique(convert.truth(params.model))
  mod.class <- sort(mod.class)
  
  summarize_by_class <- function(metrics, model.labels) {
    tab <- sapply(mod.class, function(m) {
      colMeans(metrics[which(convert.truth(model.labels) == m), , drop = FALSE], na.rm = TRUE)
    })
    rownames(tab) <- c("Precision", "Recall", "EW Precision", "EW Recall")
    return(tab)
  }
  
  edge.met.tab.mrgn <- summarize_by_class(edge.metrics.mrgn, params.model)
  edge.met.tab.mrgn.alpha01 <- summarize_by_class(edge.metrics.mrgn.alpha01, params.model)
  edge.met.tab.mrpc <- summarize_by_class(edge.metrics.mrpc, params.model)
  
  # Final table formatting
  edge.met.final <- cbind(
    `Inference Method` = c(NA, "MRGN", rep(NA, 4),
                           "MRGN", rep(NA, 4), 
                           "MRPC", rep(NA, 2)),
    `Inference Correction` = c(NA, "None: $\\alpha < 0.01$", rep(NA, 4),
                               "None: $\\alpha < 0.01$", rep(NA, 4), 
                               "ADDIS", rep(NA, 2)),
    `Confounder Selection Correction` = c(NA, "FDR $<0.05$", rep(NA, 4),
                                          "None: $\\alpha < 0.01$", rep(NA, 4), 
                                          "FDR $<0.05$", rep(NA, 2)),
    Metric = c(rownames(edge.met.tab.mrgn), NA, rownames(edge.met.tab.mrgn.alpha01), NA, rownames(edge.met.tab.mrpc)),
    rbind(round(edge.met.tab.mrgn, 4), rep(NA, 1), 
          round(edge.met.tab.mrgn.alpha01, 4), rep(NA, 1), 
          round(edge.met.tab.mrpc, 4))
  )
  
  # Optionally write to file
  if (!is.null(path_outfile)) {
    write.csv(edge.met.final, path_outfile, row.names = FALSE)
  }
  
  return(list(table=edge.met.final, edge.metrics.mrgn = edge.metrics.mrgn, edge.metrics.mrgn.alpha01 = edge.metrics.mrgn.alpha01,
              edge.metrics.mrpc = edge.metrics.mrpc))
}









generate_t1_t2_results <- function(params,
                                   params.mc,
                                   mrpc.inf,
                                   mrgn.inf,
                                   mrgn.inf.alpha01,
                                   gmac.05,
                                   gmac.01, 
                                   edge.metrics.mrgn, 
                                   edge.metrics.mrgn.alpha01, 
                                   edge.metrics.mrpc,
                                   path_supptabs,
                                   path_tabs) {
  
  
  #combine gt
  m = length(params$model)
  gt.combined = c(params$model, params.mc$model)
  
  # Decompose MRPC inference
  mrpc.inf2 <- unlist(lapply(mrpc.inf, function(x) ifelse(is.null(x$model), 'did not finish', x$model)))
  idx = which(mrpc.inf2 == 'did not finish')
  mrpc.adj <- lapply(mrpc.inf, function(x) x$Adj)[-idx]
  
  # Get MRGN and true adjacency
  mrgn.adj <- lapply(mrgn.inf, get.adj.from.class)
  mrgn.adj.alpha01 <- lapply(mrgn.inf.alpha01, get.adj.from.class)
  true.adj <- lapply(convert.truth(gt.combined), get.adj.from.class)
  
  # Extract indicators
  true.score <- unlist(lapply(true.adj, ind.med.edge))
  mrgn.edge.ind <- unlist(lapply(mrgn.adj, ind.med.edge))
  mrgn.edge.ind.alpha01 <- unlist(lapply(mrgn.adj.alpha01, ind.med.edge))
  mrpc.edge.ind <- unlist(lapply(mrpc.adj, ind.med.edge))
  #gmac.edge.ind <- apply(cbind(gmac.cis$output.table$Cis_Sig, gmac.trans$output.table$Trans_Sig), 1, ind.gmac)
  gmac.edge.ind.at05.sc <- apply(gmac.05, 1, ind.gmac)
  gmac.edge.ind.at01.sc <- apply(gmac.01, 1, ind.gmac)
  
  # Helper to build confusion matrix with PR
  build_table <- function(pred, truth, add_zero = FALSE) {
    t <- if (add_zero) rbind(c(0, 0), table(pred, truth)) else table(pred, truth)
    colnames(t) <- c("T1-T2 Absent", "T1-T2 Present")
    rownames(t) <- c("T1-T2 Pred. Absent", "T1-T2 Pred. Present")
    t2 <- rbind(t, Total = colSums(t), Recall = round(diag(t)/colSums(t), 4))
    t3 <- cbind(t2, Total = c(rowSums(t2)[1:3], NA), Precision = c(round(diag(t)/rowSums(t), 4), rep("", 2)))
    return(list(raw = t, formatted = t3))
  }
  
  # Tables
  t1 <- build_table(mrgn.edge.ind, true.score)
  t12 <- build_table(mrgn.edge.ind.alpha01, true.score)
  t2 <- build_table(mrpc.edge.ind, true.score[-idx])
  #t3 <- build_table(gmac.edge.ind, true.score)
  t31 <- build_table(gmac.edge.ind.at05.sc, true.score)
  t32 <- build_table(gmac.edge.ind.at01.sc, true.score)
  
  # Combine all into final display table
  final.t123 <- cbind(
    `Inference Method` = c(rep(NA, 2), "MRGN", 
                           rep(NA, 4), "MRGN",
                           rep(NA, 4), "MRPC", 
                           #rep(NA, 4), "GMAC",
                           rep(NA, 4), "GMAC", 
                           rep(NA, 4), "GMAC", NA),
    
    `Inference Correction` = c(rep(NA, 2), "None: $\\alpha < 0.01$", 
                               rep(NA, 4), "None: $\\alpha < 0.01$",
                               rep(NA, 4), "ADDIS", 
                               #rep(NA, 4), "FDR $< 0.1$",
                               rep(NA, 4), "None: $\\alpha < 0.05$", 
                               rep(NA, 4), "None: $\\alpha < 0.01$", NA),
    
    `Confounder Selection Correction` = c(rep(NA, 2), "FDR $<0.05$", 
                                          rep(NA, 4), "None: $\\alpha < 0.01$",
                                          rep(NA, 4), "FDR $<0.05$", 
                                          #rep(NA, 4), "FDR $<0.05$",
                                          rep(NA, 4), "FDR $<0.05$", 
                                          rep(NA, 4), "FDR $<0.05$", NA),
    Description = c(rownames(t1$formatted), NA, 
                    rownames(t12$formatted), NA,
                    rownames(t2$formatted), NA, 
                    #rownames(t3$formatted), NA,
                    rownames(t31$formatted), NA, 
                    rownames(t32$formatted)),
    rbind(t1$formatted, rep(NA, 4), 
          t12$formatted, rep(NA, 4), 
          t2$formatted, rep(NA, 4), 
          #t3$formatted, rep(NA, 4), 
          t31$formatted, rep(NA, 4), 
          t32$formatted)
  )
  
  # Save T1-T2 confusion matrix
  write.csv(final.t123, file = paste0(path_supptabs, "T1.T2.edge.results.csv"), row.names = FALSE)
  
  # Compute summary metrics
  # compute_pr_summary <- function(table, metrics = NULL) {
  #   pr <- cbind.data.frame(
  #     Recall = round(diag(table)/colSums(table), 4),
  #     Precision = round(diag(table)/rowSums(table)[1:5], 4)
  #   )
  #   t1t2_edge <- round(c(table[2,2]/sum(table[,2]), table[2,2]/sum(table[2,])), 4)
  #   if (!is.null(metrics)) {
  #     all_edges <- round(colMeans(metrics, na.rm = TRUE)[1:2], 4)
  #   } else {
  #     all_edges <- c(" - ", " - ")
  #   }
  #   return(rbind(pr, t1t2_edge, all_edges))
  # }
  # 
  # MRGN.FDRSC <- compute_pr_summary(t1$raw, edge.metrics.mrgn)
  # MRGN.noFDRSC <- compute_pr_summary(t12$raw, edge.metrics.mrgn.alpha01)
  # MRPC <- compute_pr_summary(t2$raw, edge.metrics.mrpc)
  # GMAC01 <- compute_pr_summary(t32$raw)
  # 
  # MainT.combined.res.final <- rbind(
  #   c(NA, NA, "MRGN+CSnoFDR", NA, "MRGN + CSFDR", NA, "MRPC-ADDIS+CSFDR", NA, "GMAC "),
  #   cbind(
  #     Metric = c("M0", "M1", "M2", "M3", "M4", "T1-T2 Edge", "All Edges"),
  #     MRGN.noFDRSC, MRGN.FDRSC, MRPC, GMAC01
  #   )
  # )
  # 
  # # Save summary table
  # write.csv(MainT.combined.res.final, file = paste0(path_tabs, "MT2_Combined_ALL_METICS_allconfSIMS.csv"), row.names = FALSE)
  # 
  # Return both tables
  #return(list(confusion_table = final.t123,
  #            summary_table = MainT.combined.res.final))
}























