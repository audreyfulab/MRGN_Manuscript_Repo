source("~/Documents/CausalNet/Code/classify_causal_model.R")

# Directory containing the GTEx_wholeblood_posteriorPM_trio_*.txt files.
# Change this if the files live elsewhere.
data_dir <- "."

files <- list.files(
  path = data_dir,
  pattern = "^GTEx_wholeblood_posteriorPM_mrpc_trio_[0-9]+\\.txt$",
  full.names = TRUE
)

# Pull out the trio index from each filename for labeling results.
trio_id <- as.integer(sub(".*_trio_([0-9]+)\\.txt$", "\\1", files))

results <- vector("list", length(files))
names(results) <- trio_id

for (i in seq_along(files)) {
  f <- files[i]
  post.adj <- read.delim(f, header = TRUE, row.names = 1, sep = "\t")
  results[[i]] <- classify_causal_model(post.adj)
}

# Convenience summary table: one row per trio
summary_df <- data.frame(
  trio  = trio_id,
  model = vapply(results, function(x) x$model, character(1)),
  V_T1  = vapply(results, function(x) x$edges["V_T1"], character(1)),
  V_T2  = vapply(results, function(x) x$edges["V_T2"], character(1)),
  T1_T2 = vapply(results, function(x) x$edges["T1_T2"], character(1)),
  do.call(rbind, lapply(results, function(x) x$probs)),
  row.names = NULL
)

print(summary_df)

# Optional: save summary to file
write.table(summary_df, "mediation_trio_baycn_classifications.txt", row.names = FALSE, sep="\t", quote=FALSE)

# generate summary_df_1 for MRGN mediation trios
# generate summary_df_2 for MRPC unique mediation trios
# rbind the two into summary_df
write.table (summary_df, "trios_baycn_classiciations.txt", row.names = FALSE, sep="\t", quote=FALSE)
