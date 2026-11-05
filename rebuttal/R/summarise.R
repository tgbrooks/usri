# Collect TPFPR summaries for one scheme into the tidy format used by
# code/0_fig_fdr_tpr_all.R (dataset, tool, coeff, metric, value), plus the
# per-iteration TP/FP/TN/FN table.
# Rscript rebuttal/R/summarise.R scheme=... out_tidy=... out_iter=... files=a.rds,b.rds
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "common.R"))
a <- parse_args()
fs <- strsplit(a$files, ",")[[1]]
tidy <- list(); iter <- list()
for (f in fs) {
  x <- readRDS(f)
  m <- as.data.frame(x$cm$TPFPR)
  for (metric in colnames(m)) {
    tidy[[length(tidy) + 1]] <- data.frame(scheme = a$scheme, dataset = x$dataset,
      tool = x$tool, coeff = as.numeric(rownames(m)), metric = metric, value = m[[metric]])
  }
  iter[[length(iter) + 1]] <- cbind(scheme = a$scheme, x$per_iter)
}
write.table(do.call(rbind, tidy), a$out_tidy, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, iter), a$out_iter, sep = "\t", quote = FALSE, row.names = FALSE)
