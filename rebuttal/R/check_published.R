# Compare the "published" scheme with the paper's own confusion matrices in
# analysis/confusionMats/ (per iteration, per threshold, TP/FP/TN/FN).
# DESeq2/edgeR/limma should match exactly; ALDEx2/3 match only if the exact
# package version (and RNG use) is reproduced.
# Rscript rebuttal/R/check_published.R repo=. iter_table=... out=...
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "common.R"))
a <- parse_args()
ours <- read.table(a$iter_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
out <- list()
for (key in unique(paste(ours$dataset, ours$tool))) {
  ds <- strsplit(key, " ")[[1]][1]; tl <- strsplit(key, " ")[[1]][2]
  f <- file.path(a$repo, "analysis", "confusionMats", sprintf("cm.%s.%s.Rda", ds, tl))
  if (!file.exists(f)) next
  e <- new.env(); load(f, envir = e); pub <- get(ls(e)[1], envir = e)
  for (th in names(pub$raw_coeff)) {
    p <- pub$raw_coeff[[th]]
    o <- ours[ours$dataset == ds & ours$tool == tl & ours$coeff == as.numeric(th), ]
    o <- o[order(o$iter), ]
    n <- min(nrow(p), nrow(o))
    if (n == 0) next
    d <- abs(as.matrix(p[seq_len(n), c("TP", "FP", "TN", "FN")]) -
             as.matrix(o[seq_len(n), c("TP", "FP", "TN", "FN")]))
    out[[length(out) + 1]] <- data.frame(dataset = ds, tool = tl, coeff = as.numeric(th),
      n_iter_compared = n, n_iter_exact = sum(rowSums(d) == 0), max_abs_diff = max(d),
      published_mean_FDR = mean(p$FDR[seq_len(n)], na.rm = TRUE),
      rerun_mean_FDR = mean(o$FDR[seq_len(n)], na.rm = TRUE))
  }
}
res <- do.call(rbind, out)
write.table(res, a$out, sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(res)) {
  s <- aggregate(cbind(n_iter_compared, n_iter_exact) ~ dataset + tool, res, sum)
  print(s, row.names = FALSE)
}
