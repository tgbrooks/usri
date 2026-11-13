# Number of genes called DE (BH-adjusted P < fdr) with NO injected signal:
#   ".r" runs: original data, permuted labels (every call is a false positive)
#   ".u" run : original data, original labels
# One row per scheme x dataset x tool x iteration.
# Rscript rebuttal/R/null_counts.R fdr=0.05 out=... runs=<scheme>/runs/<dataset>/<tool>,...
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(here, "common.R"))
source(file.path(here, "plots.R"))
a <- parse_args()
fdr <- num(a$fdr)
padj_col <- function(tool) switch(tool_family(tool), deseq = 6, edger = 5, edger_ql = 5, edger_ql_robust = 5,
                                  limma = 5, aldex2 = 9, aldex3 = 5)
out <- list()
for (d in strsplit(a$runs, ",")[[1]]) {
  parts <- rev(strsplit(d, "/")[[1]])
  tool <- parts[1]; dataset <- parts[2]; scheme <- parts[4]
  pc <- padj_col(tool)
  ndet <- function(tab) sum(tab[, pc] < fdr, na.rm = TRUE)
  for (f in sort(list.files(d, pattern = "^chunk_.*\\.rds$", full.names = TRUE))) {
    x <- readRDS(f)
    out[[length(out) + 1]] <- data.frame(scheme = scheme, dataset = dataset,
      tool = paper_tool_name(tool), type = "permuted (.r)", iter = x$iter,
      n_de = vapply(x$r.data, ndet, 0), n_genes = nrow(x$r.data[[1]]))
  }
  u <- file.path(d, "unpermuted.rds")
  if (file.exists(u)) {
    x <- readRDS(u)
    out[[length(out) + 1]] <- data.frame(scheme = scheme, dataset = dataset,
      tool = paper_tool_name(tool), type = "original labels (.u)", iter = NA,
      n_de = ndet(x$u.data), n_genes = nrow(x$u.data))
  }
}
res <- do.call(rbind, out)
dir.create(dirname(a$out), recursive = TRUE, showWarnings = FALSE)
write.table(res, a$out, sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(a$plot))
  ggsave(a$plot, null_counts_plot(res), width = 13, height = 9, dpi = 150)
