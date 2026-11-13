# Figures analogous to the paper's Fig 1 (FDR) and Fig 2 (TPR) for each
# scheme, plus scheme-comparison figures.
# Rscript rebuttal/R/figures.R outdir=... tidy=a.tsv,b.tsv iter=a.tsv,b.tsv
suppressPackageStartupMessages({ library(ggplot2) })
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(here, "common.R"))
source(file.path(here, "plots.R"))
a <- parse_args()
rd <- function(fs) do.call(rbind, lapply(strsplit(fs, ",")[[1]], read.table,
                                         header = TRUE, sep = "\t", stringsAsFactors = FALSE))
tidy <- rd(a$tidy)
iter <- rd(a$iter)
# optionally overlay the paper's own summary (data/all_TPR_FDR.Rda) as scheme "paper"
if (!is.null(a$paper) && file.exists(a$paper)) {
  e <- new.env(); load(a$paper, envir = e)
  p <- e$plot.df
  p <- p[p$dataset %in% unique(tidy$dataset) & p$tool %in% unique(tidy$tool), ]
  if (nrow(p)) tidy <- rbind(tidy, data.frame(scheme = "paper", p[, c("dataset", "tool", "coeff", "metric", "value")]))
}
dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)

tidy <- decorate(tidy)
iter <- decorate(iter)

for (s in unique(tidy$scheme)) {
  d <- tidy[tidy$scheme == s, ]
  ggsave(file.path(a$outdir, sprintf("fdr_%s.png", s)),
         line_plot(d[d$metric == "cFDR0", ], "FDR", target_fdr = 0.05) + ggtitle(s), width = 11, height = 8, dpi = 150)
  ggsave(file.path(a$outdir, sprintf("tpr_%s.png", s)),
         line_plot(d[d$metric == "cTPR0", ], "TPR") + ggtitle(s), width = 11, height = 8, dpi = 150)
}

if (length(unique(tidy$scheme)) > 1) {
  # scale-naive tools + ALDEx gamma = 0 and 0.5, schemes as line types
  # line type per scheme: fixed for the main ones, the rest from a pool
  fixed <- c(paper = "dotted", published = "solid", pairsplit = "dashed")
  others <- setdiff(unique(tidy$scheme), names(fixed))
  pool <- c("dotdash", "longdash", "twodash", "13", "F1", "4C88C488")
  lty <- c(fixed, setNames(rep_len(pool, length(others)), others))
  keep <- tidy$gamma %in% c("naive", "0", "0.5")
  d <- tidy[keep, ]
  d$label <- ifelse(d$gamma == "naive", d$family, paste0(d$family, " (g=", d$gamma, ")"))
  for (m in c(cFDR0 = "FDR", cTPR0 = "TPR")) {
    met <- names(which(c(cFDR0 = "FDR", cTPR0 = "TPR") == m))
    p <- ggplot(d[d$metric == met, ], aes(coeff, value, colour = label, linetype = scheme)) +
      { if (m == "FDR") geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") } +
      geom_line() + facet_wrap(~ dataset_lab) + theme_bw() +
      scale_linetype_manual(values = lty) +
      labs(x = "Threshold: model difference between groups", y = m, colour = "Tool")
    ggsave(file.path(a$outdir, sprintf("compare_%s.png", tolower(m))), p,
           width = 12, height = 8, dpi = 150)
  }
}

# Per-iteration FDR and TPR at the lowest threshold: shows how much single
# permutations vary (shared across all genes within an iteration)
for (m in c("FDR", "TPR"))
  ggsave(file.path(a$outdir, sprintf("%s_per_iteration.png", tolower(m))),
         per_iteration_plot(iter, m), width = 12, height = 9, dpi = 150)
