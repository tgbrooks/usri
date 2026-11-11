# Figures analogous to the paper's Fig 1 (FDR) and Fig 2 (TPR) for each
# scheme, plus scheme-comparison figures.
# Rscript rebuttal/R/figures.R outdir=... tidy=a.tsv,b.tsv iter=a.tsv,b.tsv
suppressPackageStartupMessages({ library(ggplot2) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "common.R"))
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

ds_lab <- c(brca = "TCGA: BRCA", kirc = "TCGA: KIRC", lihc = "TCGA: LIHC",
            luad = "TCGA: LUAD", prad = "TCGA: PRAD", thca = "TCGA: THCA",
            immuno = "Immunotherapy (PD1)", yeast = "Yeast", pseudo = "Single-cell (pseudobulk)")
tool_lab <- c(deseq = "DESeq2", edger = "edgeR (paper)", edger_ql = "edgeR QL + TMM",
              edger_ql_robust = "edgeR QL + TMM, robust",
              limma = "limma-voom")
decorate <- function(d) {
  d$family <- ifelse(grepl("^aldex", d$tool), sub("\\..*", "", d$tool), d$tool)
  d$family <- ifelse(d$family %in% names(tool_lab), tool_lab[d$family],
                     sub("aldex", "ALDEx", d$family))
  d$gamma <- ifelse(grepl("^aldex", d$tool),
                    c("0", "0.1", "0.2", "0.3", "0.4", "0.5")[as.integer(sub(".*\\.", "", d$tool)) + 1],
                    "naive")
  d$dataset_lab <- factor(ds_lab[d$dataset], levels = ds_lab)
  d
}
tidy <- decorate(tidy)
iter <- decorate(iter)

line_plot <- function(d, ylab) {
  ggplot(d, aes(coeff, value, colour = family, group = tool)) +
    geom_line(aes(linetype = gamma == "naive")) +
    geom_point(data = d[d$gamma == "naive", ], size = 1) +
    geom_text(data = d[d$gamma != "naive", ], aes(label = sub("0\\.", "", gamma)),
              size = 2, show.legend = FALSE) +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
    facet_wrap(~ dataset_lab) +
    labs(x = "Threshold: model difference between groups", y = ylab, colour = "Tool",
         caption = "ALDEx curves labelled by gamma (0 = 0.001 ... 5 = 0.5)") +
    theme_bw()
}

for (s in unique(tidy$scheme)) {
  d <- tidy[tidy$scheme == s, ]
  ggsave(file.path(a$outdir, sprintf("fdr_%s.png", s)),
         line_plot(d[d$metric == "cFDR0", ], "FDR") + ggtitle(s), width = 11, height = 8, dpi = 150)
  ggsave(file.path(a$outdir, sprintf("tpr_%s.png", s)),
         line_plot(d[d$metric == "cTPR0", ], "TPR") + ggtitle(s), width = 11, height = 8, dpi = 150)
}

if (length(unique(tidy$scheme)) > 1) {
  # scale-naive tools + ALDEx gamma = 0 and 0.5, schemes as line types
  keep <- tidy$gamma %in% c("naive", "0", "0.5")
  d <- tidy[keep, ]
  d$label <- ifelse(d$gamma == "naive", d$family, paste0(d$family, " (g=", d$gamma, ")"))
  for (m in c(cFDR0 = "FDR", cTPR0 = "TPR")) {
    met <- names(which(c(cFDR0 = "FDR", cTPR0 = "TPR") == m))
    p <- ggplot(d[d$metric == met, ], aes(coeff, value, colour = label, linetype = scheme)) +
      geom_line() + facet_wrap(~ dataset_lab) + theme_bw() +
      scale_linetype_manual(values = c(paper = "dotted", published = "solid",
        pairsplit = "dashed", pairsplit_balanced = "dotdash", paired_model = "longdash")) +
      labs(x = "Threshold: model difference between groups", y = m, colour = "Tool")
    ggsave(file.path(a$outdir, sprintf("compare_%s.png", tolower(m))), p,
           width = 12, height = 8, dpi = 150)
  }
}

# Per-iteration FDR and TPR at the lowest threshold: shows how much single
# permutations vary (shared across all genes within an iteration)
d <- iter[iter$coeff == min(iter$coeff) & iter$gamma %in% c("naive", "0", "0.5"), ]
d$label <- ifelse(d$gamma == "naive", d$family, paste0(d$family, " (g=", d$gamma, ")"))
for (m in c("FDR", "TPR")) {
  p <- ggplot(d, aes(label, .data[[m]], fill = scheme)) +
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~ dataset_lab) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = sprintf("%s per iteration (threshold = %s)", m, min(iter$coeff)))
  ggsave(file.path(a$outdir, sprintf("%s_per_iteration.png", tolower(m))), p,
         width = 12, height = 9, dpi = 150)
}
