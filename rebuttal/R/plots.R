# Plotting functions shared by the pipeline (figures.R, null_counts.R,
# pca_scree.R) and the manuscript (rebuttal/manuscript/manuscript.qmd).
# Depends only on ggplot2, and works from the pipeline's .tsv outputs.
suppressPackageStartupMessages(library(ggplot2))

ds_lab <- c(brca = "TCGA: BRCA", kirc = "TCGA: KIRC", lihc = "TCGA: LIHC",
            luad = "TCGA: LUAD", prad = "TCGA: PRAD", thca = "TCGA: THCA",
            immuno = "Immunotherapy (PD1)", yeast = "Yeast", pseudo = "Single-cell (pseudobulk)")
tool_lab <- c(deseq = "DESeq2", edger = "edgeR (paper)", edger_ql = "edgeR QL + TMM",
              edger_ql_robust = "edgeR QL + TMM, robust", limma = "limma-voom")

read_tsvs <- function(files) {
  files <- files[file.exists(files)]
  do.call(rbind, lapply(files, read.table, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, quote = "", comment.char = ""))
}

## add display columns: family (tool name), gamma, label, dataset_lab
decorate <- function(d) {
  aldex <- grepl("^aldex", d$tool)
  d$family <- ifelse(aldex, sub("\\..*", "", d$tool), d$tool)
  d$family <- ifelse(d$family %in% names(tool_lab), tool_lab[d$family],
                     sub("aldex", "ALDEx", d$family))
  d$gamma <- ifelse(aldex, c("0", "0.1", "0.2", "0.3", "0.4", "0.5")[
    suppressWarnings(as.integer(sub(".*\\.", "", d$tool))) + 1], "naive")
  d$label <- ifelse(d$gamma == "naive", d$family, paste0(d$family, " (g=", d$gamma, ")"))
  d$dataset_lab <- factor(ds_lab[d$dataset], levels = ds_lab)
  d
}

## Per-iteration FDR or TPR at the lowest threshold, boxplots by scheme.
## `iter`: rows of <scheme>/per_iteration.tsv (decorated or not).
per_iteration_plot <- function(iter, metric = c("FDR", "TPR"), target_fdr = 0.05,
                               gammas = c("naive", "0", "0.5")) {
  metric <- match.arg(metric)
  if (is.null(iter$dataset_lab)) iter <- decorate(iter)
  d <- iter[iter$coeff == min(iter$coeff) & iter$gamma %in% gammas, ]
  p <- ggplot(d, aes(label, .data[[metric]], fill = scheme))
  if (metric == "FDR" && !is.null(target_fdr))
    p <- p + geom_hline(yintercept = target_fdr, linetype = "dashed", colour = "grey30")
  p + geom_boxplot(outlier.size = 0.5) + facet_wrap(~ dataset_lab) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = sprintf("%s per iteration (threshold = %s)", metric, min(d$coeff)),
         caption = if (metric == "FDR" && !is.null(target_fdr))
           sprintf("dashed line: target FDR = %s", target_fdr) else NULL)
}

## FDR or TPR against the simulated-effect threshold, one panel per dataset.
line_plot <- function(d, ylab, target_fdr = NULL) {
  p <- ggplot(d, aes(coeff, value, colour = family, group = tool))
  if (!is.null(target_fdr))
    p <- p + geom_hline(yintercept = target_fdr, linetype = "dashed", colour = "grey30")
  p + geom_line(aes(linetype = gamma == "naive")) +
    geom_point(data = d[d$gamma == "naive", ], size = 1) +
    geom_text(data = d[d$gamma != "naive", ], aes(label = sub("0\\.", "", gamma)),
              size = 2, show.legend = FALSE) +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
    facet_wrap(~ dataset_lab) +
    labs(x = "Threshold: model difference between groups", y = ylab, colour = "Tool",
         caption = "ALDEx curves labelled by gamma (0 = 0.001 ... 5 = 0.5)") +
    theme_bw()
}

## Genes called DE with no injected signal. Points: permuted labels (.r), one
## per iteration, one column per scheme within each tool; cross bars: the
## original labels (.u) for that scheme. `res`: rows of null_counts.tsv.
null_counts_plot <- function(res, tools = NULL, schemes = NULL, jitter = 0.08, seed = 1) {
  keep <- !grepl("^aldex", res$tool) | res$tool %in% c("aldex3.0", "aldex3.5")
  res <- res[keep, ]
  if (!is.null(tools)) res <- res[res$tool %in% tools, ]
  if (!is.null(schemes)) res <- res[res$scheme %in% schemes, ]
  res <- decorate(res)
  lab <- unique(res$label[order(res$family, res$gamma)])
  sch <- if (is.null(schemes)) unique(res$scheme) else intersect(schemes, unique(res$scheme))
  w <- 0.8 / length(sch)                       # width of one scheme's column
  res$x <- match(res$label, lab) - 0.4 + w * (match(res$scheme, sch) - 0.5)
  res$scheme <- factor(res$scheme, levels = sch)
  pts <- res[res$type == "permuted (.r)", ]
  set.seed(seed)
  pts$xj <- pts$x + runif(nrow(pts), -jitter, jitter) * w / 0.2
  u <- res[res$type == "original labels (.u)", ]
  ggplot() +
    geom_point(data = pts, aes(xj, n_de + 1, colour = scheme), size = 0.9, alpha = 0.6) +
    geom_segment(data = u, aes(x = x - 0.45 * w, xend = x + 0.45 * w,
                               y = n_de + 1, yend = n_de + 1, colour = scheme),
                 linewidth = 1.1) +
    scale_x_continuous(breaks = seq_along(lab), labels = lab) +
    scale_y_log10() + facet_wrap(~ dataset_lab, scales = "free_y") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor.x = element_blank()) +
    labs(x = NULL, y = "genes with BH-adjusted P < 0.05 (+1, log scale)", colour = "scheme",
         caption = "points: permuted labels, one per iteration (all false positives); bars: original labels")
}

## PCA with paired samples joined. `scores`: pca_<dataset>.tsv
## (columns PC1, PC2, condition, pair, ve1, ve2, dataset).
pca_plot <- function(scores, title = NULL) {
  tb <- table(scores$pair)
  paired <- scores[scores$pair %in% names(tb)[tb == 2], ]
  paired <- paired[order(paired$pair), ]
  segs <- if (nrow(paired)) do.call(rbind, lapply(split(paired, paired$pair), function(p)
    data.frame(x = p$PC1[1], y = p$PC2[1], xend = p$PC1[2], yend = p$PC2[2]))) else NULL
  ds <- scores$dataset[1]
  if (is.null(title))
    title <- sprintf("%s%s", if (ds %in% names(ds_lab)) ds_lab[[ds]] else ds,
                     if (is.null(segs)) "" else sprintf(" (%d pairs joined)", nrow(segs)))
  ggplot(scores, aes(PC1, PC2)) +
    { if (!is.null(segs)) geom_segment(data = segs, aes(x = x, y = y, xend = xend, yend = yend),
                                       colour = "grey60", linewidth = 0.3) } +
    geom_point(aes(colour = condition, shape = condition), size = 1.6) +
    labs(x = sprintf("PC1: %.0f%% variance", 100 * scores$ve1[1]),
         y = sprintf("PC2: %.0f%% variance", 100 * scores$ve2[1]), title = title) +
    theme_bw()
}

## Scree plot from scree_<dataset>.tsv. metric: "variance" (absolute),
## "proportion" (of total), or "both" (two facets).
scree_plot <- function(res, metric = c("variance", "proportion", "both"), nshow = 10,
                       title = NULL) {
  metric <- match.arg(metric)
  ds <- res$dataset[1]
  conds <- unique(res$group[res$set == "original group"])
  perms <- unique(res$set[grepl("^permuted", res$set)])
  res$curve <- paste(res$set, res$group, res$iter)
  res$curve_type <- ifelse(res$set == "original group", paste("original:", res$group), res$set)
  long <- rbind(data.frame(res, metric = "variance", value = res$variance),
                data.frame(res, metric = "proportion of total variance", value = res$proportion))
  if (metric == "variance") long <- long[long$metric == "variance", ]
  if (metric == "proportion") long <- long[long$metric != "variance", ]
  long <- long[long$PC <= nshow, ]
  orig <- long$set %in% c("all samples", "original group")
  lv <- c("all samples", paste("original:", conds), perms)
  long$curve_type <- factor(long$curve_type, levels = lv)
  pal <- setNames(c("black", "#D55E00", "#0072B2")[seq_len(1 + length(conds))],
                  lv[seq_len(1 + length(conds))])
  pal[perms] <- c("grey55", "#CC79A7", "#009E73")[seq_along(perms)]
  if (is.null(title)) title <- if (ds %in% names(ds_lab)) ds_lab[[ds]] else ds
  ggplot(long, aes(PC, value, group = curve, colour = curve_type)) +
    geom_line(data = long[!orig, ], alpha = 0.4, linewidth = 0.4) +
    geom_line(data = long[orig, ], linewidth = 1) +
    geom_point(data = long[orig, ], size = 1.6) +
    scale_colour_manual(values = pal, drop = FALSE) +
    scale_x_continuous(breaks = seq_len(nshow)) +
    expand_limits(y = 0) +
    { if (metric == "both") facet_wrap(~ metric, scales = "free_y") } +
    theme_bw() +
    labs(x = "principal component",
         y = switch(metric, variance = "variance", proportion = "proportion of total variance", both = NULL),
         colour = NULL, title = title)
}
