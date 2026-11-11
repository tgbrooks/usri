# PCA (pairs connected) and scree plots for one dataset.
#
# Expression: DESeq2::vst(blind = TRUE) of the filtered counts.
#  - PCA: top 500 most variable genes, as DESeq2::plotPCA.
#  - Scree: all filtered genes. PCA of the whole dataset does not depend on
#    the labels, so it is shown once. Within-group PCAs are shown for the
#    original groups and for the permuted groups of the pipeline's first
#    n_perm iterations (same seeds, set.seed(20 + i)), under each permutation
#    scheme. Absolute variance per PC is shown as well as the proportion,
#    since the proportion hides differences in total within-group variance.
#
# Rscript rebuttal/R/pca_scree.R prep=... outdir=... n_perm=10 \
#   permutations=thin2group,pairsplit prop_null=0.95 mean=0 sd=2 npc=30
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "common.R"))
a <- parse_args()
prep <- readRDS(a$prep)
ds <- prep$dataset
npc <- as.integer(if (is.null(a$npc)) 30 else a$npc)
n_perm <- as.integer(a$n_perm)
perms <- strsplit(a$permutations, ",")[[1]]
dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)

vs <- assay(vst(DESeqDataSetFromMatrix(round(prep$counts),
                                       data.frame(row.names = colnames(prep$counts), x = rep(1, ncol(prep$counts))),
                                       ~ 1), blind = TRUE))

## ---- PCA with pairs connected --------------------------------------------
top <- head(order(matrixStats::rowVars(vs), decreasing = TRUE), 500)
pc <- prcomp(t(vs[top, ]), center = TRUE)
ve <- pc$sdev^2 / sum(pc$sdev^2)
df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], condition = prep$conditions,
                 pair = prep$pair)
paired <- df[complete_pairs(df$pair), ]
paired <- paired[order(paired$pair), ]
segs <- do.call(rbind, lapply(split(paired, paired$pair), function(p)
  data.frame(x = p$PC1[1], y = p$PC2[1], xend = p$PC1[2], yend = p$PC2[2])))
g <- ggplot(df, aes(PC1, PC2)) +
  { if (!is.null(segs)) geom_segment(data = segs, aes(x = x, y = y, xend = xend, yend = yend),
                                     colour = "grey60", linewidth = 0.3) } +
  geom_point(aes(colour = condition, shape = condition), size = 1.8) +
  labs(x = sprintf("PC1: %.0f%% variance", 100 * ve[1]),
       y = sprintf("PC2: %.0f%% variance", 100 * ve[2]),
       title = sprintf("%s: vst, top 500 genes; lines join paired samples (%d pairs)",
                       ds, if (is.null(segs)) 0 else nrow(segs))) +
  theme_bw()
ggsave(file.path(a$outdir, sprintf("pca_%s.png", ds)), g, width = 7.5, height = 6, dpi = 150)

## ---- Scree ----------------------------------------------------------------
scree <- function(cols, set, group, perm = NA, iter = NA) {
  x <- t(vs[, cols, drop = FALSE])
  v <- prcomp(x, center = TRUE, rank. = npc)$sdev^2
  k <- min(npc, length(v))
  data.frame(dataset = ds, set = set, group = group, permutation = perm, iter = iter,
             n = length(cols), PC = seq_len(k), variance = v[seq_len(k)],
             proportion = v[seq_len(k)] / sum(v), total_variance = sum(v))
}
out <- list(scree(seq_len(ncol(vs)), "all samples", "all"))
for (cnd in unique(prep$conditions))
  out[[length(out) + 1]] <- scree(which(prep$conditions == cnd), "original group", cnd)
for (perm in perms) for (i in seq_len(n_perm)) {
  set.seed(20 + i)
  th <- make_thin(prep$counts, prep$conditions, prep$pair, perm,
                  prop_null = num(a$prop_null), mean = num(a$mean), sd = num(a$sd))
  g01 <- as.vector(th$designmat)
  for (gg in 0:1)
    out[[length(out) + 1]] <- scree(which(g01 == gg), paste("permuted:", perm),
                                    paste0("group ", gg), perm, i)
}
res <- do.call(rbind, out)
write.table(res, file.path(a$outdir, sprintf("scree_%s.tsv", ds)), sep = "\t",
            quote = FALSE, row.names = FALSE)

res$curve <- paste(res$set, res$group, res$iter)
res$set <- factor(res$set, levels = c("all samples", "original group",
                                      paste("permuted:", perms)))
long <- rbind(data.frame(res, metric = "variance", value = res$variance),
              data.frame(res, metric = "proportion of total", value = res$proportion))
orig <- long$set %in% c("all samples", "original group")
g <- ggplot(long, aes(PC, value, group = curve, colour = set)) +
  geom_line(data = long[!orig, ], alpha = 0.35, linewidth = 0.4) +
  geom_line(data = long[orig, ], aes(linetype = group), linewidth = 0.9) +
  scale_y_log10() + facet_wrap(~ metric, scales = "free_y") + theme_bw() +
  labs(x = "principal component", y = NULL, linetype = "original",
       colour = NULL,
       title = sprintf("%s: scree (vst, all %d genes); permuted = groups from iterations 1-%d",
                       ds, nrow(vs), n_perm))
ggsave(file.path(a$outdir, sprintf("scree_%s.png", ds)), g, width = 11, height = 5, dpi = 150)

tv <- unique(res[, c("set", "group", "iter", "n", "total_variance")])
print(aggregate(total_variance ~ set + group, tv, mean))
