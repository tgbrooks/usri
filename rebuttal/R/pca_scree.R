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
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(here, "common.R"))
source(file.path(here, "plots.R"))
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
scores <- data.frame(dataset = ds, sample = colnames(vs), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                     condition = prep$conditions, pair = prep$pair, ve1 = ve[1], ve2 = ve[2])
write.table(scores, file.path(a$outdir, sprintf("pca_%s.tsv", ds)), sep = "\t",
            quote = FALSE, row.names = FALSE)
ggsave(file.path(a$outdir, sprintf("pca_%s.png", ds)), pca_plot(scores),
       width = 7.5, height = 6, dpi = 150)

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

nshow <- as.integer(if (is.null(a$nshow)) 10 else a$nshow)
ggsave(file.path(a$outdir, sprintf("scree_%s.png", ds)),
       scree_plot(res, "both", nshow) +
         labs(subtitle = sprintf("vst, all %d genes; thin lines: the two permuted groups of iterations 1-%d",
                                 nrow(vs), n_perm)),
       width = 11, height = 5, dpi = 150)

tv <- unique(res[, c("set", "group", "iter", "n", "total_variance")])
print(aggregate(total_variance ~ set + group, tv, mean))
