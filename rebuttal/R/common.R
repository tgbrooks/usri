# Shared code for the Fig 1 / Fig 2 re-analysis pipeline.
#
# The tool calls below are copied from the paper's code/*.fun.R files so that
# the "published" scheme reproduces the original results; the only structural
# change is that one loop iteration is run per call so that iterations can be
# spread over cluster jobs. Each iteration still starts with set.seed(20 + i)
# and runs thinning -> ".r" analysis -> ".t" analysis in the original order, so
# random-number streams (which matter for ALDEx2/ALDEx3) are unchanged.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  # (the test harness may pre-attach seqgendiff from source)
  if (!"seqgendiff" %in% sub("^package:", "", search())) library(seqgendiff)
})

## ---------------------------------------------------------------------------
## Argument parsing: Rscript script.R key=value key=value ...
## ---------------------------------------------------------------------------
parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  kv <- strsplit(a, "=", fixed = TRUE)
  out <- lapply(kv, function(x) paste(x[-1], collapse = "="))
  names(out) <- vapply(kv, `[`, "", 1)
  out
}
num <- function(x) as.numeric(x)

## ---------------------------------------------------------------------------
## Data loading: reproduces the loading/filtering in code/<dataset>_*.R,
## reading the local copies in data/ instead of raw.githubusercontent.com.
## Returns list(counts, conditions, pair).
## ---------------------------------------------------------------------------
filter_edger <- function(counts, conditions) {
  y <- DGEList(counts = counts, group = factor(conditions))
  keep <- filterByExpr(y)
  y <- y[keep, keep.lib.sizes = FALSE]
  y$counts
}

## subset_group: NULL for all samples, or 1 / 2 to keep only the samples of the
## first / second original condition (sorted alphabetically) for the
## "one group" schemes. Subsetting happens before gene filtering, so genes
## are filtered on the samples that are kept.
load_dataset <- function(dataset, variant = c("published", "corrected"),
                         datadir = "data", subset_group = NULL) {
  variant <- match.arg(variant)
  f <- function(x) file.path(datadir, x)
  keep_group <- function(cond) {
    if (is.null(subset_group)) return(seq_along(cond))
    which(cond == sort(unique(cond))[subset_group])
  }
  tcga <- c(brca = "BRCA", kirc = "KIRC", lihc = "LIHC", luad = "LUAD",
            prad = "PRAD", thca = "THCA")

  if (dataset %in% names(tcga)) {
    t <- tcga[[dataset]]
    raw <- read.table(f(sprintf("TCGA-%s.normal-tumor.pair.rawCount.tsv", t)),
                      header = TRUE, row.names = 1, sep = "\t")
    conditions <- as.vector(unlist(read.table(f(sprintf("TCGA-%s.conditions.tsv", t)),
                                              sep = "\t")))
    # TCGA barcode: TCGA.XX.XXXX.01A -> patient = first 12 characters
    pair <- substr(colnames(raw), 1, 12)
    stype <- substr(colnames(raw), 14, 15)
    stopifnot(all((stype == "01") == (conditions == "tumor")))
    k <- keep_group(conditions)
    raw <- raw[, k]; conditions <- conditions[k]; pair <- pair[k]
    counts <- filter_edger(raw, conditions)

  } else if (dataset == "immuno") {
    m <- read.table(f("imm_metadata.txt"), header = FALSE, row.names = 1, sep = "\t")
    titles <- as.character(unlist(m[1, ]))
    if (variant == "published") {
      # As in code/immuno_*.R. NB: the count file's header is on line 1, so
      # skip = 35 drops the first 35 genes and uses gene 36's counts as the
      # header. The published confusion matrices (20,478 genes) match this.
      raw <- read.table(f("imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv"),
                        header = TRUE, skip = 35, sep = "\t", row.names = 1)
    } else {
      raw <- read.table(f("imm-GSE91061_raw_counts_GRCh38.p13_NCBI.tsv"),
                        header = TRUE, sep = "\t", row.names = 1)
      acc <- as.character(unlist(m[2, ]))
      stopifnot(identical(colnames(raw), acc))
      colnames(raw) <- titles
    }
    conditions <- rep("Pre", 109)
    conditions[grep("_On", m)] <- "On"   # as in the paper's scripts
    stopifnot(all(grepl("_On_", titles) == (conditions == "On")))
    pair <- sub("_.*", "", titles)       # "Pt1_Pre_AD101148-6" -> "Pt1"
    if (variant == "corrected") {
      # Pt109 has two "On" libraries (AE527955-6 / -5): keep only the first
      drop <- which(titles == "Pt109_On_AE527955-5")
      stopifnot(length(drop) == 1)
      raw <- raw[, -drop]; conditions <- conditions[-drop]; pair <- pair[-drop]
    }
    k <- keep_group(conditions)
    raw <- raw[, k]; conditions <- conditions[k]; pair <- pair[k]
    counts <- filter_edger(raw, conditions)

  } else if (dataset == "yeast") {
    raw <- read.table(f("transcriptome.tsv"), header = TRUE, row.names = 1)
    raw <- raw[rownames(raw) != "YOR072W-B", ]
    raw[, c("SNF2.6", "SNF2.13", "SNF2.25", "SNF2.35")] <- NULL
    raw[, c("WT.21", "WT.22", "WT.25", "WT.28", "WT.34", "WT.36")] <- NULL
    conditions <- c(rep("S", 44), rep("W", 42))
    stopifnot(all(startsWith(colnames(raw), "SNF2") == (conditions == "S")))
    pair <- colnames(raw)                # unpaired: every sample its own unit
    k <- keep_group(conditions)
    raw <- raw[, k]; conditions <- conditions[k]; pair <- pair[k]
    counts <- filter_edger(raw, conditions)

  } else if (dataset == "pseudo") {
    raw <- read.table(f("pseudobulk_counts.txt"), header = TRUE, sep = "\t", row.names = 1)
    conditions <- c(rep("CRTL", 103), rep("STIM", 103))   # sic, as in paper
    stopifnot(all(startsWith(colnames(raw), "CTRL") == (conditions == "CRTL")))
    k <- keep_group(conditions)
    raw <- raw[, k]; conditions <- conditions[k]
    raw <- raw[-which(rowSums(raw) == 0), ]
    raw <- as.matrix(raw)
    set.seed(123)
    counts <- select_counts(mat = raw, nsamp = ncol(raw), ngene = 8500,
                            gselect = "max", filter_first = FALSE)
    stopifnot(identical(colnames(counts), colnames(raw)))
    # 8 donors x 13 cell types x CTRL/STIM: pair = donor x cell type
    pair <- sub("^(CTRL|STIM)_", "", colnames(counts))

  } else stop("unknown dataset: ", dataset)

  stopifnot(ncol(counts) == length(conditions), length(pair) == length(conditions))
  list(dataset = dataset, variant = variant, counts = counts,
       conditions = conditions, pair = pair)
}

## ---------------------------------------------------------------------------
## Group assignment that respects pairing.
##   pairsplit:          coin flip per pair; singletons split 50/50
##   pairsplit_balanced: exactly half the pairs flipped (orientation defined by
##                       original condition); singletons split 50/50 within
##                       each original condition
## ---------------------------------------------------------------------------
assign_groups <- function(conditions, pair, balanced = FALSE) {
  n <- length(pair)
  g <- integer(n)
  sizes <- table(pair)
  if (any(sizes > 2)) stop("pairing units with more than 2 samples: ",
                           paste(names(sizes)[sizes > 2], collapse = ", "))
  pairs <- names(sizes)[sizes == 2]
  singles <- which(pair %in% names(sizes)[sizes == 1])
  lv <- sort(unique(conditions))
  half <- function(m) sample(rep_len(sample(0:1), m))  # balanced, random odd one

  if (length(pairs)) {
    flip <- if (balanced) half(length(pairs)) else sample(0:1, length(pairs), replace = TRUE)
    for (j in seq_along(pairs)) {
      idx <- which(pair == pairs[j])
      ref <- if (length(unique(conditions[idx])) == 2) idx[conditions[idx] == lv[1]] else idx[1]
      g[ref] <- flip[j]
      g[setdiff(idx, ref)] <- 1L - flip[j]
    }
  }
  if (length(singles)) {
    if (balanced) {
      for (cnd in unique(conditions[singles])) {
        s <- singles[conditions[singles] == cnd]
        g[s] <- half(length(s))
      }
    } else {
      g[singles] <- half(length(singles))
    }
  }
  g
}

## ---------------------------------------------------------------------------
## Thinning. "thin2group" is exactly the paper's call. The pair-aware versions
## draw the signal in the same RNG order as seqgendiff::thin_2group (effect
## sizes, then which genes, then design), so for a given seed the same genes
## receive the same effects in every scheme; only the grouping differs.
## ---------------------------------------------------------------------------
make_thin <- function(counts, conditions, pair, permutation,
                      prop_null = 0.95, mean = 0, sd = 2) {
  if (permutation == "thin2group") {
    return(thin_2group(counts, prop_null = prop_null, alpha = 0,
                       signal_fun = stats::rnorm,
                       signal_params = list(mean = mean, sd = sd)))
  }
  ngene <- nrow(counts)
  coef <- matrix(0, nrow = ngene, ncol = 1)
  k <- round(ngene * (1 - prop_null))
  if (k > 0) {
    sig <- c(do.call(stats::rnorm, list(mean = mean, sd = sd, n = k)))
    nn <- sort(sample(seq_len(ngene), size = k))
    coef[nn, 1] <- sig
  }
  g <- assign_groups(conditions, pair,
                     balanced = permutation == "pairsplit_balanced")
  th <- thin_diff(mat = counts, design_fixed = matrix(g, ncol = 1),
                  coef_fixed = coef, relative = TRUE)
  colnames(th$designmat) <- "P1"
  colnames(th$coefmat) <- "P1"
  th
}

## ---------------------------------------------------------------------------
## Tools. `grp` is the group vector (0/1 numeric for permuted runs, the
## original character labels for the unpermuted ".u" run). `pair` is used only
## when model == "paired" (then only complete pairs are analysed).
## Output tables have the same columns as in the paper so that
## code/get_confusion.R column indices still apply.
## ---------------------------------------------------------------------------
tool_family <- function(tool) sub("_[0-9]+$", "", tool)

complete_pairs <- function(pair) {
  tb <- table(pair)
  pair %in% names(tb)[tb == 2]
}

run_tool <- function(tool, counts, grp, pair = NULL, model = "unpaired",
                     gamma = NULL, nsample = 128, mc_default = FALSE,
                     unpermuted = FALSE) {
  fam <- tool_family(tool)
  paired <- model == "paired"
  if (paired) {
    keep <- complete_pairs(pair)
    counts <- counts[, keep, drop = FALSE]
    grp <- grp[keep]
    pair <- factor(pair[keep])
  }

  if (fam == "deseq") {
    suppressPackageStartupMessages(library(DESeq2))
    if (unpermuted) {   # des.fun.R, unpermuted block
      cd <- data.frame(conditions = as.factor(grp))
      form <- if (paired) ~ pair + conditions else ~ conditions
    } else {            # des.fun.R loop
      cd <- data.frame(P1 = as.factor(grp))
      form <- if (paired) ~ pair + P1 else ~ P1
    }
    if (paired) cd$pair <- pair
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = cd, design = form)
    dds <- DESeq(dds, quiet = TRUE)
    return(as.data.frame(results(dds)@listData))
  }

  if (fam %in% c("edger", "edger_ql", "edger_ql_robust")) {
    group <- factor(grp)
    design <- if (paired) model.matrix(~ pair + group) else model.matrix(~ group)
    if (fam == "edger") {  # edg.fun.R: raw matrix, offset = log(library size), no TMM
      fit <- glmQLFit(counts, design)  # (edgeR >= 4 => legacy = FALSE)
    } else {               # edgeR v4 QL with TMM (+ robust prior estimation)
      y <- normLibSizes(DGEList(counts = counts))
      fit <- glmQLFit(y, design, legacy = FALSE, robust = fam == "edger_ql_robust")
    }
    qlf <- glmQLFTest(fit, coef = ncol(design))
    tt <- topTags(qlf, n = nrow(counts), adjust.method = "BH", sort.by = "none", p.value = 1)
    return(as.data.frame(tt[[1]]))
  }

  if (fam == "limma") {   # lim.fun.R
    design <- if (paired) model.matrix(~ pair + grp) else model.matrix(~ grp)
    v <- voom(counts, design)
    fit <- eBayes(lmFit(v, design))
    return(topTable(fit, coef = ncol(design), number = nrow(counts), sort.by = "none"))
  }

  if (fam == "aldex2") {  # ald2.fun.R
    suppressPackageStartupMessages(library(ALDEx2))
    if (paired) {
      # ALDEx2's paired test needs group A then group B, in matching pair order
      i0 <- which(grp == sort(unique(grp))[1]); i1 <- which(grp == sort(unique(grp))[2])
      i0 <- i0[order(pair[i0])]; i1 <- i1[order(pair[i1])]
      if (!identical(as.character(pair[i0]), as.character(pair[i1])))
        stop("ALDEx2 paired test needs every pair split across groups (use a pairsplit permutation)")
      o <- c(i0, i1)
      return(ALDEx2::aldex(counts[, o], conditions = grp[o], gamma = gamma,
                           mc.samples = nsample, paired.test = TRUE))
    }
    if (mc_default) return(ALDEx2::aldex(counts, conditions = grp, gamma = gamma))
    return(ALDEx2::aldex(counts, conditions = grp, gamma = gamma, mc.samples = nsample))
  }

  if (fam == "aldex3") {  # ald3.fun.R
    suppressPackageStartupMessages(library(ALDEx3))
    lv <- levels(factor(grp))
    x <- as.numeric(grp == lv[2])           # level 1 -> 0, level 2 -> 1
    if (unpermuted) {
      dd <- data.frame(conditions = x)
      X <- if (paired) formula(~ pair + conditions) else formula(~ conditions)
      term <- "conditions"
    } else {
      dd <- data.frame(conditions_t = x)
      X <- if (paired) formula(~ pair + conditions_t) else formula(~ conditions_t)
      term <- "conditions_t"
    }
    if (paired) dd$pair <- pair
    fit <- ALDEx3::aldex(counts, X, data = dd, nsample = nsample,
                         scale = ALDEx3::clr.sm, gamma = gamma)
    s <- summary(fit)
    if (paired) { s <- s[s$parameter == term, ]; rownames(s) <- NULL }
    stopifnot(nrow(s) == nrow(counts))
    return(s)
  }
  stop("unknown tool: ", tool)
}

## get_confusion() 'prog' argument for each tool
confusion_prog <- function(tool) {
  switch(tool_family(tool), deseq = "deseq", edger = "edger", edger_ql = "edger", edger_ql_robust = "edger",
         limma = "limma", aldex2 = "aldex2", aldex3 = "aldex3")
}
## name used in the paper's confusion-matrix files (aldex3_2 -> aldex3.2)
paper_tool_name <- function(tool) sub("^(aldex[23])_([0-9])$", "\\1.\\2", tool)

## ---------------------------------------------------------------------------
## Rebuild an object in the paper's format (as returned by des.fun() etc.)
## from the chunk files and unpermuted run of one dataset x tool.
## thin.data elements only carry coefmat/designmat (not the thinned counts).
## ---------------------------------------------------------------------------
load_run <- function(chunk_files, unpermuted_file = NULL) {
  ch <- lapply(chunk_files, readRDS)
  iters <- unlist(lapply(ch, `[[`, "iter"))
  o <- order(iters)
  get <- function(field) unlist(lapply(ch, `[[`, field), recursive = FALSE)[o]
  u <- if (!is.null(unpermuted_file) && file.exists(unpermuted_file)) readRDS(unpermuted_file) else NULL
  list(conditions = if (!is.null(u)) u$conditions else NULL,
       iter = iters[o], thin.data = get("thin"),
       u.data = if (!is.null(u)) u$u.data else NULL,
       r.data = get("r.data"), t.data = get("t.data"))
}
