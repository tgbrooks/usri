# Run iterations [from, to] of the thinning/permutation simulation for one
# dataset x tool x scheme, mirroring the loop body of code/*.fun.R:
#   set.seed(20 + i); thin; tool on original data (.r); tool on thinned (.t)
#
# Rscript rebuttal/R/run_chunk.R prep=... out=... tool=deseq from=1 to=10 \
#   permutation=thin2group model=unpaired prop_null=0.95 mean=0 sd=2 \
#   gamma=0.001 nsample=128 keep_thinned=false

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "common.R"))
a <- parse_args()

prep <- readRDS(a$prep)
iters <- seq(as.integer(a$from), as.integer(a$to))
gamma <- if (!is.null(a$gamma) && nzchar(a$gamma)) num(a$gamma) else NULL
nsample <- as.integer(a$nsample)
keep_thinned <- isTRUE(as.logical(a$keep_thinned))

res <- list(iter = iters, thin = list(), r.data = list(), t.data = list(),
            scheme = list(permutation = a$permutation, model = a$model,
                          variant = prep$variant, tool = a$tool, gamma = gamma),
            session = sessionInfo())

for (k in seq_along(iters)) {
  i <- iters[k]
  message(sprintf("[%s] %s %s iteration %d", format(Sys.time()), prep$dataset, a$tool, i))
  set.seed(20 + i)
  thin <- make_thin(prep$counts, prep$conditions, prep$pair, a$permutation,
                    prop_null = num(a$prop_null), mean = num(a$mean), sd = num(a$sd))
  grp <- as.vector(thin$designmat)

  r <- run_tool(a$tool, prep$counts, grp, prep$pair, a$model, gamma, nsample,
                mc_default = TRUE)
  t <- run_tool(a$tool, thin$mat, grp, prep$pair, a$model, gamma, nsample)

  res$thin[[k]] <- if (keep_thinned) thin else
    list(coefmat = thin$coefmat, designmat = thin$designmat)
  res$r.data[[k]] <- r
  res$t.data[[k]] <- t
}

dir.create(dirname(a$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(res, a$out)
