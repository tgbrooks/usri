# The ".u" run of code/*.fun.R: tool on the original data with the original
# labels, after set.seed(2025).
# Rscript rebuttal/R/run_unpermuted.R prep=... out=... tool=... model=... gamma=... nsample=128
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "common.R"))
a <- parse_args()
prep <- readRDS(a$prep)
gamma <- if (!is.null(a$gamma) && nzchar(a$gamma)) num(a$gamma) else NULL
set.seed(2025)
u <- run_tool(a$tool, prep$counts, prep$conditions, prep$pair, a$model, gamma,
              as.integer(a$nsample), mc_default = TRUE, unpermuted = TRUE)
dir.create(dirname(a$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(conditions = prep$conditions, u.data = u, session = sessionInfo()), a$out)
