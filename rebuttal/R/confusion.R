# Apply the paper's code/get_confusion.R to the chunk outputs of one
# dataset x tool. Saves the full confusion object (as in
# analysis/confusionMats/) plus a tidy per-iteration table.
# Rscript rebuttal/R/confusion.R tool=... dataset=... fdr=0.05 repo=. out=... chunks=a.rds,b.rds,...
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "common.R"))
a <- parse_args()
source(file.path(a$repo, "code", "get_confusion.R"))
run <- load_run(strsplit(a$chunks, ",")[[1]])
# guard against summarising a partial or stale set of chunks
if (!is.null(a$n_iter) && !identical(run$iter, seq_len(as.integer(a$n_iter))))
  stop("expected iterations 1..", a$n_iter, " but chunks contain: ", paste(range(run$iter), collapse = "-"),
       " (", length(run$iter), " iterations)")
cm <- get_confusion(input = run, prog = confusion_prog(a$tool), FDR = num(a$fdr))

per_iter <- do.call(rbind, lapply(names(cm$raw_coeff), function(th) {
  x <- cm$raw_coeff[[th]]
  data.frame(dataset = a$dataset, tool = paper_tool_name(a$tool), coeff = as.numeric(th),
             iter = run$iter, x, row.names = NULL)
}))
dir.create(dirname(a$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(cm = cm, per_iter = per_iter, dataset = a$dataset,
             tool = paper_tool_name(a$tool)), a$out)
