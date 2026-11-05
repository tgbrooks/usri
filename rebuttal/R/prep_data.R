# Load + filter one dataset exactly as the paper's scripts do (variant
# "published") or with the data-handling fixes (variant "corrected").
# Rscript rebuttal/R/prep_data.R dataset=immuno variant=published datadir=data out=...
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "common.R"))
a <- parse_args()
d <- load_dataset(a$dataset, a$variant, a$datadir)
tb <- table(d$pair)
message(sprintf("%s (%s): %d genes x %d samples; %d complete pairs, %d unpaired samples",
                a$dataset, a$variant, nrow(d$counts), ncol(d$counts),
                sum(tb == 2), sum(tb == 1)))
dir.create(dirname(a$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(d, a$out)
