# Load + filter one dataset exactly as the paper's scripts do (variant
# "published") or with the data-handling fixes (variant "corrected").
# Rscript rebuttal/R/prep_data.R dataset=immuno variant=published datadir=data out=...
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "common.R"))
a <- parse_args()
# variant may carry a group subset suffix, e.g. "corrected_g1"
m <- regmatches(a$variant, regexec("^(published|corrected)(_g([12]))?$", a$variant))[[1]]
if (!length(m)) stop("bad variant: ", a$variant)
d <- load_dataset(a$dataset, m[2], a$datadir,
                  subset_group = if (nzchar(m[4])) as.integer(m[4]) else NULL)
d$variant <- a$variant
tb <- table(d$pair)
message(sprintf("%s (%s): %d genes x %d samples (%s); %d complete pairs, %d unpaired samples",
                a$dataset, a$variant, nrow(d$counts), ncol(d$counts),
                paste(unique(d$conditions), collapse = "/"), sum(tb == 2), sum(tb == 1)))
dir.create(dirname(a$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(d, a$out)
