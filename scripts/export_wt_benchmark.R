suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
source <- args[1]; dest <- args[2]
pathway <- 'CUSTOM_CHIP__EOMES_48h_top500'
files <- sort(list.files(file.path(source,'EOMES/chunks'),pattern='rds$',full.names=TRUE))
scores <- list()
for (f in files) {
  x <- readRDS(f)
  for (cn in names(x$ES)) {
    scores[[cn]] <- c(scores[[cn]],setNames(as.numeric(x$ES[[cn]][pathway,]),x$assignment_ids))
  }
}
stopifnot(all(lengths(scores)==1566L),all(vapply(scores,function(x)!anyDuplicated(names(x)),logical(1))))
saveRDS(scores,file.path(dest,'source_benchmark_scores.rds'))
d <- fread(file.path(source,'EOMES/full_observed_results.csv'))
fwrite(d[pathway=='CUSTOM_CHIP__EOMES_48h_top500'],file.path(dest,'source_benchmark_results.csv'))
