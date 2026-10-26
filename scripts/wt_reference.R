# Extend the frozen EOMES WT assignments to this release's cell-type signatures.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma);library(fgsea);library(jsonlite);library(Rcpp)})
setDTthreads(1)
args <- commandArgs(trailingOnly=FALSE)
root <- dirname(dirname(normalizePath(sub('^--file=','',args[grepl('^--file=',args)]))))
config <- fromJSON(file.path(root,'config.json'))
inp <- file.path(root,'inputs/wt_reference')
out <- file.path(root,'results/tables/wt_reference_scores.rds')
inputs <- c(list.files(inp,full.names=TRUE),
            list.files(file.path(root,'results/tables/ranked_genes'),pattern='tsv$',full.names=TRUE),
            file.path(root,c('config.json','inputs/reference/cell_type_signatures.gmt',
                            'inputs/chip48/top1500_symbols.tsv','inputs/sensitivity/ReMap_EOMES_1446.txt',
                            'scripts/wt_reference.R','provenance/wt_reference/fit.cpp','provenance/wt_reference/score.cpp')))
fingerprint <- setNames(as.character(tools::md5sum(inputs)),substring(inputs,nchar(root)+2))
if (file.exists(out)) {
  cached <- readRDS(out)
  if (identical(cached$fingerprint,fingerprint)) {
    cat('Verified WT score cache: all inputs, rankings, configuration and scoring code match.\n')
    quit(save='no',status=0)
  }
}
sourceCpp(file.path(root,'provenance/wt_reference/score.cpp'))
sourceCpp(file.path(root,'provenance/wt_reference/fit.cpp'))
rankmat <- function(fc,symbols) {
  z <- apply(fc,2,function(f) {oo<-order(-f,symbols);r<-integer(length(f));r[oo]<-seq_along(f);r})
  if(is.null(dim(z))) z<-matrix(z,ncol=1)
  z
}
u <- fread(file.path(inp,'gene_universe.tsv'))
# R's fread can read this small gzip through base R without an R.utils dependency.
cts <- as.matrix(read.delim(gzfile(file.path(inp,'WT_counts.tsv.gz')),row.names=1,check.names=FALSE))
storage.mode(cts)<-'double'
selected <- fread(file.path(inp,'selected_WT_samples.csv'))
a <- fread(file.path(inp,'base_assignments.csv'))
assignments <- fread(file.path(inp,'all_assignments.csv'))
stopifnot(identical(rownames(cts),u$gene_id),identical(colnames(cts),selected$sample_id),
          ncol(cts)==9,nrow(a)==783,nrow(assignments)==1566,!anyNA(cts))
sets <- gmtPathways(file.path(root,'inputs/reference/cell_type_signatures.gmt'))
chip <- fread(file.path(root,'inputs/chip48/top1500_symbols.tsv'))$gene_symbol
remap <- readLines(file.path(root,'inputs/sensitivity/ReMap_EOMES_1446.txt'))
specs <- list()
for (cohort in config$contrasts$id) {
  d <- fread(file.path(root,'results/tables/ranked_genes',paste0(cohort,'.tsv')))
  for (variant in c('all_genes','outside_ChIP48_top1500','outside_ReMap')) {
    dd <- if(variant=='outside_ChIP48_top1500') d[!symbol %in% chip] else if(variant=='outside_ReMap') d[!symbol %in% remap] else d
    setorder(dd,symbol)
    members <- lapply(sets,function(s)which(dd$symbol %in% s))
    keep <- lengths(members)>=config$cell_fate$min_size & lengths(members)<=config$cell_fate$max_size
    members<-members[keep]
    ix <- match(dd$ensembl_gene_id,u$gene_id)
    stopifnot(!anyNA(ix),!anyDuplicated(dd$symbol),all(dd$eligible),all(dd$max_mean_CPM>=1))
    obs <- score_rank_matrix(rankmat(matrix(dd$logFC,ncol=1),dd$symbol),members)[,1]
    specs[[paste(cohort,variant,sep='/')]] <- list(index=ix,symbols=dd$symbol,members=members,
       observed_ES=setNames(obs,names(members)),contrast=cohort,variant=variant)
  }
}
# Reproduce an upstream pathway's scores with its original mappings, to check
# every fitted assignment independently of the new cell-type analysis.
bench <- readRDS(file.path(inp,'source_benchmark_scores.rds'))
bench_specs <- list()
source_targets <- fread(file.path(inp,'source_top500.tsv'))$gene_symbol
for (cn in names(bench)) {
  d <- fread(file.path(inp,paste0('source_observed_ranks_',cn,'.tsv')))
  ix<-match(d$gene_id,u$gene_id);stopifnot(!anyNA(ix))
  bench_specs[[cn]] <- list(index=ix,symbols=d$symbol,members=list(which(d$symbol %in% source_targets)))
}
sixkeys <- unique(a$six_indices)
batches <- split(sixkeys,ceiling(seq_along(sixkeys)/8))
scores <- lapply(specs,function(z)matrix(numeric(),nrow=length(z$members),ncol=0,dimnames=list(names(z$members),NULL)))
benchmark_errors<-numeric();native_errors<-numeric();seen<-character()
parse_ids<-function(x)as.integer(strsplit(x,';',fixed=TRUE)[[1]])
for (batch in seq_along(batches)) {
  indices<-which(a$six_indices %in% batches[[batch]])
  f<-matrix(0,nrow(cts),length(indices));prev<-''
  for(j in seq_along(indices)) {
    i<-indices[j];ii<-parse_ids(a$six_indices[i]);ko<-parse_ids(a$KO_indices[i]);cond<-as.integer(ii %in% ko)
    if(a$six_indices[i]!=prev) {prev<-a$six_indices[i];y<-suppressMessages(calcNormFactors(DGEList(cts[,ii,drop=FALSE])))}
    src<-factor(selected$group_id[ii]);design<-model.matrix(~src+cond)
    stopifnot(qr(design)$rank==ncol(design),6-ncol(design)==a$residual_df[i],a$residual_df[i]>0)
    v<-voom(y,design,plot=FALSE)
    f[,j]<-weighted_condition(v$E,v$weights,as.integer(src),cond)
    if(j==1)native_errors<-c(native_errors,max(abs(f[,j]-lmFit(v,design)$coefficients[,'cond'])))
  }
  stopifnot(all(is.finite(f)),max(native_errors)<1e-8)
  f<-cbind(f,-f)
  ids<-c(paste0(a$base_id[indices],'_forward'),paste0(a$base_id[indices],'_reverse'))
  stopifnot(!any(ids %in% seen));seen<-c(seen,ids)
  for(cn in names(bench_specs)) {
    z<-bench_specs[[cn]]
    es<-score_rank_matrix(rankmat(f[z$index,,drop=FALSE],z$symbols),z$members)[1,]
    benchmark_errors<-c(benchmark_errors,max(abs(es-bench[[cn]][ids])))
  }
  stopifnot(max(benchmark_errors)<1e-10)
  for(key in names(specs)) {
    z<-specs[[key]]
    es<-score_rank_matrix(rankmat(f[z$index,,drop=FALSE],z$symbols),z$members)
    dimnames(es)<-list(names(z$members),ids)
    scores[[key]]<-cbind(scores[[key]],es)
  }
  cat('WT reference batch',batch,'/',length(batches),'; assignments scored:',length(seen),'\n');flush.console()
}
stopifnot(setequal(seen,assignments$assignment_id),length(seen)==1566)
for(key in names(scores)) {
  scores[[key]]<-scores[[key]][,assignments$assignment_id,drop=FALSE]
  stopifnot(all(is.finite(scores[[key]])))
}
qc<-list(status='PASS',WT_samples=9,base_fits=783,assignments=1566,
         max_source_benchmark_ES_error=max(benchmark_errors),max_native_vs_limma_logFC_error=max(native_errors),
         original_q75_fitting_universe=nrow(cts),fixed_comparison_expression_masks=TRUE,
         gene_ID_matching=TRUE,benchmark_scores_checked=5*1566)
saveRDS(list(fingerprint=fingerprint,scores=scores,specs=specs,assignments=assignments,qc=qc),out)
write_json(qc,file.path(root,'qc/wt_reference_validation.json'),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(root,'qc/WT_R_sessionInfo.txt'))
cat('WT reference complete.\n')
