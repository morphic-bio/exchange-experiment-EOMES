suppressPackageStartupMessages({library(data.table);library(fgsea);library(jsonlite)})
options(warn=1)
args <- commandArgs(trailingOnly=FALSE)
root <- dirname(dirname(normalizePath(sub('^--file=','',args[grepl('^--file=',args)]))))
config <- fromJSON(file.path(root,'config.json'))
cohorts <- config$contrasts$id
labels <- setNames(config$contrasts$label,cohorts)
sets <- gmtPathways(file.path(root,'inputs/reference/cell_type_signatures.gmt'))
stopifnot(length(sets)==config$reference_signature_count)
chip <- fread(file.path(root,'inputs/chip48/top1500_symbols.tsv'))$gene_symbol
remap <- readLines(file.path(root,'inputs/sensitivity/ReMap_EOMES_1446.txt'))
wt <- readRDS(file.path(root,'results/tables/wt_reference_scores.rds'))
# Prevent direct invocation from silently using scores from obsolete inputs.
stopifnot(identical(as.character(tools::md5sum(file.path(root,names(wt$fingerprint)))),unname(wt$fingerprint)))
null_dir <- file.path(root,'results/tables/wt_reference_scores')
dir.create(null_dir,showWarnings=FALSE)
calibrate <- function(es, null) {
  same <- if(es>=0) null>=0 else null<=0
  exceed <- if(es>=0) null>=es-1e-12 else null<=es+1e-12
  n<-sum(same);mu<-mean(abs(null[same]));stopifnot(n>0,mu>0)
  list(pval=if(es==0)1 else (1+sum(exceed))/(1+n),NES=if(es==0)0 else es/mu,
       n_same_sign=n,n_exceed=sum(exceed),mean_abs_same_sign=mu)
}
results <- list()
for (cohort in cohorts) {
  cat('Cell-fate GSEA:',labels[[cohort]],'\n')
  d <- fread(file.path(root,'results/tables/ranked_genes',paste0(cohort,'.tsv')))
  for (variant in c('all_genes','outside_ChIP48_top1500','outside_ReMap')) {
    dd <- if (variant=='outside_ChIP48_top1500') d[!symbol %in% chip] else if (variant=='outside_ReMap') d[!symbol %in% remap] else d
    stats <- sort(setNames(dd$logFC,dd$symbol),decreasing=TRUE)
    set.seed(config$cell_fate$seed)
    f <- fgseaMultilevel(sets,stats,minSize=config$cell_fate$min_size,maxSize=config$cell_fate$max_size,
        eps=config$cell_fate$eps,nPermSimple=config$cell_fate$nPermSimple,gseaParam=config$cell_fate$gseaParam,
        BPPARAM=BiocParallel::SerialParam())
    f[,leadingEdge:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))]
    setnames(f,c('pval','padj','NES','log2err'),c('pval_prerank','padj_prerank','NES_prerank','log2err_prerank'))
    key<-paste(cohort,variant,sep='/')
    null<-wt$scores[[key]]
    stopifnot(setequal(f$pathway,rownames(null)),max(abs(f$ES-wt$specs[[key]]$observed_ES[f$pathway]))<1e-10)
    null<-null[f$pathway,,drop=FALSE]
    for(ref in c('WT','balanced_WT','df_ge2_WT')) {
      use<-switch(ref,WT=rep(TRUE,ncol(null)),balanced_WT=wt$assignments$balanced,df_ge2_WT=wt$assignments$residual_df>=2)
      calibrated<-rbindlist(lapply(seq_len(nrow(f)),function(i)calibrate(f$ES[i],null[i,use])))
      for(field in names(calibrated))set(f,j=paste0(field,'_',ref),value=calibrated[[field]])
      set(f,j=paste0('n_total_',ref),value=sum(use))
    }
    f[,`:=`(pval=pval_WT,NES=NES_WT,reference='WT_3v3')]
    fwrite(data.table(pathway=rownames(null),as.data.table(null)),
           file.path(null_dir,paste0(cohort,'__',variant,'.tsv.gz')),sep='\t',compress='gzip')
    f[,`:=`(contrast=cohort,method=labels[[cohort]],variant=variant,ranked_genes=length(stats))]
    results[[paste(cohort,variant)]] <- f
  }
}
out <- rbindlist(results)
out[,BH_q:=p.adjust(pval,'BH'),by=variant]
out[,BH_q_prerank:=p.adjust(pval_prerank,'BH'),by=variant]
out[,BH_q_balanced_WT:=p.adjust(pval_balanced_WT,'BH'),by=variant]
out[,BH_q_df_ge2_WT:=p.adjust(pval_df_ge2_WT,'BH'),by=variant]
stopifnot(all(is.finite(out$NES)),all(is.finite(out$pval)),all(is.finite(out$BH_q)))
# BH_q is the across-contrast WT correction used in the figures.
# padj_prerank is fgsea's within-contrast diagnostic only.
fwrite(out[variant!='outside_ReMap'],file.path(root,'results/tables/cell_fate_GSEA.tsv'),sep='\t')
fwrite(out[variant!='outside_ChIP48_top1500'],file.path(root,'sensitivity/remap/tables/cell_fate_GSEA.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(root,'qc/R_sessionInfo.txt'))
cat('Cell-fate GSEA complete.\n')
