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
    f[,`:=`(contrast=cohort,method=labels[[cohort]],variant=variant,ranked_genes=length(stats))]
    results[[paste(cohort,variant)]] <- f
  }
}
out <- rbindlist(results)
out[,BH_q:=p.adjust(pval,'BH'),by=variant]
stopifnot(all(is.finite(out$NES)),all(is.finite(out$pval)),all(is.finite(out$BH_q)))
# padj is fgsea's within-contrast diagnostic; BH_q is the release's displayed correction.
fwrite(out[variant!='outside_ReMap'],file.path(root,'results/tables/cell_fate_GSEA.tsv'),sep='\t')
fwrite(out[variant!='outside_ChIP48_top1500'],file.path(root,'sensitivity/remap/tables/cell_fate_GSEA.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(root,'qc/R_sessionInfo.txt'))
cat('Cell-fate GSEA complete.\n')
