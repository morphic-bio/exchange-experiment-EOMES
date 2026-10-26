suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma);library(fgsea);library(jsonlite);library(Rcpp)})
setDTthreads(1)
base <- 'GSEA_expanded_WT_3v3'
out <- base
targets <- c('EOMES','GCM1','ISL1','NKX2-1')
contrasts <- c('JAX','MSK','UCSF','NW','NW_parent')
cfg <- fromJSON(file.path(base,'inputs/config.json'),simplifyVector=FALSE)
sourceCpp(file.path(base,'score.cpp'),cacheDir=file.path(tempdir(),'cpm_score'))
atomic_rds <- function(x, path) {
  saveRDS(x,paste0(path,'.tmp'),compress=FALSE)
  stopifnot(file.rename(paste0(path,'.tmp'),path))
}
rankmat <- function(fc, symbols) {
  z <- apply(fc,2,function(f) {oo<-order(-f,symbols);r<-integer(length(f));r[oo]<-seq_along(f);r})
  if(is.null(dim(z))) z<-matrix(z,ncol=1)
  z
}
collection_of <- function(p) fifelse(grepl('^GOBP_',p),'GO_BP',fifelse(grepl('^HALLMARK_',p),'Hallmark',fifelse(grepl('^REMAP__',p),'ReMap',fifelse(grepl('^CUSTOM_CHIP__',p),'Custom_ChIP','Reactome'))))
calibrate_one <- function(es, null) {
  same <- if(es>=0) null>=0 else null<=0
  exceed <- if(es>=0) null>=es-1e-12 else null<=es+1e-12
  n<-sum(same);stopifnot(n>0)
  list(p=if(es==0)1 else (1+sum(exceed))/(1+n),
       NES=if(es==0)0 else es/mean(abs(null[same])),n=n,exceed=sum(exceed))
}
