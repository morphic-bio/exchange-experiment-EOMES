# Fit each original WT assignment once; score its coefficients on five fixed masks.
source('GSEA_expanded_WT_3v3/common.R')
args<-commandArgs(TRUE);target<-args[1];worker<-as.integer(args[2]);workers<-as.integer(args[3])
stopifnot(file.exists(file.path(out,'PREPARATION_COMPLETED')))
td<-file.path(out,target);spec<-cfg$targets[[target]]
u<-fread(file.path(base,'inputs',target,'gene_universe.tsv'))
raw<-fread(spec$counts);ids<-raw[[1]];cts<-as.matrix(raw[,-1]);storage.mode(cts)<-'double';rownames(cts)<-ids
meta<-fread(spec$metadata)
if(target=='GCM1') {meta[,sample:=label];idx<-match(colnames(cts),meta$file_id)} else {
  colnames(cts)<-sub('_S[0-9]+_L[0-9]+$','',colnames(cts));meta[,sample:=Name];idx<-match(colnames(cts),meta$Name)
}
stopifnot(!anyNA(idx));meta<-meta[idx];colnames(cts)<-meta$sample
cts<-cts[match(u$gene_id,rownames(cts)),,drop=FALSE]
reps<-which(u$representative);reps<-reps[order(u$symbol[reps])];symbols<-u$symbol[reps]
selected<-fread(file.path(base,'inputs/selected_WT_samples.csv'))[target_gene==target]
wi<-match(selected$sample_id,meta$sample);stopifnot(!anyNA(wi),!anyNA(cts));cts<-cts[,wi,drop=FALSE];rm(raw)
a<-fread(file.path(base,target,'base_assignments.csv'))
parse_ids<-function(x)as.integer(strsplit(x,';',fixed=TRUE)[[1]])
sourceCpp(file.path(base,'fit.cpp'),cacheDir=file.path(tempdir(),'cpm_fit'))
prepared<-setNames(lapply(contrasts,function(cn)readRDS(file.path(td,cn,'prepared.rds'))),contrasts)
# Hashes cover prepared scores, memberships and all computational sources/inputs.
manifest<-fromJSON(file.path(out,'input_manifest.json'))
fingerprint<-paste(manifest$run_fingerprint,workers,sep=':')
dir.create(file.path(td,'chunks'),showWarnings=FALSE)
sixkeys<-unique(a$six_indices);batchkeys<-split(sixkeys,ceiling(seq_along(sixkeys)/8))
batches<-lapply(batchkeys,function(k)which(a$six_indices%in%k))
logmsg<-function(...) {cat(format(Sys.time(),'%H:%M:%S'),target,worker,...,'\n');flush.console()}
qc<-list();checked<-FALSE
for(b in seq_along(batches)) {
  if((b-1)%%workers!=worker-1)next
  dest<-file.path(td,'chunks',sprintf('batch_%04d.rds',b))
  if(file.exists(dest)) {z<-readRDS(dest);stopifnot(identical(z$fingerprint,fingerprint));logmsg('cached',b);next}
  begin<-proc.time()[3];indices<-batches[[b]];f<-matrix(0,length(symbols),length(indices));prev<-'';native_error<-0
  for(j in seq_along(indices)) {
    i<-indices[j];ii<-parse_ids(a$six_indices[i]);ko<-parse_ids(a$KO_indices[i]);cond<-as.integer(ii%in%ko)
    if(a$six_indices[i]!=prev) {prev<-a$six_indices[i];y<-suppressMessages(calcNormFactors(DGEList(cts[,ii,drop=FALSE])))}
    src<-factor(selected$group_id[ii]);design<-model.matrix(~src+cond)
    stopifnot(qr(design)$rank==ncol(design),6-ncol(design)==a$residual_df[i])
    v<-voom(y,design,plot=FALSE);fc<-weighted_condition(v$E,v$weights,as.integer(src),cond)
    stopifnot(all(is.finite(fc)))
    if(j==1) {err<-max(abs(fc-lmFit(v,design)$coefficients[,'cond']));stopifnot(err<1e-8);native_error<-max(native_error,err)}
    f[,j]<-fc[reps]
  }
  # Keep only one orientation of the fits, but rerank both orientations.
  f<-cbind(f,-f);ids2<-c(paste0(a$base_id[indices],'_forward'),paste0(a$base_id[indices],'_reverse'))
  scores<-lapply(prepared,function(z) {
    r<-rankmat(f[z$keep,,drop=FALSE],z$symbols)
    es<-score_rank_matrix(r,z$members);dimnames(es)<-list(z$pathways$pathway,ids2)
    stopifnot(all(is.finite(es)))
    es
  })
  atomic_rds(list(ES=scores,assignment_ids=ids2,residual_df=rep(a$residual_df[indices],2),balanced=rep(a$balanced[indices],2),fingerprint=fingerprint,max_native_FC_error=native_error),dest)
  logmsg('completed batch',b,'/',length(batches),'assignments',length(ids2),'seconds',round(proc.time()[3]-begin,1))
}
writeLines(capture.output(sessionInfo()),file.path(td,paste0('session_worker_',worker,'.txt')))
writeLines(fingerprint,file.path(td,paste0('COMPLETED_worker_',worker)))
logmsg('COMPLETE')
