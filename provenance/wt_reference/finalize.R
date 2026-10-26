source('GSEA_expanded_WT_3v3/common.R')
requested<-commandArgs(TRUE)
if(!length(requested))stopifnot(file.exists(file.path(out,'REFERENCES_COMPLETED')))
run_targets<-if(length(requested))requested else targets
stopifnot(all(run_targets%in%targets))
manifest<-fromJSON(file.path(out,'input_manifest.json'))
finalization_fingerprint<-paste(manifest$run_fingerprint,unname(tools::md5sum(file.path(out,'finalize.R'))),sep=':')
allres<-allsens<-allsubset<-allqc<-sample_assignments<-list()
refs<-c('WT','balanced_WT','df_ge2_WT')
focus_paths<-c('GOBP_SYNCYTIUM_FORMATION','CUSTOM_CHIP__EOMES_48h_top500','REMAP__SIX2')
for(t in run_targets) {
  td<-file.path(out,t);a<-fread(file.path(base,t,'all_assignments.csv'))
  marker<-file.path(td,'FINALIZATION_COMPLETED')
  if(file.exists(marker)&&identical(readLines(marker),finalization_fingerprint)) {cat(t,'finalization cache verified\n');next}
  prepared<-setNames(lapply(contrasts,function(cn)readRDS(file.path(td,cn,'prepared.rds'))),contrasts)
  fingerprint<-paste(manifest$run_fingerprint,manifest$workers[[t]],sep=':')
  for(w in seq_len(manifest$workers[[t]]))stopifnot(identical(readLines(file.path(td,paste0('COMPLETED_worker_',w))),fingerprint))
  acc<-lapply(prepared,function(z)setNames(lapply(refs,function(ref)list(np=numeric(nrow(z$pathways)),nn=numeric(nrow(z$pathways)),sp=numeric(nrow(z$pathways)),sn=numeric(nrow(z$pathways)),ex=numeric(nrow(z$pathways)),K=0)),refs))
  focus<-setNames(lapply(contrasts,function(x)list()),contrasts)
  # Independent direct-tail checks on deterministic pathways, across ALL assignments.
  direct<-lapply(prepared,function(z) {
    ix<-sort(unique(c(1L,nrow(z$pathways),match(focus_paths,z$pathways$pathway),round(seq(1,nrow(z$pathways),length.out=8)))))
    ix<-ix[!is.na(ix)]
    list(index=ix,ex=numeric(length(ix)),ns=numeric(length(ix)),ss=numeric(length(ix)))
  })
  seen<-character();max_native<-0;bfiles<-sort(list.files(file.path(td,'chunks'),pattern='rds$',full.names=TRUE))
  for(b in seq_along(bfiles)) {
    v<-readRDS(bfiles[b]);stopifnot(identical(v$fingerprint,fingerprint),!any(v$assignment_ids%in%seen))
    seen<-c(seen,v$assignment_ids);max_native<-max(max_native,v$max_native_FC_error)
    for(cn in contrasts) {
      z<-prepared[[cn]];n<-v$ES[[cn]];obs<-matrix(z$ES[,z$full_profile],ncol=1)
      stopifnot(identical(rownames(n),z$pathways$pathway),identical(colnames(n),v$assignment_ids))
      for(ref in refs) {
        use<-switch(ref,WT=rep(TRUE,ncol(n)),balanced_WT=v$balanced,df_ge2_WT=v$residual_df>=2)
        if(!any(use))next
        nn<-n[,use,drop=FALSE];w<-calibrate_wt(obs,nn);ac<-acc[[cn]][[ref]]
        ac$np<-ac$np+rowSums(nn>=0);ac$nn<-ac$nn+rowSums(nn<=0)
        ac$sp<-ac$sp+rowSums(pmax(nn,0));ac$sn<-ac$sn+rowSums(abs(pmin(nn,0)))
        ac$ex<-ac$ex+as.numeric(w$n_exceedances);ac$K<-ac$K+ncol(nn);acc[[cn]][[ref]]<-ac
      }
      di<-direct[[cn]]
      for(k in seq_along(di$index)) {
        j<-di$index[k];ee<-obs[j,1];vv<-n[j,]
        same<-if(ee>=0)vv>=0 else vv<=0;ex<-if(ee>=0)vv>=ee-1e-12 else vv<=ee+1e-12
        di$ex[k]<-di$ex[k]+sum(ex);di$ns[k]<-di$ns[k]+sum(same);di$ss[k]<-di$ss[k]+sum(abs(vv[same]))
      }
      direct[[cn]]<-di
      for(p in intersect(focus_paths,rownames(n)))focus[[cn]][[length(focus[[cn]])+1]]<-data.table(target=t,contrast=cn,pathway=p,assignment_id=v$assignment_ids,ES=n[p,],balanced=v$balanced,residual_df=v$residual_df)
    }
  }
  stopifnot(setequal(seen,a$assignment_id),length(seen)==nrow(a))
  sa<-fread(file.path(base,'inputs',t,'observed_sample_assignments.tsv'));sa[,target:=t];sample_assignments[[t]]<-sa
  targetres<-list()
  for(cn in contrasts) {
    z<-prepared[[cn]];cd<-file.path(td,cn);mm<-z$profiles;full<-mm[kind=='full']
    sj<-sa[contrast==cn];stopifnot(sum(sj$group=='KO')==full$n_KO,sum(sj$group=='WT')==full$n_control)
    d<-copy(z$pathways);d[,`:=`(target=t,contrast=cn,profile_id=z$full_profile,n_KO=full$n_KO,n_control=full$n_control,ranked_genes=length(z$symbols),ES=as.numeric(z$ES[,z$full_profile]))]
    for(ref in refs) {
      ac<-acc[[cn]][[ref]];ns<-ifelse(d$ES>=0,ac$np,ac$nn);mu<-ifelse(d$ES>=0,ac$sp/ac$np,ac$sn/ac$nn)
      stopifnot(all(ns>0),all(mu>0))
      pv<-ifelse(d$ES==0,1,(1+ac$ex)/(1+ns))
      d[,(paste0(c('p_','NES_','n_same_sign_','n_exceed_','n_total_'),ref)):=list(pv,ifelse(ES==0,0,ES/mu),ns,ac$ex,ac$K)]
      pcol<-paste0('p_',ref);d[,(paste0('q_',ref)):=p.adjust(get(pcol),'BH'),by=collection]
    }
    pr<-fread(file.path(cd,'preranked_results.csv'));stopifnot(identical(pr$pathway,d$pathway),max(abs(pr$ES-d$ES))<1e-12)
    d[,c('p_prerank','NES_prerank','q_prerank','n_same_sign_prerank','n_exceed_prerank','leadingEdge'):=pr[,.(p_prerank,NES_prerank,q_prerank,n_same_sign_prerank,n_exceed_prerank,leadingEdge)]]
    di<-direct[[cn]];j<-di$index
    stopifnot(identical(as.numeric(d$n_exceed_WT[j]),di$ex),identical(as.numeric(d$n_same_sign_WT[j]),di$ns),max(abs(d$NES_WT[j]-ifelse(d$ES[j]==0,0,d$ES[j]/(di$ss/di$ns))))<1e-12)
    fwrite(data.table(pathway=d$pathway[j],n_exceed_WT=di$ex,n_same_sign_WT=di$ns,sum_absolute_same_sign_ES=di$ss),file.path(cd,'independent_tail_checks.csv'))
    fwrite(d,file.path(cd,'full_observed_results.csv'));targetres[[cn]]<-d
    if(length(focus[[cn]]))fwrite(rbindlist(focus[[cn]]),file.path(cd,'focused_WT_ES.csv'))
    subs<-mm[kind=='downsampled'&design=='3v3']
    if(nrow(subs)) {
      ss<-z$ES[,subs$profile_id,drop=FALSE];f<-d$ES;sens<-copy(z$pathways)
      sens[,`:=`(target=t,contrast=cn,full_ES=f,n_subsets=ncol(ss),subset_median_ES=apply(ss,1,median),subset_min_ES=apply(ss,1,min),subset_max_ES=apply(ss,1,max),subset_q025_ES=apply(ss,1,quantile,.025),subset_q975_ES=apply(ss,1,quantile,.975),fraction_same_direction=rowMeans(sign(ss)==sign(f)),mean_abs_ES_difference=rowMeans(abs(ss-f)))]
      allsens[[paste(t,cn)]]<-sens
      for(k in seq_len(nrow(subs)))for(coll in unique(d$collection)) {
        use<-d$collection==coll
        allsubset[[length(allsubset)+1]]<-data.table(target=t,contrast=cn,collection=coll,subset_profile=subs$profile_id[k],n_sets=sum(use),pearson=if(sum(use)>1)cor(f[use],ss[use,k]) else NA_real_,median_abs_ES_difference=median(abs(f[use]-ss[use,k])),spearman=if(sum(use)>1)cor(f[use],ss[use,k],method='spearman') else NA_real_,direction_agreement=mean(sign(f[use])==sign(ss[use,k])))
      }
    }
  }
  result<-rbindlist(targetres);fwrite(result,file.path(td,'full_observed_results.csv'));allres[[t]]<-result
  allqc[[t]]<-list(assignments=nrow(a),null_chunks=length(bfiles),all_assignments_scored_once=TRUE,max_native_FC_error=max_native,independent_direct_tail_checks=TRUE,universe_fixed_across_WT_assignments=TRUE)
  fwrite(rbindlist(allsens)[target==t],file.path(td,'full_vs_3v3_ES_by_pathway.csv'))
  fwrite(rbindlist(allsubset)[target==t],file.path(td,'full_vs_3v3_ES_by_subset.csv'))
  fwrite(sa,file.path(td,'full_observed_sample_assignments.csv'))
  write_json(allqc[[t]],file.path(td,'finalization_QC.json'),pretty=TRUE,auto_unbox=TRUE)
  writeLines(finalization_fingerprint,marker)
  cat(t,'finalized; WT significant',sum(result$q_WT<.05),'\n');flush.console()
}
if(!all(vapply(targets,function(t) {
  marker<-file.path(out,t,'FINALIZATION_COMPLETED')
  file.exists(marker)&&identical(readLines(marker),finalization_fingerprint)
},logical(1)))) {cat('Requested targets finalized; remaining targets still pending.\n');quit(save='no',status=0)}
allres<-setNames(lapply(targets,function(t)fread(file.path(out,t,'full_observed_results.csv'))),targets)
allsens<-lapply(targets,function(t)fread(file.path(out,t,'full_vs_3v3_ES_by_pathway.csv')))
allsubset<-lapply(targets,function(t)fread(file.path(out,t,'full_vs_3v3_ES_by_subset.csv')))
sample_assignments<-lapply(targets,function(t)fread(file.path(out,t,'full_observed_sample_assignments.csv')))
allqc<-setNames(lapply(targets,function(t)fromJSON(file.path(out,t,'finalization_QC.json'))),targets)
d<-rbindlist(allres);fwrite(d,file.path(out,'full_observed_results.csv'))
fwrite(rbindlist(allsens),file.path(out,'full_vs_3v3_ES_by_pathway.csv'))
fwrite(rbindlist(allsubset),file.path(out,'full_vs_3v3_ES_by_subset.csv'))
fwrite(rbindlist(sample_assignments),file.path(out,'full_observed_sample_assignments.csv'))
counts<-d[,.(n_sets=.N,prerank_significant=sum(q_prerank<.05),expanded_WT_significant=sum(q_WT<.05),balanced_WT_significant=sum(q_balanced_WT<.05),df_ge2_WT_significant=sum(q_df_ge2_WT<.05),expanded_WT_positive=sum(q_WT<.05&ES>0),expanded_WT_negative=sum(q_WT<.05&ES<0)),by=.(target,contrast,n_KO,n_control,collection)]
fwrite(counts,file.path(out,'discovery_counts.csv'))
fwrite(d[pathway=='CUSTOM_CHIP__EOMES_48h_top500'],file.path(out,'EOMES_48h_top500_full_results.csv'))
fwrite(d[pathway=='REMAP__SIX2'],file.path(out,'SIX2_full_results.csv'))
fwrite(d[target=='GCM1'&pathway=='GOBP_SYNCYTIUM_FORMATION'],file.path(out,'GCM1_syncytium_full_results.csv'))
write_json(allqc,file.path(out,'QC.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(manifest$run_fingerprint,file.path(out,'ANALYSIS_COMPLETED'))
