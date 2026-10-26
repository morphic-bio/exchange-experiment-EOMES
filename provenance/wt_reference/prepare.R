# Build fixed comparison-specific universes AFTER the original DE fitting.
# Do not renormalize or refit observed DE models on the reduced gene sets.
source('GSEA_expanded_WT_3v3/common.R')
dir.create(file.path(out,'inputs'),recursive=TRUE,showWarnings=FALSE)
gmt <- c(gmtPathways(file.path(base,'inputs/gene_sets.gmt')),
         gmtPathways(file.path(base,'inputs/go_bp.v2023.2.Hs.symbols.gmt')))
stopifnot(!anyDuplicated(names(gmt)))
gmt<-gmt[sort(names(gmt))]
counts<-list();validation<-list()
expected<-list(EOMES=c(15165,15009,14523,14827,15047),GCM1=c(13568,13502,13556,13140,13204),ISL1=c(14692,14627,14284,14793,14598),`NKX2-1`=c(16005,16308,15746,15859,16019))
for(t in targets) {
  inp<-file.path(base,'inputs',t);td<-file.path(out,t)
  dir.create(td,recursive=TRUE,showWarnings=FALSE)
  u<-fread(file.path(inp,'gene_universe.tsv'))
  reps<-which(u$representative);reps<-reps[order(u$symbol[reps])]
  symbols<-u$symbol[reps]
  fc<-readRDS(file.path(inp,'observed_logFC.rds'))
  stopifnot(nrow(fc)==length(symbols))
  # Original frozen matrices store rows in alphabetical representative-symbol
  # order, without row names. Verify that order against every full DE vector.
  if(!is.null(rownames(fc))) stopifnot(identical(rownames(fc),symbols))
  rownames(fc)<-symbols
  profiles<-fread(file.path(inp,'observed_profiles.tsv'))[kind=='full'|design=='3v3']
  all_coverage<-list()
  for(cn in contrasts) {
    cd<-file.path(td,cn);dir.create(cd,showWarnings=FALSE)
    de<-fread(paste0(t,'_DE_',cn,'.tsv'))
    stopifnot(!anyDuplicated(de$ensembl_gene_id),all(de$q75_pass))
    ix<-match(u$gene_id[reps],de$ensembl_gene_id);stopifnot(!anyNA(ix))
    selected<-copy(u[reps]);selected[,`:=`(mean_CPM_WT=de$mean_CPM_WT[ix],mean_CPM_KO=de$mean_CPM_KO[ix],q75_pass=de$q75_pass[ix])]
    selected[,max_mean_CPM:=pmax(mean_CPM_WT,mean_CPM_KO)]
    selected[,included:=q75_pass & max_mean_CPM>=1]
    stopifnot(sum(selected$included)==expected[[t]][match(cn,contrasts)])
    fwrite(selected,file.path(cd,'gene_filter_audit.tsv'),sep='\t')
    keep<-which(selected$included);sy<-symbols[keep]
    full<-profiles[contrast==cn & kind=='full'];stopifnot(nrow(full)==1)
    error<-max(abs(fc[,full$profile_id]-de$logFC[ix]))
    stopifnot(error<1e-7)
    mm<-profiles[contrast==cn];f<-fc[keep,mm$profile_id,drop=FALSE]
    members_all<-lapply(gmt,function(s)which(sy%in%s))
    cv<-data.table(pathway=names(gmt),collection=collection_of(names(gmt)),size=lengths(members_all),source_size=lengths(gmt))
    cv[,eligible:=size>=15 & size<=fifelse(collection%in%c('ReMap','Custom_ChIP'),2000L,500L)]
    cv[,canonical_pathway:=pathway]
    # Collapse identical GO memberships deterministically within this universe.
    go<-which(cv$collection=='GO_BP')
    keys<-vapply(members_all[go],paste,collapse=',',FUN.VALUE='')
    cv$canonical_pathway[go]<-cv$pathway[go][match(keys,keys)]
    cv[,tested:=eligible & (collection!='GO_BP' | canonical_pathway==pathway)]
    cv[,`:=`(target=t,contrast=cn)]
    fwrite(cv,file.path(cd,'gene_set_coverage.tsv'),sep='\t');all_coverage[[cn]]<-cv
    paths<-cv[tested==TRUE]$pathway;members<-members_all[paths]
    ranks<-rankmat(f,sy);es<-score_rank_matrix(ranks,members);dimnames(es)<-list(paths,mm$profile_id)
    pm<-cv[tested==TRUE,.(pathway,collection,size)]
    es_full<-es[,full$profile_id]
    oo<-order(-f[,full$profile_id],sy);ff<-f[oo,full$profile_id];names(ff)<-sy[oo]
    leading<-character(length(paths));max_error<-0
    for(k in seq_along(paths)) {
      z<-calcGseaStat(ff,ranks[members[[k]],match(full$profile_id,mm$profile_id)],gseaParam=0,returnLeadingEdge=TRUE)
      max_error<-max(max_error,abs(es_full[k]-z$res))
      leading[k]<-paste(names(ff)[z$leadingEdge],collapse=';')
    }
    stopifnot(max_error<1e-10)
    data<-list(target=t,contrast=cn,keep=keep,symbols=sy,members=members,pathways=pm,profiles=mm,ES=es,leadingEdge=leading,full_profile=full$profile_id,
               universe_hash=unname(tools::md5sum(file.path(cd,'gene_filter_audit.tsv'))),max_fgsea_ES_error=max_error,max_saved_logFC_error=error)
    atomic_rds(data,file.path(cd,'prepared.rds'))
    fwrite(data.table(gene_id=selected$gene_id[keep],symbol=sy,logFC=f[,full$profile_id],mean_CPM_WT=selected$mean_CPM_WT[keep],mean_CPM_KO=selected$mean_CPM_KO[keep]),file.path(cd,'observed_ranks.tsv'),sep='\t')
    counts[[paste(t,cn)]]<-data.table(target=t,contrast=cn,q75_representative_genes=length(symbols),ranked_genes=length(sy),removed=length(symbols)-length(sy),tested_sets=nrow(pm))
    validation[[paste(t,cn)]]<-list(max_observed_logFC_error=error,max_fgsea_ES_error=max_error)
    cat(t,cn,length(sy),'genes;',nrow(pm),'sets\n');flush.console()
  }
  fwrite(rbindlist(all_coverage),file.path(td,'gene_set_coverage.tsv'),sep='\t')
}
fwrite(rbindlist(counts),file.path(out,'gene_universe_counts.csv'))
write_json(validation,file.path(out,'preparation_QC.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines('Complete',file.path(out,'PREPARATION_COMPLETED'))
