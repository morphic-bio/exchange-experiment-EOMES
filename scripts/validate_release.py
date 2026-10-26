"""Validate input provenance, ChIP score reconstruction, DE statistics and release scope."""
import json
import platform
import numpy as np
import pandas as pd
import matplotlib
from scipy.stats import fisher_exact
from common import ROOT, CONFIG, CONTRASTS, verify_inputs, sha256

def main():
    n=verify_inputs()
    from validate_wt_gsea import validate
    wt_validation=validate()
    scores=pd.read_csv(ROOT/'inputs/chip48/all_gene_scores.tsv.gz',sep='\t').set_index('gene_id')
    links=pd.read_csv(ROOT/'inputs/chip48/peak_gene_evidence.tsv.gz',sep='\t')
    peaks=pd.read_csv(ROOT/'inputs/chip48/final.narrowPeak.gz',sep='\t',header=None)
    summary=json.loads((ROOT/'inputs/chip48/peak_summary.json').read_text())
    assert len(peaks)==summary['final_peaks']==43700
    assert len(links)==summary['gene_peak_links']==71459
    assert len(scores)==summary['gene_universe']==38899
    assert not links[['gene_id','peak_id']].duplicated().any()
    assert (links.distance_bp==(links.summit-links.tss).abs()).all()
    assert links.distance_bp.between(0,50000).all()
    assert (links.passes_idr & links.passes_pooled_q).all()
    recomputed=links.fold_enrichment*(1-links.distance_bp/50000)
    np.testing.assert_allclose(recomputed,links.weighted_evidence,rtol=1e-10,atol=1e-10)
    summed=recomputed.groupby(links.gene_id).sum().reindex(scores.index,fill_value=0)
    np.testing.assert_allclose(summed,scores.weighted_score,rtol=1e-10,atol=1e-8)
    peak_ids=set(peaks[0]+':'+peaks[1].astype(str)+'-'+peaks[2].astype(str))
    assert set(links.peak_id)<=peak_ids
    targets=pd.read_csv(ROOT/'results/tables/chip48_directional_DE_enrichment.tsv',sep='\t')
    assert len(targets)==20 and set(targets.timepoint)=={'48h'}
    assert set(targets.top_n)=={500,1500} and set(targets.direction)=={'up','down'}
    assert not targets[['contrast','top_n','direction']].duplicated().any()
    for r in targets.itertuples():
        _,p=fisher_exact([[r.target_DE,r.target_other],[r.nontarget_DE,r.nontarget_other]],alternative='greater')
        assert np.isclose(p,r.pvalue,rtol=1e-8,atol=1e-300)
        assert r.target_DE+r.target_other+r.nontarget_DE+r.nontarget_other==r.background_genes
    gsea=pd.read_csv(ROOT/'results/tables/cell_fate_GSEA.tsv',sep='\t')
    assert set(gsea.variant)=={'all_genes','outside_ChIP48_top1500'}
    assert gsea[['NES','pval','BH_q']].notna().all().all()
    assert not gsea[['contrast','pathway','variant']].duplicated().any()
    summary=pd.read_csv(ROOT/'results/tables/contrast_summary.tsv',sep='\t').set_index('contrast')
    for c in CONTRASTS:
        x=gsea[gsea.contrast==c]
        assert set(x.loc[x.variant=='all_genes','ranked_genes'])=={summary.loc[c,'eligible_symbols']}
        assert set(x.loc[x.variant=='outside_ChIP48_top1500','ranked_genes'])=={
            summary.loc[c,'eligible_symbols']-summary.loc[c,'eligible_top1500']}
    remap=pd.read_csv(ROOT/'sensitivity/remap/tables/cell_fate_GSEA.tsv',sep='\t')
    a=gsea[gsea.variant=='all_genes'].sort_values(['contrast','pathway']).reset_index(drop=True)
    b=remap[remap.variant=='all_genes'].sort_values(['contrast','pathway']).reset_index(drop=True)
    pd.testing.assert_frame_equal(a,b)
    ev=pd.read_csv(ROOT/'results/tables/outside_chip48_DE_program_genes.tsv',sep='\t')
    assert ev.is_DE.all() and not ev.ChIP48_top1500.any()
    assert (ev.signature_BH_q<.05).all()
    check=dict(status='PASS',frozen_inputs=n,final_48h_peaks=43700,peak_gene_links=71459,
        score_formula_reconstructed=True,directional_main_tests=20,
        cell_fate_tests_by_variant=gsea.groupby('variant').size().to_dict(),
        excluded_gene_counts_match_top1500=True,primary_contains_no_ReMap_variant=True,
        sensitivity_all_gene_panel_matches_primary=True,python=platform.python_version(),matplotlib=matplotlib.__version__)
    check['WT_reference_validation']=wt_validation['status']
    (ROOT/'qc/release_validation.json').write_text(json.dumps(check,indent=2)+'\n')
    files=[p for directory in ['inputs','scripts','provenance','results','sensitivity']
           for p in (ROOT/directory).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    files += [ROOT/x for x in ['config.json','run.py','README.md','requirements.txt','R_packages.txt'] if (ROOT/x).exists()]
    manifest={str(p.relative_to(ROOT)):dict(sha256=sha256(p),bytes=p.stat().st_size) for p in sorted(files)}
    (ROOT/'release_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(check,indent=2))

if __name__=='__main__':main()
