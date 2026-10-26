"""Independently reconstruct observed ES, empirical WT tails, NES and BH in Python."""
import json
import numpy as np
import pandas as pd
from common import ROOT, CONFIG, load_sets
from analyze import bh


def validate():
    primary = pd.read_csv(ROOT/'results/tables/cell_fate_GSEA.tsv', sep='\t')
    remap = pd.read_csv(ROOT/'sensitivity/remap/tables/cell_fate_GSEA.tsv', sep='\t')
    results = pd.concat([primary, remap[remap.variant=='outside_ReMap']], ignore_index=True)
    sets = load_sets()
    chip = set(pd.read_csv(ROOT/'inputs/chip48/top1500_symbols.tsv', sep='\t').gene_symbol)
    remap_targets = set((ROOT/'inputs/sensitivity/ReMap_EOMES_1446.txt').read_text().split())
    assignments = pd.read_csv(ROOT/'inputs/wt_reference/all_assignments.csv').set_index('assignment_id')
    audit = []
    for (contrast, variant), d in results.groupby(['contrast', 'variant']):
        ranked = pd.read_csv(ROOT/f'results/tables/ranked_genes/{contrast}.tsv', sep='\t')
        assert ranked.eligible.all() and ranked.q75_pass.all() and ranked.max_mean_CPM.ge(1).all()
        excluded = chip if variant=='outside_ChIP48_top1500' else remap_targets if variant=='outside_ReMap' else set()
        ranked = ranked[~ranked.symbol.isin(excluded)].sort_values(['logFC','symbol'], ascending=[False,True])
        n = len(ranked)
        assert set(d.ranked_genes)=={n} and ranked.symbol.is_unique
        null = pd.read_csv(ROOT/f'results/tables/wt_reference_scores/{contrast}__{variant}.tsv.gz', sep='\t', index_col=0)
        assert set(null.columns)==set(assignments.index) and null.shape[1]==1566
        assert set(null.index)==set(d.pathway)
        max_es_error = 0
        for r in d.itertuples():
            hits = ranked.symbol.isin(sets[r.pathway]).to_numpy()
            m = hits.sum()
            assert m==r.size and CONFIG['cell_fate']['min_size']<=m<=CONFIG['cell_fate']['max_size']
            walk = np.cumsum(np.where(hits,1/m,-1/(n-m)))
            hi, lo = max(0,walk.max()), min(0,walk.min())
            es = hi if hi>-lo else lo
            max_es_error = max(max_es_error, abs(es-r.ES))
            np.testing.assert_allclose(es,r.ES,atol=1e-10,rtol=0)
            for ref in ['WT','balanced_WT','df_ge2_WT']:
                a = assignments.loc[null.columns]
                mask = np.ones(len(a),dtype=bool) if ref=='WT' else a.balanced.to_numpy() if ref=='balanced_WT' else a.residual_df.ge(2).to_numpy()
                values = null.loc[r.pathway].to_numpy()[mask]
                same = values>=0 if r.ES>=0 else values<=0
                exceed = values>=r.ES-1e-12 if r.ES>=0 else values<=r.ES+1e-12
                p = 1 if r.ES==0 else (1+exceed.sum())/(1+same.sum())
                nes = 0 if r.ES==0 else r.ES/np.abs(values[same]).mean()
                assert same.sum()==getattr(r,'n_same_sign_'+ref)
                assert exceed.sum()==getattr(r,'n_exceed_'+ref)
                assert len(values)==getattr(r,'n_total_'+ref)
                np.testing.assert_allclose(p,getattr(r,'pval_'+ref),rtol=1e-12)
                np.testing.assert_allclose(nes,getattr(r,'NES_'+ref),rtol=1e-12)
            assert r.reference=='WT_3v3' and r.NES==r.NES_WT and r.pval==r.pval_WT
            assert set(str(r.leadingEdge).split(';'))<=set(ranked.symbol[hits])
        audit.append(dict(contrast=contrast,variant=variant,ranked_genes=n,signatures=len(d),max_observed_ES_error=max_es_error))
    for variant, d in results.groupby('variant'):
        for p,q in [('pval','BH_q'),('pval_prerank','BH_q_prerank'),
                    ('pval_balanced_WT','BH_q_balanced_WT'),('pval_df_ge2_WT','BH_q_df_ge2_WT')]:
            np.testing.assert_allclose(bh(d[p]),d[q],rtol=1e-12,atol=1e-14)
    old = pd.read_csv(ROOT/'archive/preranked_2026-10-24/results/tables/cell_fate_GSEA.tsv',sep='\t')
    keys = ['contrast','variant','pathway']
    joined = primary.merge(old,on=keys,suffixes=('', '_original'),validate='one_to_one')
    assert len(joined)==len(primary)==len(old)
    for col in ['ES','size','ranked_genes']:
        np.testing.assert_allclose(joined[col],joined[col+'_original'],rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(joined.NES_prerank,joined.NES_original,rtol=1e-12)
    np.testing.assert_allclose(joined.BH_q_prerank,joined.BH_q_original,rtol=1e-12)
    rows=[]
    for (c,v), d in results.groupby(['contrast','variant']):
        rows.append(dict(contrast=c,variant=v,tested=len(d),significant_WT=int(d.BH_q.lt(.05).sum()),
                         significant_prerank=int(d.BH_q_prerank.lt(.05).sum()),
                         significant_balanced_WT=int(d.BH_q_balanced_WT.lt(.05).sum()),
                         significant_df_ge2_WT=int(d.BH_q_df_ge2_WT.lt(.05).sum())))
    pd.DataFrame(rows).to_csv(ROOT/'results/tables/cell_fate_reference_comparison.tsv',sep='\t',index=False)
    report=dict(status='PASS',tests=len(results),null_scores_checked=int(sum(len(d)*1566 for _,d in results.groupby(['contrast','variant']))),
                observed_rankings_and_ES_unchanged=True,original_preranked_results_reproduced=True,
                empirical_tails_NES_and_BH_independently_recomputed=True,comparisons=audit)
    (ROOT/'qc/wt_gsea_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    print(json.dumps(validate(),indent=2))
