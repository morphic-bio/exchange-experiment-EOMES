"""48h directional DE enrichment and preparation of fixed cell-fate rankings."""
import json
import platform
import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact, hypergeom
from common import ROOT, CONFIG, CONTRASTS, LABELS, verify_inputs, load_sets

def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    out = np.empty(len(p))
    out[order] = np.minimum(1, np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1])
    return out

def test_overlap(universe, targets, de):
    selected = targets & universe
    a = len(selected & de)
    b = len(selected - de)
    c = len(de - selected)
    d = len(universe - (selected | de))
    odds, p = fisher_exact([[a,b],[c,d]], alternative='greater')
    independent_p = hypergeom.sf(a-1, len(universe), len(selected), len(de))
    assert np.isclose(p, independent_p, rtol=1e-8, atol=1e-300)
    expected = len(selected) * len(de) / len(universe)
    return dict(background_genes=len(universe), eligible_targets=len(selected),
                directional_DE_genes=len(de), overlap=a, expected_overlap=expected,
                fold_enrichment=a/expected if expected else np.nan, odds_ratio=odds,
                pvalue=p, target_DE=a, target_other=b, nontarget_DE=c, nontarget_other=d,
                overlap_genes=';'.join(sorted(selected & de)))

def main():
    checked = verify_inputs()
    tabledir = ROOT / 'results/tables'
    rankdir = tabledir / 'ranked_genes'
    rankdir.mkdir(exist_ok=True)
    ranked = pd.read_csv(ROOT/'inputs/chip48/top1500_symbols.tsv', sep='\t')
    assert len(ranked)==1500 and ranked.gene_symbol.is_unique
    assert ranked.symbol_rank.tolist()==list(range(1,1501))
    raw = pd.read_csv(ROOT/'inputs/chip48/all_gene_scores.tsv.gz', sep='\t')
    raw['gene_symbol'] = raw.gene_symbol.str.upper()
    reproduced = raw.groupby('gene_symbol',as_index=False).weighted_score.max().sort_values(
        ['weighted_score','gene_symbol'],ascending=[False,True]).head(1500)
    assert reproduced.gene_symbol.tolist()==ranked.gene_symbol.tolist()
    np.testing.assert_allclose(reproduced.weighted_score,ranked.weighted_score,rtol=1e-9)
    sets = {n:set(ranked.head(n).gene_symbol) for n in CONFIG['target_sizes']}
    remap = set((ROOT/'inputs/sensitivity/ReMap_EOMES_1446.txt').read_text().split())
    assert len(remap)==CONFIG['remap_reference_size']
    refs = load_sets()
    assert len(refs)==CONFIG['reference_signature_count']
    membership = {}
    for signature, genes in refs.items():
        for gene in genes: membership.setdefault(gene,[]).append(signature)
    main_rows, sensitivity_rows, annotations, target_membership, qc, name_audit = [],[],[],[],[],[]
    for contrast in CONTRASTS:
        x=pd.read_csv(ROOT/f'inputs/de/{contrast}.tsv.gz',sep='\t')
        assert x.ensembl_gene_id.is_unique
        x.insert(0,'contrast',contrast)
        x.insert(1,'method',LABELS[contrast])
        x['symbol']=x.gene_symbol.fillna(x.ensembl_gene_id).astype(str).str.upper()
        x['eligible']=x.q75_pass & (x.max_mean_CPM>=CONFIG['eligibility']['min_max_mean_CPM'])
        x['is_DE']=x.threshold
        assert (x.is_DE==(x.eligible & (x.qvalue<CONFIG['DE']['max_storey_q']) &
                                  (x.logFC.abs()>CONFIG['DE']['min_abs_log2FC']))).all()
        for n,genes in sets.items(): x[f'ChIP48_top{n}']=x.symbol.isin(genes)
        x['reference_cell_signatures']=x.symbol.map(lambda z:';'.join(membership.get(z,[])))
        annotations.append(x)
        representative=x[x.eligible].sort_values(['AveExpr','ensembl_gene_id'],ascending=[False,True]).drop_duplicates('symbol')
        representative.sort_values(['logFC','symbol'],ascending=[False,True]).to_csv(
            rankdir/f'{contrast}.tsv',sep='\t',index=False)
        universe=set(x.loc[x.eligible,'symbol'])
        all_de=set(x.loc[x.is_DE,'symbol'])
        directions={direction:set(x.loc[x.is_DE & (x.logFC.gt(0) if direction=='up' else x.logFC.lt(0)),'symbol']) for direction in ['up','down']}
        for n,genes in sets.items():
            for direction,de in directions.items():
                main_rows.append(dict(contrast=contrast,method=LABELS[contrast],timepoint='48h',top_n=n,direction=direction,
                                      **test_overlap(universe,genes,de)))
        for direction,de in directions.items():
            sensitivity_rows.append(dict(contrast=contrast,method=LABELS[contrast],reference='ReMap',direction=direction,
                                         **test_overlap(universe,remap,de)))
        ids=x.groupby('symbol').ensembl_gene_id.agg(';'.join).to_dict()
        present=set(x.symbol)
        for row in ranked.itertuples(index=False):
            gene=row.gene_symbol
            status='eligible' if gene in universe else 'expression_filter_failed' if gene in present else 'not_in_export'
            target_membership.append(dict(contrast=contrast,method=LABELS[contrast],symbol=gene,rank=row.symbol_rank,
                weighted_score=row.weighted_score,in_top500=row.symbol_rank<=500,status=status,
                DE_up=gene in directions['up'],DE_down=gene in directions['down'],gene_ids=ids.get(gene,'')))
        if contrast==CONTRASTS[0]:
            rna_map=x.set_index('ensembl_gene_id').symbol
            for gene in ranked.gene_symbol:
                source_ids=raw.loc[raw.gene_symbol==gene,'gene_id']
                mapped=sorted(set(source_ids.map(rna_map).dropna()))
                if mapped and gene not in mapped:
                    name_audit.append(dict(ChIP_symbol=gene,RNA_symbols_for_same_gene_IDs=';'.join(mapped),
                                           source_gene_IDs=';'.join(source_ids),RNA_exact_symbol_present=gene in present))
        qc.append(dict(contrast=contrast,method=LABELS[contrast],input_gene_IDs=len(x),eligible_symbols=len(universe),
            DE_symbols=len(all_de),up_DE=len(directions['up']),down_DE=len(directions['down']),
            mixed_direction_symbols=len(directions['up']&directions['down']),
            eligible_top500=len(universe&sets[500]),eligible_top1500=len(universe&sets[1500]),
            DE_top500=len(all_de&sets[500]),DE_top1500=len(all_de&sets[1500]),
            DE_outside_top1500=len(all_de-sets[1500])))
    main=pd.DataFrame(main_rows)
    assert len(main)==CONFIG['target_test']['main_BH_family']
    main['BH_q']=bh(main.pvalue)
    main.to_csv(tabledir/'chip48_directional_DE_enrichment.tsv',sep='\t',index=False)
    sens=pd.DataFrame(sensitivity_rows)
    assert len(sens)==CONFIG['target_test']['sensitivity_BH_family']
    sens['BH_q']=bh(sens.pvalue)
    sens.to_csv(ROOT/'sensitivity/remap/tables/directional_DE_enrichment.tsv',sep='\t',index=False)
    allgenes=pd.concat(annotations,ignore_index=True)
    allgenes.to_csv(tabledir/'gene_evidence.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
    allgenes[allgenes.is_DE & allgenes.ChIP48_top1500].to_csv(tabledir/'chip48_target_DE_genes.tsv',sep='\t',index=False)
    pd.DataFrame(target_membership).to_csv(tabledir/'chip48_target_eligibility.tsv',sep='\t',index=False)
    pd.DataFrame(qc).to_csv(tabledir/'contrast_summary.tsv',sep='\t',index=False)
    pd.DataFrame(name_audit,columns=['ChIP_symbol','RNA_symbols_for_same_gene_IDs','source_gene_IDs','RNA_exact_symbol_present']).to_csv(ROOT/'qc/gene_name_audit.tsv',sep='\t',index=False)
    (ROOT/'qc/analysis_validation.json').write_text(json.dumps(dict(status='PASS',frozen_inputs_verified=checked,
        target_ranking_matches_frozen_source=True,Fisher_matches_hypergeometric=True,
        original_DE_calls_preserved=True,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),indent=2)+'\n')
    print(pd.DataFrame(qc).to_string(index=False))

if __name__=='__main__': main()
