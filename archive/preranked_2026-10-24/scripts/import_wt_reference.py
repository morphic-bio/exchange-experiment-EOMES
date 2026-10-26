"""Freeze EOMES WT inputs and reference checks from the expanded WT analysis.

Read-only with respect to the source folder. This is an explicit import step;
normal reproduction uses only the resulting local frozen inputs.
"""
import argparse
import json
import shutil
from pathlib import Path
import pandas as pd
from common import ROOT, sha256


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    args = ap.parse_args()
    source = args.source.resolve()
    dest = ROOT / 'inputs/wt_reference'
    prov = ROOT / 'provenance/wt_reference'
    dest.mkdir(exist_ok=True)
    prov.mkdir(exist_ok=True)
    spec = json.loads((source / 'inputs/config.json').read_text())['targets']['EOMES']
    records = []

    def record(src, target, storage):
        records.append(dict(source=str(src), source_sha256=sha256(src),
                            release_path=str(target.relative_to(ROOT)), sha256=sha256(target),
                            bytes=target.stat().st_size, storage=storage))

    for src, name in [('inputs/EOMES/gene_universe.tsv', 'gene_universe.tsv'),
                      ('EOMES/base_assignments.csv', 'base_assignments.csv'),
                      ('EOMES/all_assignments.csv', 'all_assignments.csv'),
                      ('inputs/EOMES_48h_top500_targets.tsv', 'source_top500.tsv')]:
        shutil.copy2(source / src, dest / name)
        record(source / src, dest / name, 'unchanged source bytes')
    selected = pd.read_csv(source / 'inputs/selected_WT_samples.csv')
    selected = selected[selected.target_gene == 'EOMES'].copy()
    assert len(selected) == 9 and selected.sample_id.is_unique
    selected.to_csv(dest / 'selected_WT_samples.csv', index=False)
    record(source / 'inputs/selected_WT_samples.csv', dest / 'selected_WT_samples.csv', 'EOMES rows, original order')
    counts = pd.read_csv(spec['counts'], index_col=0)
    counts.columns = counts.columns.str.replace(r'_S[0-9]+_L[0-9]+$', '', regex=True)
    universe = pd.read_csv(dest / 'gene_universe.tsv', sep='\t')
    counts = counts.loc[universe.gene_id, selected.sample_id]
    assert counts.notna().all().all() and counts.index.is_unique
    counts.index.name = 'gene_id'
    counts.to_csv(dest / 'WT_counts.tsv.gz', sep='\t', compression={'method': 'gzip', 'mtime': 0})
    record(Path(spec['counts']), dest / 'WT_counts.tsv.gz', 'original q75 fitting universe; nine selected EOMES WT columns')
    for cn in ['JAX', 'MSK', 'NW', 'NW_parent', 'UCSF']:
        target = dest / f'source_observed_ranks_{cn}.tsv'
        shutil.copy2(source / 'EOMES' / cn / 'observed_ranks.tsv', target)
        record(source / 'EOMES' / cn / 'observed_ranks.tsv', target, 'unchanged source bytes')
    for name in ['METHODS.md', 'common.R', 'prepare.R', 'run.R', 'finalize.R', 'score.cpp', 'fit.cpp', 'input_manifest.json', 'INDEPENDENT_QC.json']:
        shutil.copy2(source / name, prov / name)
        record(source / name, prov / name, 'unchanged upstream provenance; not executed by the release')
    # A small exact-score benchmark from the existing WT results, all assignments
    # and all five contrasts. Extraction is delegated to R to preserve RDS values.
    import subprocess
    subprocess.run(['Rscript', str(ROOT / 'scripts/export_wt_benchmark.R'), str(source), str(dest)], check=True)
    record(source / 'EOMES/full_observed_results.csv', dest / 'source_benchmark_results.csv', 'custom ChIP benchmark rows')
    records.append(dict(source=str(source / 'EOMES/chunks'), release_path='inputs/wt_reference/source_benchmark_scores.rds',
                        sha256=sha256(dest / 'source_benchmark_scores.rds'), bytes=(dest / 'source_benchmark_scores.rds').stat().st_size,
                        storage='custom ChIP WT ES extracted from all EOMES source chunks'))
    chunks = [{'path': str(p), 'sha256': sha256(p)} for p in sorted((source / 'EOMES/chunks').glob('*.rds'))]
    (prov / 'import.json').write_text(json.dumps(dict(source=str(source), files=records, benchmark_chunks=chunks), indent=2)+'\n')
    manifest_path = ROOT / 'provenance/source_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    paths = {r['release_path'] for r in records}
    manifest = [r for r in manifest if r['release_path'] not in paths] + records
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Frozen {len(counts)} genes × {len(selected)} WT samples and source benchmark.')


if __name__ == '__main__':
    main()
