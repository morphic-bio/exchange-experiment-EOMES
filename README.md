# EOMES analysis release — 24 October 2026

This package reproduces the figures and tables comparing MorPhiC differential expression (DE) with **48-hour EOMES ChIP-seq targets** and reference cell-type expression programs. It uses frozen processed ChIP inputs, original MorPhiC DE exports and a fixed signature library. A separate ReMap analysis assesses sensitivity to the target definition.

## Main outputs

- **[Main figure: one Degron comparison](results/figures/figure_1_chip48_and_cell_fate_single_degron.pdf)**, also available as PNG and SVG. Columns show Whole-gene KO, Knock-in KO, Degron and CRISPRi. Degron compares IAA-treated tagged clones with their matched DMSO controls.
- [Directional target enrichment](results/tables/chip48_directional_DE_enrichment.tsv): overlap counts, expected counts, fold enrichment, Fisher p-values and BH-adjusted q-values for all 20 tests.
- [Cell-type program GSEA](results/tables/cell_fate_GSEA.tsv): every tested signature, including signatures omitted from the figure.
- [Target eligibility](results/tables/chip48_target_eligibility.tsv) and [target DE genes](results/tables/chip48_target_DE_genes.tsv): identify which targets are represented, expression-eligible and differentially expressed in each contrast.
- [DE genes supporting significant programs outside the targets](results/tables/outside_chip48_DE_program_genes.tsv): leading-edge DE genes after excluding the top 1,500 ChIP targets.
- [ReMap sensitivity figure](sensitivity/remap/figures/figure_s1_remap_sensitivity.pdf), with its tables under `sensitivity/remap/tables/`.

Panel A shows enrichment separately for upregulated and downregulated DE genes, using the top 500 and top 1,500 ChIP targets. Panels B and C show cell-type GSEA before and after excluding the top 1,500 ChIP targets. The figure displays 16 reference signatures, including “Anterior-enriched definitive endoderm (DE2).”

The original five-column figure, including the parental-control degron comparison, is retained as `results/figures/figure_1_chip48_and_cell_fate.{pdf,png,svg}`. **The four-column main figure hides that comparison but preserves the original five-contrast statistics and multiple-testing corrections.** Tables and the ReMap sensitivity figure retain all five contrasts.

## Installation and execution

Run these commands from this release directory. Python 3.9 or later and R are required. Python dependencies are specified in `requirements.txt`; R dependencies are listed in `R_packages.txt`.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("fgsea", "BiocParallel", "data.table", "jsonlite"), ask=FALSE, update=FALSE)'
python run.py
```

Installation requires network access. Once dependencies are installed, the analysis runs entirely from the packaged inputs, without downloads or access to the parent exploration folders. It does not realign reads, rerun peak calling, select an input control or refit MorPhiC DE models.

The tested environment used Python 3.9.6, NumPy 1.26.3, pandas 2.2.0, SciPy 1.12.0, matplotlib 3.8.2, R 4.6.1 and fgsea 1.38.0. Consult `qc/analysis_validation.json`, `qc/release_validation.json` and `qc/R_sessionInfo.txt` for complete recorded versions. Dependency ranges permit compatible installations; use the recorded versions when matching numerical results exactly, especially GSEA results.

The launcher runs the analysis, GSEA, plotting and validation stages in sequence, overwriting generated scientific outputs and writing stage logs to `qc/`. It can also be invoked from another directory using an absolute path. To select an R installation or validate existing outputs:

```sh
python run.py --rscript /path/to/Rscript
python run.py --validate-only
```

For individual stages, run these commands in order:

```sh
python scripts/analyze.py
Rscript scripts/cell_fate.R
python scripts/render_results.py
python scripts/render_results.py --single-degron
python scripts/validate_release.py
```

| Code | Purpose |
| --- | --- |
| `run.py` | Execute the complete figure-and-table workflow and stop on failure. |
| `scripts/common.py` | Resolve release-local paths, load configuration and signatures, and verify input checksums. |
| `scripts/analyze.py` | Verify target rankings and original DE calls; calculate directional enrichment; export gene annotations, eligibility and GSEA rankings. |
| `scripts/cell_fate.R` | Run GSEA for all genes, outside ChIP targets and outside ReMap targets; adjust p-values and record the R environment. |
| `scripts/render_results.py` | Plot the five-column main and sensitivity figures; export leading-edge evidence tables. With `--single-degron`, plot the four-column main figure from the same statistics. |
| `scripts/validate_release.py` | Check input provenance and numerical consistency, then refresh the release checksum manifest. |

To change only figure formatting, edit `scripts/render_results.py` and rerun the two rendering commands. Changing target definitions, thresholds or signatures changes the analysis and requires a full rerun and corresponding validation changes. Keep the five contrasts in `config.json` to reproduce the released multiple-testing families.

## Packaged inputs and provenance

| Location | Contents |
| --- | --- |
| `inputs/chip48/` | Final peaks, peak–gene links, all-locus scores, peak summary and fixed top-1,500 symbol ranking; the first 500 symbols define the smaller set. |
| `inputs/de/` | Five original DE exports: `whole_gene_ko`, `knock_in_ko`, `degron_matched`, `degron_parental` and `crispri`. |
| `inputs/reference/` | Frozen `cell_type_signatures.gmt` and `signature_sources.tsv`. |
| `inputs/sensitivity/` | The retained 1,446-symbol ReMap EOMES list. |
| `provenance/` | Source locations, input checksums and selected upstream sample/scoring information. |
| `config.json` | Contrast definitions, thresholds, GSEA parameters and displayed signatures. |

The 48-hour ChIP data are from HUES8 cells differentiating toward definitive endoderm, [Luo et al. (2023)](https://doi.org/10.1038/s41588-023-01450-7), GEO **GSE213394**: **GSM6585561 / SRR21586161** and **GSM6585562 / SRR21586115**. The selected external input is **GSM3381961 / SRR7814534**, from [Li et al. (2019)](https://doi.org/10.1038/s41588-019-0408-9); it is not an experimentally matched 48-hour control.

Upstream processing used GRCh38, GENCODE v44 basic, blacklist filtering, MACS2 with a fixed 200-bp extension, global irreproducible discovery rate (IDR) ≤ 0.05 and pooled peak q ≤ 0.05. The retained inputs contain **43,700 peaks, 71,459 peak–gene links and 38,899 scored gene loci**. Details are recorded in `provenance/selected_input_and_scoring.json` and `provenance/selected_upstream_samples.tsv`. Raw FASTQ/BAM files are not required for this release.

`provenance/source_manifest.json` defines the checksums of the 14 frozen input files. Input mismatches stop execution. Original absolute paths in provenance are historical records, not runtime dependencies.

## Target definition and DE enrichment

Each retained peak summit is linked to every gene-level transcription start site (TSS) within 50 kb. The locus score is:

```text
sum over linked peaks [peak fold enrichment × (1 − |summit − TSS| / 50,000)]
```

The score sums evidence across peaks. For repeated gene symbols, the maximum locus score is retained; distinct symbols are ranked by decreasing score, with alphabetical ordering for tied scores. The fixed top 500 and top 1,500 unique uppercase symbols are candidate targets based on reproducible ChIP binding enrichment and TSS proximity. Targets absent from the RNA export or failing expression eligibility are not replaced by lower-ranked genes. Validation independently reconstructs the scores and verifies the ranking.

The original DE exports contain `ensembl_gene_id`, `gene_symbol`, `logFC`, `AveExpr`, `qvalue`, `q75_pass`, `max_mean_CPM` and `threshold`, among other annotations. Expression eligibility requires a positive 75th percentile of counts across the original count matrix (`q75_pass`) and maximum group-mean TMM-normalized CPM ≥ 1. DE additionally requires **Storey q < 0.05 and |log2FC| > 1**. The exported `threshold` calls are checked against these rules and preserved. Positive logFC means higher expression after perturbation.

For each contrast, the Fisher-test background is **all unique expression-eligible RNA symbols**. Matching uses exact uppercase symbols without alias inference; missing symbols are represented by Ensembl IDs. A symbol is directionally DE if any eligible underlying gene ID meets the DE threshold in that direction. Symbols with qualifying loci in opposite directions can enter both tests; their counts appear in `contrast_summary.tsv`.

Each one-sided Fisher exact test compares target/non-target membership with directional DE/other status. The expected overlap is `(eligible targets × directional DE symbols) / background symbols`, and the plotted fold enrichment is observed/expected overlap. Benjamini–Hochberg (BH) adjustment covers **20 tests: five contrasts × two target sizes × two directions**. Values below one indicate smaller-than-expected overlap; the test is for enrichment, not depletion.

## Cell-type program enrichment

The fixed library contains **96 signatures**: 21 human gastrula signatures using the first 50 published markers from supplementary tables 1, 8, 11 and 13 of [Tyser et al. (2021)](https://doi.org/10.1038/s41586-021-04158-y), and 75 main human fetal cell-type signatures from [Cao et al. (2020)](https://doi.org/10.1126/science.aba7721), distributed in MSigDB C8 v2025.1.Hs.

All expression-eligible symbols are ranked by signed logFC, not only DE genes. Where several eligible loci share a symbol, GSEA uses the locus with highest `AveExpr`, breaking ties by Ensembl ID. This representative-locus rule differs from the any-qualifying-locus rule used for directional DE overlap.

GSEA uses `fgseaMultilevel` with `gseaParam=0` (unweighted), `minSize=10`, `maxSize=1000`, `nPermSimple=1000`, `eps=1e-7`, seed `20261024` and `BiocParallel::SerialParam()`. Panel C reruns GSEA after removing the top-1,500 ChIP symbols from the ranking, which also removes them from each signature's intersection with the ranking. Positive NES indicates enrichment toward increased expression; negative NES indicates enrichment toward decreased expression.

BH correction includes **all tested signatures across all five contrasts**, separately for each analysis variant, rather than only the 16 displayed signatures. There are 441 tests for all genes and 436 after excluding ChIP targets. Use the table's **`BH_q`** column for the displayed significance; `padj` is fgsea's within-contrast adjustment. Expected warnings about the `1e-7` p-value floor are retained in `qc/cell_fate.log`; values at this floor are not precise estimates of smaller p-values.

The separate ReMap analysis uses the fixed 1,446-symbol list with the same eligibility and directional DE rules, BH correction across ten Fisher tests, and GSEA after removing ReMap targets. Its all-gene GSEA results are identical to the main analysis; its exclusion results have their own BH family. ReMap does not define the main target sets or panel C.

## Reading tables and checking results

- `contrast_summary.tsv` reports eligible and DE symbol counts, target overlap and mixed-direction symbols per contrast.
- `chip48_target_eligibility.tsv` contains one row per fixed top-1,500 target per contrast, with rank, score, top-500 membership, eligibility status and directional DE flags.
- `gene_evidence.tsv.gz` retains gene-ID-level RNA evidence with target/signature annotations. `chip48_target_DE_genes.tsv` is its DE subset within the top 1,500 targets; rows are gene IDs, not necessarily unique symbols.
- `ranked_genes/` contains the eligible representative-locus rankings passed to GSEA.
- `cell_fate_GSEA.tsv` includes the signature, contrast, analysis variant, NES, p-values, `BH_q` and semicolon-separated leading-edge symbols.
- `outside_chip48_leading_edge_gene_evidence.tsv.gz` links outside-target leading edges to RNA gene-ID evidence. `outside_chip48_DE_program_genes.tsv` retains DE genes supporting signatures with `BH_q < 0.05`. A gene may appear for multiple signatures or loci; count unique symbols when reporting numbers of genes.

Validation checks frozen inputs, score reconstruction, peak links, fixed target ranking, preservation of DE calls, Fisher/hypergeometric agreement, target-exclusion counts and identical all-gene panels between the main and ReMap analyses. Reports and logs are under `qc/`. `release_manifest.json` records current file checksums; validation refreshes this manifest rather than comparing generated outputs against a historical manifest.

Both ChIP target sets are enriched among downregulated DE genes in all five contrasts (2.60–3.72-fold); no upregulated enrichment test is significant. Mesoderm and endothelial programs remain enriched after excluding the top 1,500 targets. Reference similarity describes expression programs, not proof of a replacement cell identity, and absence from a target list does not establish indirect regulation. Whole-gene and knock-in KO share WT controls; the two degron contrasts share perturbed samples. The 48-hour target definition was selected during earlier exploration and is fixed for this release.
