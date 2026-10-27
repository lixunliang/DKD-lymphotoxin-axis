# Lymphotoxin-β axis in advanced diabetic kidney disease

Analysis code for:

> Li X, Ji Z, Zhao L, Bao T, Wang D, Yuan L. *Lymphotoxin-β and non-canonical NF-κB signalling mark lymphoid infiltration in advanced diabetic kidney disease: an integrative analysis of kidney transcriptomes and single-cell atlases with genetic triangulation.* (manuscript submitted)

All inputs are public (GEO, CELLxGENE, Human Protein Atlas, g:Profiler, UKB-PPP, deCODE, FinnGen R12, CKDGen, eQTLGen). No individual-level patient data are included in this repository.

## Layout

```
scripts/   analysis scripts (run from anywhere; paths are resolved relative to the repository root)
raw/       downloaded data (created by the scripts; not tracked)
results/   tables written by the scripts
figures/   intermediate figures (e.g. figures/umap_hires/)
paper_figures/  publication figures (Figures 1-6, Figure S1)
```

## Run order

| Step | Script | What it does | Main outputs (`results/`) |
|---|---|---|---|
| 1 | `01_bulk_meta_analysis.py` | Downloads GSE30122, GSE104948, GSE104954, GSE96804, GSE142025; moderated t-tests; random-effects meta-analysis (glomeruli, tubulointerstitium); advanced vs early DKD | `META_*.csv`, `DE_*.csv`, `ALL_genes_integrated.csv` |
| 1b | `01b_define_progression_tiers.py` | Tier-1 / Tier-2 progression genes (79 / 166) | `progression_genes_tiered.csv` (Additional file 1) |
| 2 | `02_pathway_enrichment.py` | g:Profiler enrichment of progression genes | `ENRICH_progression_*.csv` |
| 3 | `03_singlecell_localisation.py` | Downloads the CELLxGENE DKD snRNA-seq dataset; cell-type localisation | `SC_*.csv` |
| 4 | `04_snRNA_discovery_cohort.py` | Discovery cohort: Harmony/UMAP, donor-level composition, module scores, LIANA communication, decoupleR/CollecTRI TF activity | `S2_*.csv` |
| 4b | `04b_umap_hires.py` | High-resolution UMAP panels (Figures 2A, 4A); run after steps 4 and 5 | `figures/umap_hires/` |
| 5 | `05_validation_cohort_and_virtual_knockout.py` | GSE209781 validation; GRN virtual knockout / over-expression (`--vko-only` reruns the perturbation step) | `S3_*.csv` |
| 6 | `06_mendelian_randomization_extract.py` | Extracts cis-pQTL MR estimates (UKB-PPP Olink, deCODE SomaScan → FinnGen R12) for candidate proteins | `S4_MR_candidates_kidney_disease.csv` |
| 7 | `09_genetic_triangulation.py` | LTB and lymphotoxin-axis instruments vs diabetic nephropathy, type 1 / type 2 diabetes, eGFR with/without diabetes, rapid eGFR decline, incident CKD, UACR; MHC-protein vs type 1 diabetes comparison; eQTLGen cis-eQTL MR. Set the environment variable `MR_PROJ` to the folder holding the outcome GWAS (`raw/outcome/...`) and the Olink cis-instrument table (`extracted/olink_cis_instruments.tsv`); FinnGen T1D/T2D and eQTLGen are downloaded automatically | `S5_*.tsv` |
| 8 | `10_singlecell_robustness.py` | TLS score without lymphotoxin-axis genes, donor-level values, pre-specified FDR, virtual-knockout stability (seeds, network size, PCs, cells, expression-matched null) | `S6_*.csv` |
| 9 | `07_figures_1_to_4.py`, `08_figures_5_6.py`, `08b_figure_S1_workflow.py` | Publication figures | `paper_figures/` |

## Key methodological note

FinnGen diabetic nephropathy cases are compared with population controls, so MHC variants that confer type 1 diabetes risk produce apparent protein–nephropathy associations. `09_genetic_triangulation.py` tests this explicitly; in our data, MR estimates for diabetic nephropathy tracked those for type 1 diabetes across MHC-region proteins (Spearman ρ = 0.81), and the LTB association was not supported by type 2 diabetes, a gene-proximal instrument or diabetes-stratified kidney-function outcomes.

## Virtual knockout

Gene-regulatory networks are inferred per cell type by principal-component regression (3 PCs; 3 subsamples of 600 cells; 1,000 highly variable genes expressed in ≥5% of cells). A −1 SD (knockout) or +2 SD (over-expression) perturbation is propagated for three steps (CellOracle-style). Significance is assessed against 1,000 knockouts of random genes in the same network; opposition to the bulk progression signature is tested by Spearman correlation with advanced-vs-early DKD log2 fold changes.

## Environment

See `requirements.txt`. Bulk scripts were run with Python 3.9 (Anaconda); single-cell scripts in a separate virtual environment with the pinned versions listed.

## License

MIT
