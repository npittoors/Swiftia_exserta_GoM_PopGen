# Population Structure & Differentiation Analysis — `02_popgen_pop_structure_p.Rmd`

This script is **Part 2** of a multi-part pipeline. It picks up after sample/locus
filtering and clone correction (Part 1: `01_Data_Setup`) and runs the
downstream population-level analyses: population structure (PCA/DAPC/UMAP), pairwise
FST, LD pruning, AMOVA, and isolation-by-distance (Mantel) tests.

ADMIXTURE, fineRADstructure, and FEEMS analyses are **not** in this file — they live in
separate scripts.

---

## Sample groupings

All population-level analyses in this script are run on some or all of four nested
groupings of the clone-corrected, neutral SNP dataset:

| Dataset  | N individuals | N populations | Definition                                             |
|----------|---------------|----------------|---------------------------------------------------------|
| **All**  | 192           | 13             | Every population, including Edisto                     |
| **GoMa** | 176           | 12             | All Gulf populations regardless of size (no Edisto)     |
| **GoM10**| 158           | 9              | Gulf populations with >10 individuals (excl. WFGB, Bright, Geyer, Edisto) |
| **Sub**  | 174           | 10             | GoM10 + Edisto                                          |

In manuscript text:
**GoM10** = Northern Gulf subset 
**Sub** = Gulf-Carolinian subset

FST, LD pruning, AMOVA, and the Mantel tests are only run on **GoM10** and **Sub**,
which are the two groupings used in the manuscript. Population structure (PCA/DAPC/UMAP
+ geographic correlations) is run on all four.

---

## What this script does

1. **Setup** — loads packages and Part 1's filtered/clone-corrected objects, defines
   site color palettes and west-to-east population orderings, and validates each
   dataset's metadata table against its genind object's individual order before use.
2. **Population structure** (All → GoMa → GoM10 → Sub) — PCA, PCO, and UMAP for
   unsupervised exploration; DAPC without a site prior (data-driven K via
   `find.clusters()`); DAPC with a site prior (used for all downstream figures and
   interpretation); Spearman correlations of retained discriminant functions against
   longitude, latitude, and depth.
3. **FST** — pairwise Weir & Cockerham (1984) FST via SNPRelate, for GoM10 and Sub,
   with west-to-east ordered heatmaps and dendrograms.
4. **Linkage disequilibrium pruning** — SNPRelate-based LD pruning (500 kb sliding
   window, |r| ≤ 0.1, MAF ≥ 0.05, missing rate ≤ 0.1) ahead of population structure
   analyses, plus a separate, more lenient pruning pass on the clones-only dataset for
   relatedness (IBD/IBS) work.
5. **AMOVA** — `poppr::poppr.amova()` (ade4 engine): single-level AMOVA for GoM10 and
   Sub, plus a two-level hierarchical AMOVA for Sub (Region [GoM vs. Atlantic] / Population
   / Individual). Significance assessed via 9,999-permutation `randtest()`.
6. **Mantel tests** — isolation-by-distance: linearized FST [FST/(1-FST)] regressed
   against ln(great-circle distance in km), via `vegan::mantel()` (Pearson, 9,999
   permutations). GoM10 is the primary inference; the Sub result is interpreted with
   caution since Atlantic–Gulf divergence likely inflates the IBD signal.

---

## Requirements

R packages used across the script:

```
adegenet, poppr, ade4, pegas, umap, vegan, geosphere, hierfstat,
SNPRelate, gdsfmt, ggplot2, gridExtra, patchwork, tidyverse (dplyr, tidyr),
viridis, reshape2, ggdendro
```

## Required input files

Produced upstream by `01_Data_Setup` and expected in the working directory:

- `snps.Rdata`, `pc.Rdata`, `snps_mlg.Rdata` — filtered / MLG-filtered / clone-corrected genind objects
- `snps_mlg_he.Rdata`, `snps_mlg_he.hw.Rdata` — Ho- and HWE-filtered neutral datasets
- `nonclone_subsetted_data.Rdata` — All/GoMa/GoM10/Sub genind objects, pop code tables,
  metadata, and population order vectors
- `clean_clone_analysis_final.RData`, `initial_clone_analysis_results.Rdata` — clone detection results
- `distgenEUCL_se.Rdata`, `distgenDISS_se.Rdata` — distance matrices
- `basicstat_all.Rdata` — BayeScan outlier results
- `se_metadata_neutral_all.csv`, `se_metadata_neutral_GoMa.csv`,
  `se_metadata_neutral_GoM10.csv`, `se_metadata_neutral_sub.csv` — per-dataset sample metadata (SITE, LONG, LAT, DEPTH)
- `swiftia_snps_90_cov_maf_bi.recode.vcf` — filtered biallelic SNP VCF (also archived on Figshare)
- `pop_code_GoM10.Rdata`, `pop_code_sub.Rdata` — population code tables used by the AMOVA section
- `snps_mlg_GoM10.Rdata` — GoM10 genind object used by the AMOVA section

## Outputs

- **Figures** (PDF/PNG): PCA/PCO/UMAP and DAPC scatterplots per dataset; geographic
  correlation panels; FST heatmaps and dendrograms; LD heatmaps and post-pruning PCAs;
  Mantel IBD scatterplots (individual and combined panel).
- **Tables** (CSV): geographic correlation results per dataset; FST matrices; AMOVA
  variance components, Phi statistics, and combined summaries (in `AMOVA_Output/`);
  Mantel results summary (`mantel_results_summary.csv`).
- **R objects** (`.Rdata`): saved DAPC objects per dataset for reuse without rerunning.

## Notes

- The `optim.a.score()` PC-selection step for DAPC is interactive and its result
  depends on a PC ceiling entered at a separate prompt that is never recorded, so it
  is **not reproducible** on rerun. These exploratory chunks are commented out and kept
  for reference only; the confirmed `n.pca`/`n.da` values used for each dataset's final
  DAPC are hardcoded directly.
