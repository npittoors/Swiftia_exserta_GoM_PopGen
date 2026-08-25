# ADMIXTURE Analysis — 03_Admixture_p.Rmd

This script is Part 3 of a multi-part pipeline. It picks up after sample/locus
filtering and clone correction (Part 1: `01_Data_Setup`) and population
structure analyses (Part 2: `02_popgen_pop_structure_p.Rmd`), and runs
ADMIXTURE ancestry estimation on three of the four nested sample groupings.

## Sample groupings

ADMIXTURE is run on three LD-pruned datasets derived from the clone-corrected,
BayeScan-filtered neutral SNP panel (16,746 loci; see `01_Data_Setup`):

| Dataset | N individuals | N populations | Pruned loci | Definition |
|---|---|---|---|---|
| GoM10 | 158 | 9 | 5,419 | Gulf populations with >10 individuals (excl. WFGB, Bright, Geyer, Edisto) |
| GoMa | 176 | 12 | 5,634 | All Gulf populations regardless of size (no Edisto) |
| Sub | 174 | 10 | 5,629 | GoM10 + Edisto |

In manuscript text: GoM10 = Northern Gulf subset; Sub = Gulf-Carolinian subset.
The **All** (192 ind, 13 pop) grouping used elsewhere in the pipeline is not
run through ADMIXTURE.

LD pruning uses PLINK's SNP-count window (`--indep-pairwise 50 5 0.1`) rather
than a physical-distance window, since RAD loci lack real linkage-map
coordinates — this differs from the 500 kb sliding-window approach used for LD
pruning ahead of population structure analyses in Part 2, which is appropriate
there but violates ADMIXTURE's marker-independence assumption if applied here
unpruned.

## What this script does

1. **Export unpruned PLINK trios** — converts each clone-corrected, BayeScan-filtered
   genind object (GoM10, GoMa, Sub) to genlight, then to PLINK `.ped`/`.map`
   format via `dartR::gl2plink()`, with `stopifnot()` checks confirming
   individual and locus counts before export.
2. **LD pruning (server, PLINK)** — `--indep-pairwise 50 5 0.1` on each
   dataset's converted `.bed` file, producing a `.prune.in` locus list per
   dataset.
3. **Extraction** — subsets each `.bed`/`.bim`/`.fam` trio to the pruned locus
   list, producing the final ADMIXTURE-ready inputs.
4. **ADMIXTURE runs (server, sequential)** — K = 1–10, 10 replicates each
   (300 runs total across the three datasets), run one dataset at a time via
   `nohup`/`disown` to avoid oversubscribing server threads.
5. **CV error extraction and summary** — parses CV error from all 300 log
   files, computes per-K mean/SD, identifies the best-supported K per dataset,
   and identifies the lowest-CV replicate at each K.
6. **Best-run file copying** — copies the winning replicate's `.Q`/`.P` files
   into a `best_run_files/` folder per dataset.
7. **Structure plots** — faceted, west-to-east ancestry bar plots via a shared
   `plot_admixture_lineages()` function, with cluster-to-lineage identity
   (NW Gulf / NE Gulf / Carolinian) verified against per-population mean
   ancestry before any name or color is assigned.
8. **Reproducibility archive and supplementary tables** — bundles CV results,
   verified lineage mappings, and population data into a single `.RData`
   archive; exports supplementary CSV tables (CV error, pruning parameters,
   per-dataset ancestry proportions).
9. **Sanity checks** — confirms ancestry proportions sum to 1 per individual
   and that sample counts/no duplicates hold in each output table.

## Requirements

R packages used across the script:

```
dartR, dplyr, tidyr, ggplot2
```

External tools (server-side, not R):

```
PLINK v1.9, ADMIXTURE
```

## Required input files

Produced upstream by `01_Data_Setup` / `02_popgen_pop_structure_p.Rmd` and
expected to be available in the R environment before running this script:

- `snps_mlg_GoM10_n`, `snps_mlg_gom_all_n`, `snps_mlg_sub_n` — clone-corrected,
  BayeScan-filtered genind objects (158 / 176 / 174 individuals, 16,746 loci
  each)
- `pop_data_GoM10`, `pop_data_GoMa`, `pop_data_sub` — individual-to-population
  data frames (columns: `individual`, `population`)
- `se_metadata_neutral_GoM10`, `se_metadata_neutral_GoMa`,
  `se_metadata_neutral_sub` — metadata data frames including a `LONG` column,
  row-order matched to the genind objects above
- `WE_order_GoM10`, `WE_order_GoMa`, `WE_order_sub` — west-to-east population
  name vectors used to order structure-plot facets

## Outputs

**Figures (PDF):** CV error curve (all three datasets, one panel);
west-to-east faceted ADMIXTURE structure bar plots for GoM10 (K=2), Sub (K=3),
and GoMa (K=2, supplementary).

**Tables (CSV):** CV error summary (`SuppTable_ADMIXTURE_CVerror.csv`);
pruning/run parameters (`SuppTable_ADMIXTURE_parameters.csv`); best replicate
per K per dataset (`ADMIXTURE_best_reps_pruned_all3.csv`); per-individual
ancestry proportions for GoM10, Sub, and GoMa
(`SuppTable_ADMIXTURE_ancestry_*.csv`).

**R objects (.RData):** `ADMIXTURE_reproducibility_archive.RData` — bundles CV
results, verified cluster-to-lineage mappings, pruning parameters, and
population data/order vectors for all three datasets. Intended for public
data deposit alongside the pruned `.bim` files and best-run `.Q`/`.P` files
(not embedded in the archive, to keep file size manageable).

## Notes

- **K=1 is the best-supported result in all three datasets** — CV error
  increases monotonically through K=10 with no exceptions. This is treated as
  a genuine robustness finding (not a null result to work around), since it
  shows the same pattern holds regardless of dataset scope. K=2/K=3 panels are
  retained for interpretability, not because they represent a better-supported
  clustering solution.
- **Column order is not guaranteed to be consistent across independent
  ADMIXTURE runs.** Cluster-to-lineage identity must be re-verified via
  per-population mean ancestry any time ADMIXTURE is rerun — never assume a
  prior run's `column_order` mapping still applies.
- **LD pruning here is ADMIXTURE-specific** and uses a different method
  (PLINK SNP-count window) than the LD pruning used ahead of population
  structure analyses in Part 2 (SNPRelate, 500 kb physical window). The two
  should not be conflated or treated as interchangeable pruned sets.
- Two earlier, now-superseded scripts (`01_data_prep.R`, `02_process_results.R`)
  ran ADMIXTURE on unpruned data for GoM10/GoMa, and used a mislabeled Q file
  for Sub (the all-populations 192-ind run cosmetically filtered to look like
  a 174-ind Sub run). Those scripts are not part of the public pipeline; this
  script is the corrected, sole version used for the manuscript.
