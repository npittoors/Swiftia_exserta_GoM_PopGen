# Data Import, Filtering, Clone Correction & BayeScan Outlier Removal — `01_Data_Setup.Rmd`

This script is **Part 1** of a multi-part pipeline. It starts from the raw SNP
VCF and sample metadata and produces the filtered, clone-corrected, neutral
SNP dataset — and its four population groupings — used by Part 2
(`02_popgen_pop_structure_p.Rmd`) and all downstream analyses.

Population structure, FST, LD pruning, AMOVA, and Mantel tests are **not** in
this file — they live in `02_popgen_pop_structure_p.Rmd`. ADMIXTURE,
fineRADstructure, and FEEMS analyses are also **not** in this file — they live
in separate scripts.

---

## Sample groupings

This script produces four nested groupings of the clone-corrected, neutral SNP
dataset, used throughout the rest of the pipeline:

| Dataset  | N individuals | N populations | Definition                                             |
|----------|---------------|----------------|---------------------------------------------------------|
| **All**  | 192           | 13             | Every population, including Edisto                     |
| **GoMa** | 176           | 12             | All Gulf populations regardless of size (no Edisto)     |
| **GoM10**| 158           | 9              | Gulf populations with >10 individuals (excl. WFGB, Bright, Geyer, Edisto) |
| **Sub**  | 174           | 10             | GoM10 + Edisto                                          |

In manuscript text:
**GoM10** = Northern Gulf subset
**Sub** = Gulf-Carolinian subset

All four groupings are built from `neutral_dataset_final` at the end of this
script and carried forward into Part 2.

---

## What this script does

1. **Import & initial filtering** — Reads the raw VCF, converts it to a genind
   object, and removes loci with >10% missing data and individuals with >30%
   missing data.
2. **Clone detection & correction** — Identifies multilocus genotypes via
   `poppr::mlg.filter` on pairwise genetic distance, manually verifies
   ambiguous cross-site clone groups against raw genotype differences, and
   produces the clone-corrected dataset (`snps_mlg`, 192 individuals).
3. **Further locus filtering** — Removes loci with excess observed
   heterozygosity (Ho > 0.5) and loci deviating from Hardy-Weinberg
   equilibrium (`pegas::hw.test`, p < 0.01).
4. **BayeScan outlier removal** — Runs BayeScan externally on the server to
   identify loci under putative selection, then removes high-confidence
   outliers (FST > 0.35, q < 0.05) to produce the final neutral dataset.
5. **Final dataset assembly** — Builds `neutral_dataset_final` (all 13
   populations, 192 individuals, 16,746 loci) and the four population subsets
   described above (All/GoMa/GoM10/Sub).
6. **Summary & QC tables** — Generates a filtering table (individuals/loci/
   populations retained at each step) and a per-site diversity table (Ho, He,
   FIS, allelic richness, % polymorphic, private alleles), and exports
   population-specific metadata files for each dataset subset.

---

## Requirements

R packages used across the script:

```
vcfR, adegenet, poppr, dartR, qgraph, hierfstat, ggplot2, pegas, radiator,
umap, assigner, dplyr, tidyr, tidyverse, ggdendro, pheatmap, viridis, geosphere
```

External software: [BayeScan v2.1](http://cmpg.unibe.ch/software/BayeScan/)
(run outside R; this script prepares BayeScan's input files and processes its
output, but the BayeScan run itself happens on a separate server/environment)

## Required input files

- `swiftia_snps_90_cov_maf_bi.recode.vcf` — raw, biallelic SNP VCF (also
  archived on Figshare)
- `popmap.txt` — tab-delimited sample-to-population map (also archived on
  Figshare)
- `se_metadata` — ⚠ sample metadata (SITE, LAT, LONG, DEPTH per individual) is
  used throughout this script (clone analysis, coordinate correction,
  diversity table) but is not currently loaded from a file anywhere in the
  script itself. Confirm the source filename/object before this is run
  standalone by someone without it already in their environment.

## Outputs

- **Data objects** (`.Rdata`): `all_snps_genind.Rdata`, `snps.Rdata`,
  `snps_mlg.Rdata`, `snps_mlg_he.Rdata`, `neutral_dataset_final.Rdata`,
  `snps_mlg_GoM10_n.Rdata`, `snps_mlg_sub_n.Rdata`,
  `snps_mlg_gom_all_n.Rdata`, `nonclone_subsetted_data.Rdata`, and
  intermediate clone-analysis `.RData` files
- **Tables** (CSV): `filtering_table_FINAL.csv`, `site_diversity_table_v2.csv`,
  per-dataset sample metadata (`se_metadata_neutral_all.csv`, `_sub`,
  `_GoM10`, `_GoMa`), `se_metadata_all_snps.csv`
- **Figures** (PDF): per-site genetic distance histograms and clone networks,
  clone-distance scatterplots, BayeScan outlier plots and heatmap

## Notes

- File paths at the top of the script (`main_dir`, `output_dir`) should be set
  to your local working and output directories before running.
- The BayeScan section requires access to a server/environment with BayeScan
  installed; the R chunks prepare the input files and process the output, but
  the BayeScan run itself is external.
- Two open items are flagged inline and still need author confirmation before
  this is treated as final: the seeded Hardy-Weinberg re-run
  (`snps_mlg.neutral_192`) is not yet wired into the downstream BayeScan/
  outlier-removal steps, and the final BayeScan outlier vector
  (`outliers_to_remove`) was reconstructed to match the documented threshold
  and should be re-run to confirm it evaluates to 18 loci.
