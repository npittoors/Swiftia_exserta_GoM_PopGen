# Swiftia exserta Population Genomics — Part 1: Data Import & Neutral Dataset Setup

This script (`popgen_01_data_import_processed.R`) is the first stage of the analysis
pipeline for population genomic structure of the mesophotic octocoral
*Swiftia exserta* across the northern Gulf of Mexico and one Atlantic site (Edisto, SC).
It takes raw SNP genotype calls through filtering, clone identification and correction,
Hardy–Weinberg and heterozygosity filtering, BayeScan outlier removal, and finally
builds the neutral datasets and population-map objects used by all downstream scripts.

All file paths in the script are replaced with the placeholder `/File/Path`. Set
`main_dir` and `output_dir` at the top of the script to your own working directory
before running.

---

## Pipeline overview

```
raw VCF (288 ind)
   └─ missingness filtering (loci ≤10% missing, individuals ≤30% missing)
       └─ clone identification & correction (MLG collapse, manual group review)
           └─ heterozygosity (Ho) filter
               └─ HWE filter
                   └─ BayeScan outlier removal
                       └─ neutral_dataset_final (192 ind)  ──►  4 subset datasets
```

---

## 1. Required input files (NOT produced by this script)

These must be present in the working directory before running. **Raw sequence data
and the VCF are archived separately** (see "Raw data" below); everything else here is
a small text/metadata file you should include in the repository alongside the code.

| File | Description |
|------|-------------|
| `swiftia_snps_90_cov_maf_bi.recode.vcf` | Filtered biallelic SNP genotype calls (90% coverage, MAF-filtered). Input to the whole pipeline. **Archive with raw data, not in repo.** |
| `popmap.txt` | Tab-delimited population map: column 1 = individual ID, column 2 = site/strata. No header. |
| `se_metadata` (R object / source table) | Sample metadata: site (`SITE`), depth (`DEPTH`), coordinates (`LAT`, `LONG`), used for clone spatial analysis. Include the source table (e.g. `se_metadata.csv`) the object is built from. |



---

## 2. Software / environment

R packages: `vcfR`, `dplyr`, `adegenet`, `poppr`, `dartR`, `qgraph`, `hierfstat`,
`ggplot2`, `pegas`, `radiator`, `umap`, `assigner`, `tidyr`, `ggdendro`, `pheatmap`,
`tidyverse`, `geosphere`, `viridis`, `gridExtra`, `reshape2`, `SNPRelate`.

External tools: **BayeScan v2.1** (run on a compute server — the bash chunk shows the
file-formatting and run command; adjust the executable path to your environment).

---

## 3. R data objects produced

### Intermediate objects (full dataset, pre-subsetting)

| Object / file | Contents |
|---------------|----------|
| `snps.Rdata` (`snps`) | genind object after missingness filtering, with population assignments (288 ind, 13 pops, includes clones). |
| `pc.Rdata` (`pc`) | genclone object; MLGs collapsed at diss.dist threshold 2000. |
| `snps_mlg.Rdata` (`snps_mlg`) | Clone-corrected genind (one individual per confirmed clone group; manual handling of groups 182, 221, 195). |
| `snps_mlg_he.Rdata` (`snps_mlg_he`) | After heterozygosity (Ho) filter. |
| `snps_mlg_he.hw.Rdata` (`snps_mlg_he.hw`) | After HWE filter. |
| `neutral_dataset_byscn.Rdata` (`neutral_dataset_byscn`) | After initial BayeScan outlier removal. |
| `neutral_dataset_final.Rdata` (`neutral_dataset_final`) | **Final neutral dataset, 192 ind** (non-clonal singletons re-added). Master object for all downstream analyses. |

### Clone-analysis objects

| File | Contents |
|------|----------|
| `initial_clone_analysis_results.Rdata` | `clone_groups_initial`, `clone_details_initial`, `clone_stats_initial`. |
| `clone_analysis_final.RData` | `final_clone_details`, `all_results`, `site_summary`, `snps_clones`. |
| `clean_clone_analysis_final.RData` | Cleaned `snps_clones` (cross-site clone groups resolved). |
| `distgenDISS_*.Rdata` | Per-site pairwise allelic-difference distance matrices (one per site: AAR, ALD, BOU, BRI, DIA, EDI, EFGB, ELV, FTR, GEY, MCG, RTR, WFGB). |
| `distgenEUCL_se.Rdata`, `distgenDISS_se.Rdata` | Full-dataset Euclidean / dissimilarity distance matrices. |
| `basicstat_all.Rdata` | Basic population statistics (hierfstat). |

### Final subset datasets (the deliverables for downstream scripts)

All bundled in **`nonclone_subsetted_data.Rdata`** (and individually saved):

| Dataset | n | genind object | pop map | metadata | pop order | W→E order |
|---------|---|---------------|---------|----------|-----------|-----------|
| **All** (all sites) | 192 | `neutral_dataset_final` | `pop_code_all` | `se_metadata_neutral_all` | `pop_order_all` | `WE_order_all` |
| **Sub** (Edisto + GoM pops >10 ind) | 174 | `snps_mlg_sub_n` | `pop_code_sub` | `se_metadata_neutral_sub` | `pop_order_sub` | `WE_order_sub` |
| **GoMa** (all GoM, no Edisto) | 176 | `snps_mlg_gom_all_n` | `pop_code_GoMa` | `se_metadata_neutral_GoMa` | `pop_order_GoMa` | `WE_order_GoMa` |
| **GoM10** (GoM pops >10 ind, no Edisto) | 158 | `snps_mlg_GoM10_n` | `pop_code_GoM10` | `se_metadata_neutral_GoM10` | `pop_order_GoM10` | `WE_order_GoM10` |

Helper function `create_pop_code()` (builds the pop_code objects) is also saved in the bundle.

---

## 4. BayeScan input/output files

| File | Description |
|------|-------------|
| `swiftia_snps_bayescan_neutral.txt` | BayeScan-format genotype input written from R. |
| `swiftia_snps_bayescan_neutral_comprehensive.txt` | Reformatted input with complete locus coverage (server bash step). |
| `swiftia_population_assignments.txt`, `swiftia_population_assignments_precise.txt` | Population-prior files. |
| `swiftia_bayescan_output_fst.txt` | BayeScan FST results (read back into R for outlier ID). |
| `swiftia_bayescan_output_Verif.txt` | BayeScan verification/log output. |

---

## 5. CSV / text outputs

| File | Contents |
|------|----------|
| `all_sites_clone_analysis.csv` | Clone groups detected per site. |
| `samples_with_clones.txt` | Samples belonging to clone groups. |
| `Swiftia_popgen_all_samples_postfilter.txt` | All retained samples after filtering. |
| `clone_samples_with_ID.txt` | Clone samples with group IDs. |
| `outlier_summary.csv` | BayeScan outlier loci. |
| `bayescan_summary_stats.csv` | Summary FST/alpha stats across loci. |
| `top_30_outliers.csv`, `high_fst_loci.csv` | Outlier locus subsets. |

---

## 6. Figures produced (PDF)

`genetic_distances.pdf`, `genetic_distance_histograms_allsites.pdf`,
`genetic_distance_histograms_ggplot.pdf`, `clone_network_<site>.pdf`,
`clone_distances_depth.pdf`, `clone_distances_histo.pdf`,
`BayeScan_Results_Outliers.pdf`.


**Raw data** link the DOI/accessions here once completed

---

## Population codes & ordering

Sites west → east: **WFGB, EFGB, Bright, Geyer, Elvers, McGrail, Bouma, Alderdice,
Diaphus, FTR, AAR, RTR, Edisto.**

Subset definitions:
- **GoM10** = GoM populations with >10 individuals, excludes WFGB, Bright, Geyer, Edisto.
- **Sub** = GoM10 + Edisto.
- **GoMa** = all GoM populations, excludes Edisto only.
