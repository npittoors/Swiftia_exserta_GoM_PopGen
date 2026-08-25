**`01_Data_Setup.Rmd`** — Part 1 of the pipeline (this document). Raw SNP
  import and quality filtering, clone detection and correction, further locus
  filtering (heterozygosity, Hardy-Weinberg equilibrium), BayeScan outlier
  removal, and assembly of the final neutral dataset and its population
  subsets.


## Data (available on Figshare https://doi.org/10.6084/m9.figshare.33292446)

The following files accompany this repository on Figshare and are the inputs
and key intermediate objects for `01_Data_Setup.Rmd`:

1. **`swiftia_snps_90_cov_maf_bi.recode.vcf`**
   Raw, unfiltered SNP dataset for *Swiftia exserta* in VCF format, called from
   RADseq/ddRAD data with a 90% coverage threshold, minor allele frequency
   filter, and restricted to biallelic sites. Contains 311 individuals prior to
   any missing-data filtering, clone correction, or outlier removal.

2. **`popmap.txt`**
   Tab-delimited population map assigning each sequenced individual to its
   sampling population/site. Two columns: `INDIVIDUALS` (sample ID, matching
   genind object individual names) and `STRATA` (population/site assignment).
   Used to assign population membership when building genind objects from the
   raw VCF.

3. **`all_snps_genind.Rdata`**
   R genind object (adegenet) generated directly from
   `swiftia_snps_90_cov_maf_bi.recode.vcf` via `vcfR::vcfR2genind()`, with no
   additional filtering applied. Contains 311 individuals and represents the
   full, unfiltered SNP set. Load in R with `load("all_snps_genind.Rdata")`;
   the object name is `genind_all`.

4. **`se_metadata_all_snps.csv`**
   Sample metadata for the 311 individuals in `all_snps_genind.Rdata` and
   `swiftia_snps_90_cov_maf_bi.recode.vcf`. Includes sample ID, collection
   site, depth, and geographic coordinates for each *Swiftia exserta* colony.

5. **`neutral_dataset_final.Rdata`**
   R genind object (adegenet) containing the final neutral SNP dataset used in
   population genomic analyses. Derived from the raw SNP dataset (above) after
   loci/individual missing-data filtering, clone correction, Hardy-Weinberg
   equilibrium filtering, and removal of BayeScan-identified FST outlier loci.
   Contains 192 individuals and 16,746 neutral loci across 13 populations
   spanning the northern Gulf of Mexico and one Atlantic site (Edisto, SC).
   Load in R with `load("neutral_dataset_final.Rdata")`; the object name is
   `neutral_dataset_final`.

## What `01_Data_Setup.Rmd` does

1. **Import & initial filtering** — Reads the raw VCF, converts to a genind
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
   populations, 192 individuals, 16,746 loci) and three population subsets
   used throughout Part 2:
   - **All** — all 13 populations
   - **Sub** (Gulf–Atlantic) — 10 populations, including Edisto
   - **GoM10** — 9 Gulf of Mexico populations with ≥10 individuals
   - **GoMa** — 12 Gulf of Mexico populations (all sizes), Edisto excluded
6. **Summary & QC tables** — Generates a filtering table (individuals/loci/
   populations retained at each step) and a per-site diversity table (Ho, He,
   FIS, allelic richness, % polymorphic, private alleles), and exports
   population-specific metadata files for each dataset subset.

## Requirements

**R packages:** vcfR, adegenet, poppr, dartR, qgraph, hierfstat, ggplot2,
pegas, radiator, umap, assigner, dplyr, tidyr, tidyverse, ggdendro, pheatmap,
viridis, geosphere

**External software:** [BayeScan v2.1](http://cmpg.unibe.ch/software/BayeScan/)
(run outside R; see the BayeScan section of `01_Data_Setup.Rmd` for file
preparation and invocation)

## Notes

- File paths at the top of each script (`main_dir`, `output_dir`) should be set
  to your local working and output directories before running.
- The BayeScan section requires access to a server/environment with BayeScan
  installed; the R chunks prepare the input files and process the output, but
  the BayeScan run itself is external.
