# Stepping-Stone Population Structure and Widespread Clonality of the Mesophotic Octocoral *Swiftia exserta* in the Warm Temperate Northwest Atlantic

R code accompanying:
> Pittoors, N. C., Lopera, L., Vohsen, S. A., Galaska, M. P., West, D.,
> Quattrini, A. M., Bracco, A., & Herrera, S. **Stepping-Stone Population
> Structure and Widespread Clonality of the Mesophotic Octocoral *Swiftia
> exserta* in the Warm Temperate Northwest Atlantic.** In prep.
> *[citation to be updated upon bioRxiv posting]*

This repository contains the R Markdown pipeline used to filter, clone-correct,
and analyze RADseq-derived SNP data for *Swiftia exserta* across mesophotic
banks in the northern Gulf of Mexico and one Atlantic site (Edisto, SC), and
to compare the resulting population structure against a biophysical larval
dispersal model.

---

## Pipeline overview

```
                ┌────────────────────────────────────────┐
raw VCF ───►    │ 01_Data_Setup                          │
                │  Filtering, clone correction (MLG),    │
                │  BayeScan outlier removal               │
                └───────────────────┬─────────────────────┘
                                    │ neutral_dataset_final.Rdata
                                    │ snps_mlg_*.Rdata
                ┌───────────────────┴─────────────────────┐
                ▼                                          ▼
┌────────────────────────────────┐      ┌─────────────────────────────────┐
│ 02_popgen_pop_structure        │      │ 03_ADMIXTURE                    │
│  PCA / DAPC / UMAP             │      │  LD pruning, ancestry           │
│  FST, LD pruning, AMOVA,       │      │  estimation (K=1-10)            │
│  Mantel (IBD) tests            │      │                                  │
└─────────────────────────────────┘      └─────────────────────────────────┘
                ┌───────────────────────────────────────────┐
                ▼                                            
┌────────────────────────────────┐      ┌─────────────────────────────────┐
│ 04_fineRADstructure             │      │ 05_Plot_Physical_Connectivity   │
│  Co-ancestry matrices           │      │  Larval dispersal model         │
│  (RADpainter + fineSTRUCTURE)   │      │  (standalone — no genetic data) │
└─────────────────────────────────┘      └─────────────────────────────────┘
```

Scripts `03` and `04` both branch from `01`'s outputs (and reuse population
grouping objects from `02`). Script `05` is independent — it plots the
biophysical larval dispersal model (Ichthyop/CROCO) on its own and requires
no genetic data.

---

## Repository layout

```
.
├── 01_Data_Setup/
│   ├── 01_Data_Setup_p.Rmd            Raw import, filtering, clone correction, BayeScan
│   └── README.md
├── 02_popgen_pop_structure/
│   ├── 02_popgen_pop_structure_p.Rmd  PCA/DAPC/UMAP, FST, LD pruning, AMOVA, Mantel
│   └── README.md
├── 03_ADMIXTURE/
│   ├── 03_Admixture_p.Rmd             LD-pruned ADMIXTURE (K=1-10), structure plots
│   └── README.md
├── 04_fineRADstructure/
│   ├── 04_fineRADstructure_p.Rmd      Co-ancestry matrices (population + individual)
│   ├── fineRADstructurePlot.R         Plotting script (Malinsky et al. 2018), customized
│   ├── FinestructureLibrary.R         Unmodified upstream function library (Lawson)
│   └── README.md
└── 05_Plot_Physical_Connectivity/
    ├── 05_Plot_Physical_Connectivity.Rmd   Larval dispersal connectivity matrix + map
    └── README.md
```

Each folder has its own README with a full parameter list, required input
files, generated outputs, and step-by-step reproduction instructions —
start there for a given analysis.

**This repository contains scripts only, not the underlying sequence data
or intermediate analysis outputs.** See Data Availability below.

---

## Datasets

All population-level analyses run on some or all of four nested groupings
of the clone-corrected, neutral SNP dataset:

| Dataset | N individuals | N populations | Definition |
|---|---|---|---|
| **All** | 192 | 13 | Every population, including Edisto |
| **GoMa** | 176 | 12 | All Gulf populations regardless of size (no Edisto) |
| **GoM10** | 158 | 9 | Gulf populations with >10 individuals (excl. WFGB, Bright, Geyer, Edisto) |
| **Sub** | 174 | 10 | GoM10 + Edisto |

In manuscript text: GoM10 = Northern Gulf subset; Sub = Gulf-Carolinian subset.

---

## Software requirements

R packages used across the pipeline (see individual READMEs for the exact
subset each script needs):

```
adegenet, poppr, ade4, pegas, umap, vegan, geosphere, hierfstat,
SNPRelate, gdsfmt, dartR, vcfR, ggplot2, gridExtra, patchwork,
tidyverse (dplyr, tidyr, readr, purrr, forcats), viridis, reshape2,
ggdendro, pheatmap, scales, marmap, maps
```

External tools (not R packages):

```
PLINK v1.9
ADMIXTURE
fineRADstructure (RADpainter, fineSTRUCTURE) — Malinsky et al. 2018
```

---

## How to run

Each script's README documents its own reproduction steps in detail. At a
high level:

1. Run `01_Data_Setup_p.Rmd` first — every other script depends on its
   outputs.
2. Run `02_popgen_pop_structure_p.Rmd` for population structure, FST, AMOVA,
   and IBD results.
3. Run `03_Admixture_p.Rmd` and/or `04_fineRADstructure_p.Rmd` in either
   order — both depend on `01` (and reuse population-grouping objects from
   `02`), but not on each other.
4. `05_Plot_Physical_Connectivity.Rmd` can be run independently at any
   point — it only requires the larval connectivity CSVs, not any genetic
   data.

---

## Data availability

Raw VCFs, intermediate ADMIXTURE/fineRADstructure output files, and
generated figures are not tracked in this repository. Data are archived on
Figshare: https://doi.org/10.6084/m9.figshare.33292446 *(embargoed —
will be made public upon manuscript submission)*.

---


## Citation

If you use this code, please cite the manuscript above.

## Contact

Nicole C. Pittoors — <npittoors@mbl.edu>
Santiago Herrera - santiago.herrera@lehigh.edu
