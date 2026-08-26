# fineRADstructure Co-ancestry Analysis — *Swiftia exserta*

---

This directory contains the code required to reproduce the fineRADstructure
co-ancestry analysis for *Swiftia exserta*, using RADpainter and
fineSTRUCTURE (Malinsky et al. 2018, *Mol Biol Evol* 35:1284–1290).

Raw data can be found at https://doi.org/10.6084/m9.figshare.33292446

Two datasets are analyzed:
- **GoM10**: Gulf of Mexico populations with ≥10 individuals (n=158, 9 pops, 16,609 neutral SNPs)
- **Sub**: GoM10 + Edisto (South Atlantic outgroup) (n=174, 10 pops, 16,688 neutral SNPs)

---

## Files

| File | Description |
|------|-------------|
| `04_fineRADstructure_p.Rmd` | Main analysis script. Exports GoM10/Sub VCF subsets for RADpainter, computes and plots population-averaged and individual-level co-ancestry matrices. |
| `fineRADstructurePlot.R` | Plotting script from the fineRADstructure GitHub repo (Malinsky et al. 2018), customized for this project. Must be manually edited before each `source()` call in `04_fineRADstructure_p.Rmd` — see Usage below. |
| `FinestructureLibrary.R` | Unmodified upstream R function library from the fineSTRUCTURE distribution (Daniel Lawson). Required by `fineRADstructurePlot.R`; not edited for this project. |

Everything else referenced by `04_fineRADstructure_p.Rmd` — the input VCF,
`popmap.txt`, and the RADpainter/fineSTRUCTURE intermediate and figure
files it generates — is produced by running the pipeline yourself rather than shipped in this folder.

### Figures produced by this script

| File | Description |
|------|-------------|
| `figures/Fig_coancestry_GoM10_population_matrix.pdf/png` | Population-averaged co-ancestry matrix, GoM10 (9×9), west-to-east order. Diagonal excluded. |
| `figures/Fig_coancestry_Sub_individual_heatmap.pdf/png` | Individual co-ancestry heatmap, Sub (174×174), west-to-east order. Color scale restricted to 1st–99th percentile off-diagonal range. |
| `figures/Fig_coancestry_Sub_population_matrix.pdf/png` | Population-averaged co-ancestry matrix, Sub (10×10), west-to-east order. Diagonal excluded. |
| `figures/fineRAD_ind_heatmap_GoM10_clustered_annotated_final.pdf` | Individual co-ancestry heatmap, GoM10, ordered by the fineRADstructure hierarchical clustering tree rather than geography. |
| `figures/fineRAD_ind_heatmap_sub_clustered_annotated_final.pdf` | Individual co-ancestry heatmap, Sub, tree-ordered. |

---

## Software Requirements

### fineRADstructure binaries (compile from source)

```bash
# Clone repository
git clone https://github.com/millanek/fineRADstructure.git
cd fineRADstructure

# On Apple Silicon (M1/M2/M3) Mac: patch SSE flags before compiling
sed -i '' 's/-mfpmath=sse -msse -msse2 //g' Makefile

# Compile (requires GSL: brew install gsl)
make CXXFLAGS="-O3 -Wall -funroll-loops -fomit-frame-pointer \
  -ftree-vectorize -funsafe-math-optimizations \
  -I/opt/homebrew/Cellar/gsl/2.8/include" \
  LDFLAGS="-L/opt/homebrew/Cellar/gsl/2.8/lib" \
  LIBS="-lgsl -lgslcblas -lz"
```

Produces two binaries: `RADpainter` and `finestructure` in the
`fineRADstructure/` directory.

**On Intel Mac or Linux:** standard `make` should work after
`brew install gsl` or `apt install libgsl-dev`.

### Create symlink to avoid spaces-in-path issues (Mac)

```bash
ln -s "/path/to/Swiftia_PopGen2/fineRADstructure" ~/fineRAD
```

### R packages required

```r
install.packages(c("vcfR", "tidyverse", "reshape2", "viridis", "pheatmap", "scales", "XML"))
# XML is required by FinestructureLibrary.R
```

---

## Usage — Step-by-step reproduction

### Step 1: Generate input VCF files

Run Section 1 of `04_fineRADstructure_p.Rmd`. Requires
`swiftia_snps_90_cov_maf_bi.recode.vcf`, `pop_code_GoM10.Rdata`, and
`pop_code_sub.Rdata` (not included in this repo — see Overview above).

**Note:** `write.vcf()` from vcfR always produces gzipped output. Use
`vcfR::write.vcf()` explicitly to avoid conflict with `pegas::write.vcf()`,
which expects a different object class.

### Step 2: Run fineRADstructure pipeline (Terminal)

```bash
cd "/path/to/Swiftia_PopGen2"

# --- GoM10 ---
~/fineRAD/RADpainter hapsFromVCF swiftia_snps_GoM10.vcf.gz > swiftia_GoM10_haps.txt

Rscript ~/fineRAD/sampleLD.R -s 1 -n 500 \
  swiftia_GoM10_haps.txt swiftia_GoM10_haps_reordered.txt

~/fineRAD/RADpainter paint swiftia_GoM10_haps_reordered.txt
mv swiftia_GoM10_haps_reordered_chunks.out swiftia_GoM10_chunks.out

~/fineRAD/finestructure -x 100000 -y 100000 -z 1000 \
  swiftia_GoM10_chunks.out swiftia_GoM10_chunks.mcmc.xml

~/fineRAD/finestructure -m T -x 10000 \
  swiftia_GoM10_chunks.out \
  swiftia_GoM10_chunks.mcmc.xml \
  swiftia_GoM10_chunks.mcmcTree.xml

# --- Sub ---
~/fineRAD/RADpainter hapsFromVCF swiftia_snps_sub.vcf.gz > swiftia_sub_haps.txt

Rscript ~/fineRAD/sampleLD.R -s 1 -n 500 \
  swiftia_sub_haps.txt swiftia_sub_haps_reordered.txt

~/fineRAD/RADpainter paint swiftia_sub_haps_reordered.txt
mv swiftia_sub_haps_reordered_chunks.out swiftia_sub_chunks.out

~/fineRAD/finestructure -x 100000 -y 100000 -z 1000 \
  swiftia_sub_chunks.out swiftia_sub_chunks.mcmc.xml

~/fineRAD/finestructure -m T -x 10000 \
  swiftia_sub_chunks.out \
  swiftia_sub_chunks.mcmc.xml \
  swiftia_sub_chunks.mcmcTree.xml
```

### Step 3: Edit `fineRADstructurePlot.R` before each `source()` call

Open `fineRADstructurePlot.R` and set these lines for **GoM10**:

```r
chunkfile    <- "swiftia_GoM10_chunks.out"
mcmcfile     <- "swiftia_GoM10_chunks.mcmc.xml"
treefile     <- "swiftia_GoM10_chunks.mcmcTree.xml"
plotsFolder  <- "/path/to/Swiftia_PopGen2/figures/"
analysisName <- "swiftia_GoM10"; maxIndv <- 1000; maxPop <- 1000
```

For **Sub**, change to:

```r
chunkfile    <- "swiftia_sub_chunks.out"
mcmcfile     <- "swiftia_sub_chunks.mcmc.xml"
treefile     <- "swiftia_sub_chunks.mcmcTree.xml"
analysisName <- "swiftia_sub"; maxIndv <- 1000; maxPop <- 1000
```

`04_fineRADstructure_p.Rmd` calls `source("fineRADstructurePlot.R")`
multiple times (once per dataset, per figure section) — each call picks up
whatever config is currently active in the file, so re-editing and
re-sourcing between GoM10 and Sub sections is required, not optional.

### Step 4: Run the R Markdown pipeline

Run `04_fineRADstructure_p.Rmd` section by section, alternating with Step 3
edits as indicated by the "BEFORE RUNNING" notes in the script.

**Critical known issue — Edisto sample name mangling:** R's `read.table()`
converts sample IDs beginning with digits and containing hyphens (e.g.
`14-v-18-1-001`) to R-safe names (e.g. `X14.v.18.1.001`). Every Sub-dataset
chunk in `04_fineRADstructure_p.Rmd` corrects this immediately after
`source("fineRADstructurePlot.R")`:

```r
rownames(dataraw) <- gsub("^X", "", gsub("\\.", "-", rownames(dataraw)))
colnames(dataraw) <- gsub("^X", "", gsub("\\.", "-", colnames(dataraw)))
```

This fix is already incorporated into the script. Do not remove it.

---

## Parameter decisions

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `sampleLD.R -n 500` | 500 loci sampled for LD estimation | Standard recommendation for datasets without genomic coordinates |
| `sampleLD.R -s 1` | Random seed = 1 | Set for reproducibility |
| fineSTRUCTURE `-x 100000` | 100,000 MCMC iterations | Standard for sample sizes of 158–174 individuals |
| fineSTRUCTURE `-y 100000` | 100,000 burn-in iterations | Equal burn-in to sampling, conservative |
| fineSTRUCTURE `-z 1000` | Thinning interval | Retains 100 samples from MCMC chain |
| Color scale (Sub individual heatmap) | 1st–99th percentile | Prevents diagonal dominance in individual heatmap |
| Diagonal exclusion | NA | Within-population self-comparisons excluded from population matrices |

---

## Citation

Malinsky, M., Trucchi, E., Lawson, D.J. & Falush, D. (2018) RADpainter and fineRADstructure: Population Inference from RADseq Data. *Molecular Biology and Evolution* 35:1284–1290. https://doi.org/10.1093/molbev/msy023

Lawson, D.J., Hellenthal, G., Myers, S. & Falush, D. (2012) Inference of population structure using dense haplotype data. *PLOS Genetics* 8:e1002453.
