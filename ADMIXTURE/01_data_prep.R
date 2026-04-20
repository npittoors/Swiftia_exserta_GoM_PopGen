# =============================================================================
# 01_data_prep.R
# Swiftia exserta population genomics — ADMIXTURE data preparation
# Converts clone-corrected neutral genind objects to PLINK format for ADMIXTURE
#
# Input:  neutral_dataset_final.Rdata  (192 ind, ~18,616 loci)
#         snps_mlg_GoM10_n.Rdata       (158 ind, GoM only)
#         snps_mlg_gom.Rdata           (176 ind, GoMa — includes re-added singletons)
# Output: PLINK .ped/.map files for each dataset (convert to .bed with PLINK v1.9)
#
# Run after: clone correction, Ho filter, HWE filter, BayeScan outlier removal
# =============================================================================

library(dartR)
library(adegenet)

# ---- Set working directory --------------------------------------------------
setwd("/set/to/working/directory")

# =============================================================================
# DATASET 1: Sub (Gulf-Atlantic, 174 ind)
# =============================================================================

load("neutral_dataset_final.Rdata")  # loads: neutral_dataset_byscn

# Confirm composition
table(pop(neutral_dataset_byscn))
nInd(neutral_dataset_byscn)  # 192 (full neutral dataset; Sub subsets to 174)

# Convert genind → genlight
gl_sub <- gi2gl(neutral_dataset_byscn)

# Add required chromosome and position metadata
gl_sub$chromosome <- as.factor(rep("1", nLoc(gl_sub)))
gl_sub$position   <- as.integer(1:nLoc(gl_sub))

# Verify
nInd(gl_sub)   # 192
nLoc(gl_sub)   # ~18,616

# Export to PLINK .ped/.map
gl2plink(gl_sub, outfile = "swiftia_neutral", outpath = "./")


# =============================================================================
# DATASET 2: GoM10 (Gulf of Mexico only, 158 ind)
# =============================================================================

load("snps_mlg_GoM10_n.Rdata")  # loads: snps_mlg_GoM10_n

# Confirm composition
table(pop(snps_mlg_GoM10_n))
nInd(snps_mlg_GoM10_n)  # 158
nLoc(snps_mlg_GoM10_n)  # ~16,609–18,616

# Convert genind → genlight
gl_GoM10 <- gi2gl(snps_mlg_GoM10_n)
gl_GoM10$chromosome <- as.factor(rep("1", nLoc(gl_GoM10)))
gl_GoM10$position   <- as.integer(1:nLoc(gl_GoM10))

# Verify
nInd(gl_GoM10)  # 158
nLoc(gl_GoM10)

# Export
gl2plink(gl_GoM10, outfile = "swiftia_GoM10", outpath = "./")


# =============================================================================
# DATASET 3: GoMa (GoM with non-clonal singletons re-added, 176 ind)
# =============================================================================

load("snps_mlg_gom.Rdata")  # loads: snps_mlg_gom_all_n

# Confirm composition
table(pop(snps_mlg_gom_all_n))
nInd(snps_mlg_gom_all_n)  # 176
nLoc(snps_mlg_gom_all_n)

# Convert genind → genlight
gl_GoMa <- gi2gl(snps_mlg_gom_all_n)
gl_GoMa$chromosome <- as.factor(rep("1", nLoc(gl_GoMa)))
gl_GoMa$position   <- as.integer(1:nLoc(gl_GoMa))

# Verify
nInd(gl_GoMa)  # 176
nLoc(gl_GoMa)

# Export
gl2plink(gl_GoMa, outfile = "swiftia_GoMa", outpath = "./")


# =============================================================================
# PLINK CONVERSION (run in terminal after this script)
# Converts .ped/.map to binary .bed/.bim/.fam required by ADMIXTURE
# =============================================================================
# plink --file swiftia_neutral --make-bed --out swiftia_neutral --allow-extra-chr
# plink --file swiftia_GoM10   --make-bed --out swiftia_GoM10   --allow-extra-chr
# plink --file swiftia_GoMa    --make-bed --out swiftia_GoMa    --allow-extra-chr


# =============================================================================
# BUILD POPULATION METADATA FOR DOWNSTREAM PLOTTING
# Run after loading the genind objects above
# =============================================================================

# --- Sub ---
pop_data_sub <- data.frame(
  individual = indNames(neutral_dataset_byscn),
  population = as.character(pop(neutral_dataset_byscn)),
  stringsAsFactors = FALSE
)

# --- GoM10 ---
pop_data_GoM10 <- data.frame(
  individual = indNames(snps_mlg_GoM10_n),
  population = as.character(pop(snps_mlg_GoM10_n)),
  stringsAsFactors = FALSE
)

# --- GoMa ---
pop_data_GoMa <- data.frame(
  individual = indNames(snps_mlg_gom_all_n),
  population = as.character(pop(snps_mlg_gom_all_n)),
  stringsAsFactors = FALSE
)

# Confirm individual order matches metadata before ADMIXTURE plotting:
# all(pop_data_GoM10$individual == rownames(se_metadata_neutral_GoM10))  # must be TRUE
# all(pop_data_GoMa$individual  == rownames(se_metadata_neutral_GoMa))   # must be TRUE

# Save pop_data objects for use in 02_process_results.R
save(pop_data_sub, pop_data_GoM10, pop_data_GoMa,
     file = "pop_data_admixture.Rdata")

cat("Data prep complete. Run PLINK conversion, then shell scripts, then 02_process_results.R\n")
