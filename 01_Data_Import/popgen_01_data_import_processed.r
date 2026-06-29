---
title: "popgen_01_data_import"
author: "Nicole Pittoors"
date: "2026-03-30"
output: html_document
---
Set Working Directory
setwd('/File/Path')

```{r}
main_dir <- "/File/Path"
output_dir <- "/File/Path/Feb18Output"
setwd(main_dir)
```


```{r, setup, include=FALSE}
knitr::opts_knit$set(root.dir = '/File/Path')
```


#Load Libraries
```{r}
library(vcfR)
library(dplyr)
library(adegenet)
library(poppr)
library(dartR)
library(qgraph)
library(hierfstat)
library(ggplot2)
library(pegas)
library(radiator)
library(umap)
library(assigner)
library(tidyr)
library(ggdendro)
library(pheatmap)
library(tidyverse)
library(geosphere)
library(viridis)
```

```{r}
# Define colors and sites
# Matches site for ARMS paper as best as possible
c13_sites <- c("Alderdice" = "#60714E",      
               "Geyer" = "#015951",
               "Elvers" = "#3B8C6E",
               "Bright" = "#B0DBBB",
               "Bouma" = "#5BBFB0",
               "Diaphus" = "#9DC9D3",
               "FTR" = "#408BA8",
               "WFGB" = "#0D3D59",
               "EFGB" = "#A8B1D6",
               "AAR" = "#FEBC2C",
               "McGrail" = "#FC8D2B",
               "RTR" = "#F1531B",
               "Edisto"= "#802A5A")

c_13sites_WE <- c("WFGB" = "#F1531B",  
  "EFGB" = "#FC8D2B",  
  "Bright" = "#FEBC2C",  
  "Geyer" = "#802A5A",  
  "McGrail" = "#A8B1D6",  
  "Elvers" = "#9DC9D3",  
  "Bouma" = "#408BA8",  
  "Alderdice" = "#0D3D59",  
  "Diaphus" = "#5BBFB0",  
  "FTR" = "#B0DBBB",  
  "AAR" = "#3B8C6E",  
  "RTR" = "#015951",  
  "Edisto" = "#60714E")

site_colors <- c(
  "WFGB" = "#700500",
  "EFGB" = "#FF2EBF",
  "Bright" = "#FF2E00",
  "Geyer" = "#FF7514",
  "McGrail" = "#F78C61",
  "Elvers" = "#FFB812",
  "Bouma" = "#6F8E00",
  "Alderdice" = "#1B6758",
  "Diaphus" = "#B0FFAD",
  "FTR" = "#45C2E3",
  "AAR" = "#000872",
  "RTR" = "#F5C7FF",
  "Edisto" = "#6E0A7C"
)

# Create a vector of site names in west-to-east order
W2E <- c('WFGB','EFGB','Bright','Geyer','McGrail','Elvers','Bouma','Alderdice','Diaphus','FTR','AAR','RTR','Edisto')

# Create color mapping for sites
sitecol <- data.frame(
  W2E = names(c13_sites),
  mypal = unname(c13_sites),
  stringsAsFactors = FALSE
)

# Arrange colors in west-to-east order
myCol.df <- sitecol[match(W2E, sitecol$W2E), ]
myCol <- myCol.df$mypal
```

=======================================================================
---------------------------Initial Filtering---------------------------
========================================================================
```{r}
vcf <- vcfR::read.vcfR("swiftia_snps_90_cov_maf_bi.recode.vcf")
vcf

# Transform vcf to genind
genind_all <- vcfR::vcfR2genind(vcf) #all SNP

# Remove loci with too many missing (10% allowable missing data in this case)
filt.loci <- missingno(genind_all, type = "loci", cutoff = 0.10, quiet = FALSE, freq = FALSE)
filt.loci

# Remove individuals with too many missing (30% allowable missing data in this case)
filt.loci_ind <- missingno(filt.loci, type =  "geno", cutoff = 0.3, quiet = FALSE, freq = FALSE)

# Load pop info
pop <- read.table("popmap.txt", header=F, sep="\t", stringsAsFactors = TRUE)
colnames(pop)<-c('INDIVIDUALS','STRATA')
rownames(pop)<-pop$INDIVIDUALS

# Check for differences bewteen snp file and popmap file
ind_snps<-as.vector(indNames(filt.loci_ind), mode = 'character')
ind_pop<-as.vector(pop[,1], mode = 'character')
print(ind_snps %in% ind_pop)
setdiff(ind_snps,ind_pop)

# Exclude Edisto
#pop<-filter(pop, STRATA != 'Edisto')

#Add population info
ind<-pop[indNames(filt.loci_ind),]
ind<-na.omit(ind)
snps <- filt.loci_ind[ind$INDIVIDUALS]
snps@pop <- ind$STRATA
summary(snps@pop)
# AAR Alderdice     Bouma    Bright   Diaphus    Edisto      EFGB    Elvers       FTR     Geyer   McGrail       RTR 
# 14        31        19        29        28        16        20        29        19        26        25        26 
# WFGB 
# 6
save(snps,file='snps.Rdata')
```

========================================================================
----------------------------Clone Analysis-----------------------------
========================================================================
TO DO: Measure relatedness
Install Xcode and the recommended development tools for macOS systems. https://cran.r-project.org/bin/macosx/tools/
In R Install related https://rdrr.io/rforge/related/
Use the function coancestry with the ritland estimator as recommended by Attard et al 2018 Mol Ecol Res
Also check out the R package SNPRelate (already installed) for pca based relatedness
"snps" have been filtered and also contain clones
```{r}
load('snps.Rdata')
snps.gl<-gi2gl(snps)
gl2gds(snps.gl, outfile = "snps.gds", outpath = '/Users/nicolepittoors/Documents/Research/Swiftia_PopGen')
#The total number of samples: 288 
#The total number of SNPs: 23869 
#SNP genotypes are stored in individual-major mode (SNP X Sample).
#The number of valid samples: 288 
#The number of biallelic unique SNPs: 0 
```

# Check visually for clones
First pass look before more stringent classification of clones
```{r}
pdf('./genetic_distances.pdf')
#AAR
#pairwise number of allelic differences between two individuals (diss.dist {poppr})
distgenDISS_AAR <- diss.dist(snps[pop='AAR'], percent = FALSE, mat = FALSE)
hist(distgenDISS_AAR)
save(distgenDISS_AAR, file = 'distgenDISS_AAR.Rdata')
#Visualize genetic distances  as network to detect clusters of highly related individuals (potential clones)
qgraph(1/distgenDISS_AAR, layout='spring', vsize=5, minimum = 0.00020)
#Alderdice
distgenDISS_ALD <- diss.dist(snps[pop='Alderdice'], percent = FALSE, mat = FALSE)
hist(distgenDISS_ALD)
save(distgenDISS_ALD, file = 'distgenDISS_ALD.Rdata')
qgraph(1/distgenDISS_ALD, layout='spring', vsize=5, minimum = 0.00020)
#Bouma
distgenDISS_BOU <- diss.dist(snps[pop='Bouma'], percent = FALSE, mat = FALSE)
hist(distgenDISS_BOU)
save(distgenDISS_BOU, file = 'distgenDISS_BOU.Rdata')
qgraph(1/distgenDISS_BOU, layout='spring', vsize=5, minimum = 0.00020)
#Bright
distgenDISS_BRI <- diss.dist(snps[pop='Bright'], percent = FALSE, mat = FALSE)
hist(distgenDISS_BRI)
save(distgenDISS_BRI, file = 'distgenDISS_BRI.Rdata')
qgraph(1/distgenDISS_BRI, layout='spring', vsize=5, minimum = 0.00020)
#Diaphus
distgenDISS_DIA <- diss.dist(snps[pop='Diaphus'], percent = FALSE, mat = FALSE)
hist(distgenDISS_DIA)
save(distgenDISS_DIA, file = 'distgenDISS_DIA.Rdata')
qgraph(1/distgenDISS_DIA, layout='spring', vsize=5, minimum = 0.00020)
#EFGB
distgenDISS_EFGB <- diss.dist(snps[pop='EFGB'], percent = FALSE, mat = FALSE)
hist(distgenDISS_EFGB)
save(distgenDISS_EFGB, file = 'distgenDISS_EFGB.Rdata')
qgraph(1/distgenDISS_EFGB, layout='spring', vsize=5, minimum = 0.00020)
#Edisto
distgenDISS_EDI <- diss.dist(snps[pop='Edisto'], percent = FALSE, mat = FALSE)
hist(distgenDISS_EDI)
save(distgenDISS_EDI, file = 'distgenDISS_EDI.Rdata')
qgraph(1/distgenDISS_EDI, layout='spring', vsize=5, minimum = 0.00020)
#Elvers
distgenDISS_ELV <- diss.dist(snps[pop='Elvers'], percent = FALSE, mat = FALSE)
hist(distgenDISS_ELV)
save(distgenDISS_ELV, file = 'distgenDISS_ELV.Rdata')
qgraph(1/distgenDISS_ELV, layout='spring', vsize=5, minimum = 0.00020)
#FTR
distgenDISS_FTR <- diss.dist(snps[pop='FTR'], percent = FALSE, mat = FALSE)
hist(distgenDISS_FTR)
save(distgenDISS_FTR, file = 'distgenDISS_FTR.Rdata')
qgraph(1/distgenDISS_FTR, layout='spring', vsize=5, minimum = 0.00020)
#Geyer
distgenDISS_GEY <- diss.dist(snps[pop='Geyer'], percent = FALSE, mat = FALSE)
hist(distgenDISS_GEY)
save(distgenDISS_GEY, file = 'distgenDISS_GEY.Rdata')
qgraph(1/distgenDISS_GEY, layout='spring', vsize=5, minimum = 0.00020)
#McGrail
distgenDISS_MCG <- diss.dist(snps[pop='McGrail'], percent = FALSE, mat = FALSE)
hist(distgenDISS_MCG)
save(distgenDISS_MCG, file = 'distgenDISS_MCG.Rdata')
qgraph(1/distgenDISS_MCG, layout='spring', vsize=5, minimum = 0.00020)
#RTR
distgenDISS_RTR <- diss.dist(snps[pop='RTR'], percent = FALSE, mat = FALSE)
hist(distgenDISS_RTR)
save(distgenDISS_RTR, file = 'distgenDISS_RTR.Rdata')
qgraph(1/distgenDISS_RTR, layout='spring', vsize=5, minimum = 0.00020)
#WFGB
distgenDISS_WFGB <- diss.dist(snps[pop='WFGB'], percent = FALSE, mat = FALSE)
hist(distgenDISS_WFGB)
save(distgenDISS_WFGB, file = 'distgenDISS_WFGB.Rdata')
qgraph(1/distgenDISS_WFGB, layout='spring', vsize=5, minimum = 0.00020)
dev.off()
```

# PLot all Histograms in one figure
```{r}
# Method 1: Using base R with par()
# This sets up a 4x3 grid of plots (adjust rows and columns as needed)
pdf("genetic_distance_histograms_allsites.pdf", width=12, height=10)  # Create a PDF file
par(mfrow=c(5,3), mar=c(4,4,3,1))  # 4 rows, 3 columns, adjust margins

# Plot each histogram with its title
hist(distgenDISS_AAR, main="AAR", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_ALD, main="Alderdice", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_BOU, main="Bouma", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_BRI, main="Bright", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_DIA, main="Diaphus", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_ELV, main="Elvers", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_FTR, main="FTR", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_GEY, main="Geyer", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_MCG, main="McGrail", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_RTR, main="RTR", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_EFGB, main="EFGB", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_WFGB, main="WFGB", xlab="Genetic Distance", col="lightblue")
hist(distgenDISS_EDI, main="Edisto", xlab="Genetic Distance", col="lightblue")
# Reset the plot parameters
#par(mfrow=c(1,1))

dev.off()

# Method 2: Using ggplot2 for more customization (recommended)
# Install required packages if not already installed
# install.packages(c("ggplot2", "gridExtra", "reshape2"))

library(ggplot2)
library(gridExtra)
library(reshape2)

# Create a function to convert distance matrix to dataframe for ggplot
matrix_to_df <- function(dist_matrix, site_name) {
  # Convert distance matrix to vector
  dist_vec <- as.vector(as.matrix(dist_matrix))
  # Remove diagonal elements (zeros) if it's a square matrix
  if(is.matrix(dist_matrix)) {
    dist_vec <- dist_vec[dist_vec > 0]
  }
  return(data.frame(site = site_name, distance = dist_vec))
}

# Combine all distance data into one dataframe
all_distances <- rbind(
  matrix_to_df(distgenDISS_AAR, "AAR"),
  matrix_to_df(distgenDISS_ALD, "Alderdice"),
  matrix_to_df(distgenDISS_BOU, "Bouma"),
  matrix_to_df(distgenDISS_BRI, "Bright"),
  matrix_to_df(distgenDISS_DIA, "Diaphus"),
  matrix_to_df(distgenDISS_EDI, "Edisto"),
  matrix_to_df(distgenDISS_EFGB, "EFGB"),
  matrix_to_df(distgenDISS_ELV, "Elvers"),
  matrix_to_df(distgenDISS_FTR, "FTR"),
  matrix_to_df(distgenDISS_GEY, "Geyer"),
  matrix_to_df(distgenDISS_MCG, "McGrail"),
  matrix_to_df(distgenDISS_RTR, "RTR"),
  matrix_to_df(distgenDISS_WFGB, "WFGB")
)

# Create the faceted plot
ggplot(all_distances, aes(x = distance)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  facet_wrap(~ site, scales = "free_y", ncol = 3) +
  labs(title = "Genetic Distance Distributions by Site",
       x = "Genetic Distance",
       y = "Count") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#3B8C6E"),
    strip.text = element_text(color = "white", face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16)
  )

# Save the plot
ggsave("genetic_distance_histograms_ggplot.pdf", width = 12, height = 10)
```


Convert SNP data into a genclone object
```{r}
pc <- as.genclone(snps)
pc
```

#Perform clone filtering
Collapse multilocus genotypes (MLGs) based on genetic distance, using a threshold of 2000. 
The result reduces the 288 original genotypes to 189 multilocus genotypes
After further investigations into clone groups the clone corrected object will change
```{r}
mlg.filter(pc, threshold = 2000, missing = "asis", distance = "diss.dist")
mlg.filter(pc,  missing = "asis", distance = "diss.dist") <- 2000
pc
save(pc, file='./pc.Rdata')

#Genotype information:

#     189 contracted multilocus genotypes
#         (2000) [t], (diss.dist) [d], (farthest) [a] 
#     288 diploid individuals
#   23869 codominant loci
# 13 populations defined
```

# Look into initial identified clones
These clones will be further investigated and final clone samples will be defined through subsequent analyses
```{r}
# Get the original multilocus genotype assignments before filtering
original_mlg <- mlg.vector(snps)

# Get the current filtered multilocus genotype assignments
filtered_mlg <- mll(pc)

# Create a data frame to compare
comparison <- data.frame(
  sample_names = indNames(pc),
  original_mlg = original_mlg,
  filtered_mlg = filtered_mlg
)

# Find samples that were merged as clones
clone_groups_initial <- split(comparison$sample_names, comparison$filtered_mlg)

# Print groups where multiple samples were considered clones
for (group in clone_groups_initial) {
  if (length(group) > 1) {
    print("These samples were considered clones:")
    print(group)
  }
}

# Count original vs. filtered MLGs
print(paste("Original number of MLGs:", length(unique(original_mlg))))
print(paste("After filtering (at threshold 2000): 189 MLGs"))
```

# Look into clone group stats and summaries
```{r}
# 1. Initial setup and metadata preparation
se_metadata$sample_id <- rownames(se_metadata)

# Create initial clone IDs from your mlg.filter results
initial_clone.ids <- mlg.id(pc)

# 2. Create base clone analysis dataframe with initial clone IDs
clone_groups_initial <- data.frame(
  sample_id = rownames(se_metadata),
  se_metadata
) %>%
  left_join(
    data.frame(
      sample_id = unlist(initial_clone.ids),
      clone_group = rep(names(initial_clone.ids), sapply(initial_clone.ids, length))
    ),
    by = "sample_id"
  )

# 3. Create detailed clone summaries
clone_details_initial <- clone_groups_initial %>%
  group_by(clone_group) %>%
  filter(n() > 1) %>%
  arrange(clone_group, SITE) %>%
  summarise(
    samples = paste(sample_id, collapse=", "),
    sites = paste(unique(SITE), collapse=", "),
    depths = paste(sort(DEPTH), collapse=", ")
  )

# 4. Create comprehensive statistics
clone_stats_initial <- list(
  # Group level statistics
  group_summary = clone_groups_initial %>%
    group_by(clone_group) %>%
    summarize(
      n_samples = n(),
      n_sites = n_distinct(SITE),
      sites = paste(sort(unique(SITE)), collapse=", "),
      depth_range = paste(min(DEPTH), "-", max(DEPTH), "m"),
      sample_ids = paste(sample_id, collapse=", "),
      depths = paste(sort(DEPTH), collapse=", ")
    ) %>%
    filter(n_samples > 1),
  
  # Overall statistics
  overall_stats = list(
    total_groups = length(unique(clone_groups_initial$clone_group)),
    cloned_samples = sum(table(clone_groups_initial$clone_group) > 1),
    size_distribution = table(table(clone_groups_initial$clone_group)),
    original_samples = nInd(pc),
    corrected_samples = length(unique(clone_groups_initial$clone_group)),
    samples_per_pop = table(pop(pc))
  )
)

# 5. Save all results with "initial" in the filename
save(clone_groups_initial, clone_details_initial, clone_stats_initial, 
     file='initial_clone_analysis_results.Rdata')

# 6. Print summary information
print("Initial Clone Group Details:")
print(clone_details_initial, n=Inf)
print("\nInitial Clone Statistics:")
with(clone_stats_initial$overall_stats, {
  cat("Total initial clone groups:", total_groups, "\n")
  cat("Number of cloned samples in initial analysis:", cloned_samples, "\n")
  cat("Initial clone group size distribution:\n")
  print(size_distribution)
  cat("\nSamples before correction:", original_samples, "\n")
  cat("Unique multilocus genotypes after initial correction:", corrected_samples, "\n")
  cat("\nSamples per population in initial dataset:\n")
  print(samples_per_pop)
})
```



# Investigating clone groups with multiple sites
Clone group 182 included sites Geyer and EFGB
Clone group 221 included sites Geyer and EFGB
Clone group 195 included sites Elvers and Bouma
Examine whether samples that were clustered as clones but appear at different sites are truly genetically identical or if they might have been erroneously grouped due to the distance threshold (was originally 5000)
```{r}
# Create target_groups object for specific clone groups of interest (182, 2211, 195)
target_groups <- clone_groups_initial %>%
  filter(clone_group %in% c(182, 221, 195))

# Look at the structure of the distance matrix
dim(as.matrix(distgenDISS))
rownames(as.matrix(distgenDISS))[1:10]  # Check first few sample names

# Define matrix_samples as the rownames of the distance matrix
matrix_samples <- rownames(as.matrix(distgenDISS))

# Filter for only samples that exist in distance matrix
available_samples <- target_groups %>%
  filter(sample_id %in% matrix_samples)

# Available samples by clone group:
print(available_samples %>% 
      select(clone_group, sample_id, SITE, DEPTH) %>%
      arrange(clone_group))

# Get distances just for these samples
clone_distances <- as.matrix(distgenDISS)[available_samples$sample_id, available_samples$sample_id]
print("Genetic distances between available samples:")
print(clone_distances)

# Get context by comparing to overall distance distribution
dist_summary <- summary(as.vector(distgenDISS))
print("Overall distance distribution summary:")
print(dist_summary)
```

#Looking into pairs SH10553/SH10653 (distance=951) and SH10781/SH10817 (distance=1475)
```{r}
# Extract SNP data for these specific samples
pairs_to_check <- c("SH10553", "SH10653", "SH10781", "SH10817")

# Get the SNP data - first let's look at the format of our data
str(snps)

# Extract genotypes for our pairs
snp_data <- tab(snps)[pairs_to_check,]

# Compare SNPs between first pair
pair1_diff <- snp_data["SH10553",] != snp_data["SH10653",]
pair2_diff <- snp_data["SH10781",] != snp_data["SH10817",]

# Summary of differences
cat("Number of SNP differences between SH10553/SH10653:", sum(pair1_diff), "\n")
cat("Number of SNP differences between SH10781/SH10817:", sum(pair2_diff), "\n")

# Look at the actual differing positions
cat("First few SNP differences between SH10553/SH10653:\n")
head(which(pair1_diff))
print(snp_data[c("SH10553","SH10653"), which(pair1_diff)[1:10]])

cat("First few SNP differences between SH10781/SH10817:\n")
head(which(pair2_diff))
print(snp_data[c("SH10781","SH10817"), which(pair2_diff)[1:10]])
```

# Count SNP differences between the above pairs
```{r}
# Extract SNP data for these specific samples
pairs_to_check <- c("SH10553", "SH10653", "SH10781", "SH10817")
snp_data <- tab(snps)[pairs_to_check,]

# Function to count meaningful genotype differences
count_genotype_differences <- function(sample1, sample2) {
  # Convert the data to a format that makes comparison easier
  geno1 <- snp_data[sample1,]
  geno2 <- snp_data[sample2,]
  
  # Count differences - considering the actual genotype states
  differences <- sum(geno1 != geno2, na.rm=TRUE)
  
  # Get the locations of differences for examination
  diff_positions <- which(geno1 != geno2)
  
  return(list(
    total_differences = differences,
    diff_positions = diff_positions,
    sample_differences = cbind(
      geno1[diff_positions[1:min(10, length(diff_positions))]],
      geno2[diff_positions[1:min(10, length(diff_positions))]]
    )
  ))
}

# Compare both pairs
pair1_diffs <- count_genotype_differences("SH10553", "SH10653")
pair2_diffs <- count_genotype_differences("SH10781", "SH10817")

# Print results
cat("Differences between SH10553/SH10653:", pair1_diffs$total_differences, "\n")
cat("First few genotype differences:\n")
print(pair1_diffs$sample_differences)

cat("Differences between SH10781/SH10817:", pair2_diffs$total_differences, "\n")
cat("First few genotype differences:\n")
print(pair2_diffs$sample_differences)
```

# Checking SH10781 (Elvers) against all Bouma clones
```{r}
# Get all samples from clone group containing SH10781
clone_group_samples <- clone_analysis %>%
  filter(clone_group == clone_analysis$clone_group[clone_analysis$sample_id == "SH10781"]) %>%
  pull(sample_id)

print("All samples in this clone group:")
print(clone_group_samples)

# Get the data from original snps object
snp_data_full <- tab(snps)[clone_group_samples,]

# Modify our comparison function for the full dataset
count_genotype_differences <- function(sample1, sample2) {
  geno1 <- snp_data_full[sample1,]
  geno2 <- snp_data_full[sample2,]
  
  differences <- sum(geno1 != geno2, na.rm=TRUE)
  diff_positions <- which(geno1 != geno2)
  
  return(list(
    total_differences = differences,
    diff_positions = diff_positions,
    sample_differences = cbind(
      geno1[diff_positions[1:min(10, length(diff_positions))]],
      geno2[diff_positions[1:min(10, length(diff_positions))]]
    )
  ))
}

# Compare SH10781 to all other samples in the group
for(sample in clone_group_samples[clone_group_samples != "SH10781"]) {
  diffs <- count_genotype_differences("SH10781", sample)
  
  cat("\nDifferences between SH10781 and", sample, ":", diffs$total_differences, "\n")
  cat("First few genotype differences:\n")
  print(diffs$sample_differences[1:5,])
}
```

# Create final clone-corrected snps dataset
Includes individuals confirmed not to be true clones between sites
Retains only one individual per clone group to avoid bias from identical genotype
Removing Clone group 182 and 221 because their sample numbers are similar and likely got switched during sample processing/extraction. Removing Bouma sample from Elvers clone group 195 because it is a separate individual than the Elvers individuals in this group
```{r}
#library(poppr)
# Start with original data
pc <- as.genclone(snps)

# Set distance threshold for initial filtering
mlg.filter(pc, distance = "diss.dist") <- 2000

# Get current MLG assignments
current_mlgs <- mlg.vector(pc)
custom_mlgs <- current_mlgs

# Get the specific Elvers sample from group 195
elvers_195 <- clone_analysis %>%
  filter(clone_group == "195", SITE == "Elvers") %>%
  pull(sample_id)

# Get all samples from groups 221 and 182
samples_221_182 <- clone_analysis %>%
  filter(clone_group %in% c("221", "182")) %>%
  pull(sample_id)

# Combine all samples to separate
samples_to_separate <- c(elvers_195, samples_221_182)

# Give these samples unique MLGs
for(i in seq_along(samples_to_separate)) {
  custom_mlgs[which(indNames(pc) == samples_to_separate[i])] <- max(as.numeric(current_mlgs)) + i
}

# Set strata and perform clone correction
strata(pc) <- data.frame(pop = pop(pc))
pc.corr <- clonecorrect(pc, strata = ~pop)

# Convert to genind object
snps_mlg <- genclone2genind(pc.corr)

# Save the results
save(snps_mlg, file='./snps_mlg.Rdata')

# Verify the reduction in samples
print("Number of samples before correction:")
print(nInd(pc))
print("Number of MLGs before correction:")
print(length(unique(mlg.vector(pc))))
print("Number of samples after correction:")
print(nInd(pc.corr))
print("Number of MLGs after correction:")
print(length(unique(mlg.vector(pc.corr))))

# Verify the conversion
print(snps_mlg)
print("Samples per population:")
print(summary(snps_mlg$pop))
#NEW CLONE CORRECTED DATA
```

#Attempting to simplify clone analysis
This code performs a comprehensive analysis of clonal structure across multiple sampling sites:
 1. Prepares helper functions for genetic distance analysis and network visualization
 2. Corrects metadata errors in specific sample coordinates
 3. Performs site-by-site analysis to identify potential clones based on genetic similarity
 4. Creates a consolidated dataframe linking samples to their clone groups
 5. Generates detailed summaries of clone group characteristics
 6. Creates a filtered dataset containing only clonal samples for downstream analysis
 7. Cleans problematic cross-site clone groups
 8. Analyzes depth distribution patterns of clones within each site
 9. Recalculates spatial distribution of clones using corrected coordinates
 10. Saves all results for downstream analyses

Outputs:
 - Network visualizations of genetic similarity within each site
 - Identification of clone groups and their member samples
 - Analysis of clone distribution by depth and geographic distance
 - Clean datasets for further population genetic analyses
Required existing objects (from your environment):
 snps: Master genind object containing SNP data
 se_metadata: Sample metadata with coordinates, depth, and site information
 clone_details: Existing details about clone groups
 all_sites: Vector of site names to analyze
 spatial_analysis: Existing spatial analysis object referenced in prints
 
 Object Outputs:
 all_results: List containing clone analysis results for each site
 final_clone_details: Combined dataframe of clone groups across all sites
 site_summary: Table showing count of clone groups per site
 clone_analysis_updated: Complete dataframe with sample metadata and assigned clone groups
 comparison: Simple dataframe comparing original and updated clone analyses
 clone_group_summary: Detailed statistics for each clone group
 clone_samples: Vector of sample IDs that belong to clone groups
 snps_clones: Genind object containing only samples that are part of clone groups
 clone_metadata: Subset of metadata for clonal samples only
 bouma_samples_182_221: Filtered dataset with just Bouma samples from problematic groups
 clean_clone_details: Cleaned version of clone details after removing problematic groups
 depth_analysis: Statistics on depth distribution of clones per site and group
 clone_spatial: Spatial data for clone groups with corrected coordinates
 all_clone_distances: Calculated geographic distances between clone samples
 clone_site_summary: Site-level summary of clone spatial distribution
```{r}
#===== IMPROVED CLONE ANALYSIS WORKFLOW =====

# 1. Setup and functions
#-----------------------------------------
# Helper function to prepare site data and visualization parameters
prepare_site_data <- function(site_name, threshold = 2000) {
  # Get site data
  site_snps <- snps[pop=site_name]
  n_samples <- length(indNames(site_snps))
  
  if(n_samples < 2) {
    return(list(valid = FALSE, message = paste("Skipping", site_name, "- not enough samples")))
  }
  
  # Calculate distances
  dist_mat <- as.matrix(diss.dist(site_snps, percent = FALSE, mat = FALSE))
  
  # Setup visualization parameters
  edge_colors <- matrix("grey80", n_samples, n_samples)
  edge_widths <- matrix(0.5, n_samples, n_samples)
  
  # Find close pairs
  close_pairs <- which(dist_mat < threshold & dist_mat > 0, arr.ind = TRUE)
  
  # Set colors and widths for close pairs if they exist
  if(length(close_pairs) > 0) {
    for(i in 1:nrow(close_pairs)) {
      edge_colors[close_pairs[i,1], close_pairs[i,2]] <- "darkgreen"
      edge_colors[close_pairs[i,2], close_pairs[i,1]] <- "darkgreen"
      edge_widths[close_pairs[i,1], close_pairs[i,2]] <- 2
      edge_widths[close_pairs[i,2], close_pairs[i,1]] <- 2
    }
  }
  
  return(list(
    valid = TRUE,
    site_snps = site_snps,
    dist_mat = dist_mat,
    n_samples = n_samples,
    close_pairs = close_pairs,
    edge_colors = edge_colors,
    edge_widths = edge_widths
  ))
}

# Function to create network plots
create_network_plot <- function(site_data, site_name, threshold = 2000, save_to_file = FALSE) {
  if(save_to_file) {
    pdf(file = paste0("clone_network_", site_name, ".pdf"))
  }
  
  qgraph(1/site_data$dist_mat, 
         layout='spring', 
         vsize=5, 
         minimum = 1/threshold,
         edge.color = site_data$edge_colors,
         edge.width = site_data$edge_widths,
         labels = rownames(site_data$dist_mat),
         title = paste("Clone Network -", site_name))
  
  if(save_to_file) {
    dev.off()
  }
}

# Main function to analyze site clones
analyze_site_clones <- function(site_name, threshold = 2000, save_plot = FALSE) {
  tryCatch({
    # Prepare site data
    site_data <- prepare_site_data(site_name, threshold)
    
    if(!site_data$valid) {
      cat("\n", site_data$message, "\n")
      return(NULL)
    }
    
    # Create network visualization
    cat("\nCreating network for", site_name, "\n")
    create_network_plot(site_data, site_name, threshold, save_plot)
    
    # Initialize clone groups
    clone_groups <- list()
    processed <- c()
    
    # Process close pairs if they exist
    if(length(site_data$close_pairs) > 0 && nrow(site_data$close_pairs) > 0) {
      for(i in 1:nrow(site_data$close_pairs)) {
        s1 <- rownames(site_data$dist_mat)[site_data$close_pairs[i,1]]
        s2 <- colnames(site_data$dist_mat)[site_data$close_pairs[i,2]]
        
        if(all(c(s1, s2) %in% processed)) next
        
        group_found <- FALSE
        for(j in seq_along(clone_groups)) {
          if(any(c(s1, s2) %in% clone_groups[[j]])) {
            clone_groups[[j]] <- unique(c(clone_groups[[j]], s1, s2))
            group_found <- TRUE
            break
          }
        }
        
        if(!group_found) {
          clone_groups[[length(clone_groups) + 1]] <- c(s1, s2)
        }
        
        processed <- unique(c(processed, s1, s2))
      }
    }
    
    # Create results dataframe
    if(length(clone_groups) > 0) {
      results <- data.frame(
        site = site_name,
        clone_group = paste(site_name, seq_along(clone_groups), sep="_"),
        samples = sapply(clone_groups, paste, collapse=", "),
        n_samples = sapply(clone_groups, length),
        stringsAsFactors = FALSE
      )
    } else {
      results <- data.frame(
        site = site_name,
        clone_group = character(0),
        samples = character(0),
        n_samples = integer(0),
        stringsAsFactors = FALSE
      )
    }
    
    # Print summary
    cat("\nResults for", site_name, ":\n")
    cat("Total samples:", nrow(site_data$dist_mat), "\n")
    cat("Number of clone groups:", length(clone_groups), "\n")
    
    if(length(clone_groups) > 0) {
      cat("Samples in clone groups:", length(unique(unlist(clone_groups))), "\n")
      cat("\nClone groups:\n")
      print(results)
      
      if(length(site_data$close_pairs) > 0) {
        cat("\nDistances between clone pairs:\n")
        print(summary(site_data$dist_mat[site_data$close_pairs]))
      }
    }
    
    return(results)
    
  }, error = function(e) {
    cat("\nError processing", site_name, ":", conditionMessage(e), "\n")
    return(NULL)
  })
}

# 2. Correct metadata errors
#-----------------------------------------
# Correct Diaphus sample SH11757 latitude
cat("Current value for SH11757 latitude:", se_metadata["SH11757", "LAT"], "\n")
se_metadata["SH11757", "LAT"] <- 28.08812
cat("New value for SH11757 latitude:", se_metadata["SH11757", "LAT"], "\n")

# 3. Initial site-by-site clone analysis
#-----------------------------------------
# Extract unique site names from the metadata
all_sites <- unique(se_metadata$SITE)
# Process all sites
all_results <- list()
cat("\nProcessing", length(all_sites), "sites...\n")

for(site in all_sites) {
  cat("\n=====================================")
  cat("\nAnalyzing", site)
  cat("\n=====================================\n")
  result <- analyze_site_clones(site)
  if(!is.null(result) && nrow(result) > 0) {
    all_results[[site]] <- result
  }
}

# Combine results
if(length(all_results) > 0) {
  final_clone_details <- do.call(rbind, all_results)
  
  # Print summary
  cat("\nFinal summary across all sites:\n")
  site_summary <- table(final_clone_details$site)
  print(site_summary)
  
  cat("\nTotal clone groups found:", nrow(final_clone_details), "\n")
  
  # Save results
  write.csv(final_clone_details, "all_sites_clone_analysis.csv", row.names = FALSE)
} else {
  cat("\nNo clone groups found in any site.\n")
}

# 4. Create consolidated clone analysis dataframe
#-----------------------------------------
# Create updated clone_analysis with correct coordinates
clone_analysis_updated <- data.frame(
  sample_id = rownames(se_metadata),
  se_metadata
) %>%
  left_join(
    data.frame(
      sample_id = unlist(strsplit(final_clone_details$samples, ", ")),
      clone_group = rep(final_clone_details$clone_group, 
                       sapply(strsplit(final_clone_details$samples, ", "), length))
    ),
    by = "sample_id"
  )

# Verify the update
cat("Dimensions of updated clone_analysis:", dim(clone_analysis_updated), "\n")
cat("Number of clone groups:", length(unique(na.omit(clone_analysis_updated$clone_group))), "\n")

# Check specific coordinates that were corrected
cat("\nVerifying corrected coordinates (e.g., SH11757):\n")
print(clone_analysis_updated[clone_analysis_updated$sample_id == "SH11757", 
                           c("sample_id", "LONG", "LAT", "DEPTH", "clone_group")])

# Instead of comparing with a non-existent object, just show the stats of the updated object
summary_stats <- data.frame(
  Dataset = "Updated",
  Samples = nrow(clone_analysis_updated),
  Clone_Groups = length(unique(na.omit(clone_analysis_updated$clone_group)))
)
print(summary_stats)

# 5. Create clone group summary
#-----------------------------------------
# Create detailed clone group summary
clone_group_summary <- clone_analysis_updated %>%
  filter(!is.na(clone_group)) %>%
  group_by(clone_group) %>%
  summarise(
    n_samples = n(),
    sites = paste(unique(SITE), collapse=", "),
    depth_range = paste(range(DEPTH), collapse="-"),
    .groups = "drop"
  )

cat("\nClone group summary:\n")
print(clone_group_summary, n=Inf)

# 6. Create filtered clone dataset for downstream analysis
#-----------------------------------------
# Extract sample IDs that have a valid clone group (not NA)
clone_samples <- clone_analysis_updated %>%
  filter(!is.na(clone_group)) %>%
  pull(sample_id)

# Check how many clone samples we found
cat("Number of clone samples identified:", length(clone_samples), "\n")

# Subset the original genind object to include only clones
snps_clones <- snps[clone_samples]

# Check how many samples are in the result
cat("Number of samples in snps_clones:", nInd(snps_clones), "\n")

# Add clone group information to the genind object
clone_metadata <- clone_analysis_updated %>%
  filter(!is.na(clone_group)) %>%
  select(sample_id, clone_group, SITE, DEPTH, LAT, LONG)

# Add metadata to the genind object
snps_clones@other$metadata <- clone_metadata

# 7. Clean-up problematic clone groups
#-----------------------------------------
# First get the Bouma samples from groups 182 and 221
bouma_samples_182_221 <- clone_details_initial %>%
  filter(clone_group %in% c("182", "221")) %>%
  separate_rows(samples, sep = ", ") %>%
  separate_rows(sites, sep = ", ") %>%
  filter(sites == "Bouma") %>%
  group_by(clone_group) %>%
  summarise(
    samples = paste(samples, collapse = ", "),
    sites = "Bouma",
    depths = paste(sort(unique(depths)), collapse = ", ")
  )

# Now create our clean clone dataset
clean_clone_details <- clone_details_initial %>%
  filter(!clone_group %in% c("182", "221")) %>%  # Remove these groups initially
  filter(sites != "Bouma, Elvers") %>%  # Remove cross-site group
  bind_rows(bouma_samples_182_221)  # Add back just the Bouma samples

# 8. Analyze depth distribution of clones
#-----------------------------------------
# Create depth analysis
depth_analysis <- clean_clone_details %>%
  separate_rows(depths, sep = ", ") %>%
  mutate(depth = as.numeric(depths)) %>%
  group_by(clone_group, sites) %>%
  summarise(
    mean_depth = mean(depth),
    depth_range = max(depth) - min(depth),
    n_samples = n()
  )

# Create visualization
ggplot(depth_analysis, aes(x=mean_depth)) +
  geom_histogram(bins=20) +
  facet_wrap(~sites) +
  labs(title="Clone Group Depth Distribution by Site",
       x="Mean Depth (m)",
       y="Number of Clone Groups") +
  theme_minimal()

print("Clone distribution by site:")
print(spatial_analysis)

print("\nDepth distribution of clone groups:")
print(depth_analysis)

# 9. Re-analyze spatial distribution with corrected coordinates
#-----------------------------------------
# Remake clone_spatial with corrected metadata
clone_spatial <- clean_clone_details %>%
 separate_rows(samples, sep = ", ") %>%
 left_join(se_metadata, by = c("samples" = "sample_id")) %>%
 group_by(clone_group, sites)

# Recalculate distances for all clone groups
all_clone_distances <- clone_spatial %>%
 group_by(clone_group, sites) %>%
 summarise(
   n_samples = n(),
   mean_distance = if(n() > 1) {
     coords <- unique(cbind(LONG, LAT))
     if(nrow(coords) > 1) {
       dist_matrix <- distm(coords, fun = distHaversine)
       mean(dist_matrix[upper.tri(dist_matrix)])/1000  # Convert to km
     } else {
       0
     }
   } else {
     0
   },
   .groups = "drop"
 ) %>%
 # Convert to meters for easier interpretation
 mutate(mean_distance_m = mean_distance * 1000)

# Create site summary
clone_site_summary <- all_clone_distances %>%
 group_by(sites) %>%
 summarise(
   n_clone_groups = n(),
   mean_clone_distance_m = mean(mean_distance_m),
   max_clone_distance_m = max(mean_distance_m),
   total_clones = sum(n_samples)
 )

print("Updated clone distribution summary by site:")
print(clone_site_summary, n=Inf)

# 10. Save final results
#-----------------------------------------
# Save all important dataframes
save(final_clone_details, all_results, site_summary, 
     clone_analysis_updated, clean_clone_details, 
     snps_clones, file = "clone_analysis_final.RData")

# Save datasets
save(snps_clones, 
     clone_analysis_updated, # sample_id, READS, SITE, LAT, LONG, DEPTH, sample_id.1, clone_group
     clean_clone_details, # clone_group, samples, sites, depths
     clone_spatial, #clone_group, samples, sites, depths, READS, SITE, LAT, LONG, DEPTH
     all_clone_distances, # clone_group, sites, n_samples, mean_distance, mean_distance_m
     file = "clean_clone_analysis_final.RData")
```

# Scatterplot of mean distances within clone groups, with points colored by depth
```{r}
# First, prepare the data
library(tidyverse)
library(geosphere)  # For distm function
library(viridis)    # For viridis color scales

get_clone_distances <- function(clone_group, metadata) {
  samples <- unlist(strsplit(clone_group, ", "))
  
  coords <- metadata %>%
    filter(sample_id %in% samples) %>%
    select(LONG, LAT)
  
  if(nrow(coords) > 1) {
    dist_matrix <- distm(as.matrix(coords), fun = distHaversine)  # Using distm with Haversine
    distances <- dist_matrix[upper.tri(dist_matrix)]
    if(length(distances) > 0) {
      return(list(
        mean_distance = mean(distances),
        max_distance = max(distances),
        sample_depths = metadata$DEPTH[metadata$sample_id %in% samples]  # Changed from 'depths' to 'sample_depths'
      ))
    }
  }
  return(list(mean_distance = 0, max_distance = 0, sample_depths = metadata$DEPTH[metadata$sample_id %in% samples]))
}

# Process the data - exclude NA clone groups and composite sites
plot_data <- clean_clone_details %>%
  # Filter out rows with NA clone_group
  filter(!is.na(clone_group)) %>%
  # Also filter out any row where sites contains a comma
  filter(!grepl(",", sites)) %>%
  rowwise() %>%
  mutate(
    sample_list = list(unlist(strsplit(samples, ", "))),
    n_samples = length(sample_list),
    # Get depth information from clone_analysis_updated
    mean_depth = mean(clone_analysis_updated$DEPTH[clone_analysis_updated$sample_id %in% sample_list], na.rm = TRUE),
    # Calculate distances
    distances = list(get_clone_distances(samples, clone_analysis_updated))
  ) %>%
  unnest_wider(distances, names_sep = "_") %>%
  ungroup()

# Plot scatter plot
# Use the correct column names after unnest_wider
ggplot(plot_data, aes(x = reorder(sites, distances_mean_distance), y = distances_mean_distance, color = mean_depth)) + 
    geom_point(size = 5, alpha = 0.8) +
    scale_color_viridis(option = "mako", direction = -1, 
                      name = "Mean Depth (m)") +
    
    # Add the mean line across the entire plot
    geom_hline(yintercept = mean(plot_data$distances_mean_distance, na.rm = TRUE),
              color = "red", linetype = "dashed", size = 1.0) +
    # Add text label for the mean
    annotate("text", 
             x = 1, # Position at first site
             y = mean(plot_data$distances_mean_distance, na.rm = TRUE) + 1.5,
             label = paste("Mean =", round(mean(plot_data$distances_mean_distance, na.rm = TRUE), 2), "m"),
             color = "red", hjust = 0) +
    scale_y_continuous(
         limits = c(0, 30),  # Set limits to 0-30 meters
         breaks = seq(0, 30, 5),  # Breaks every 5 meters
         labels = function(x) format(x, digits = 1)
     ) +
     coord_flip() +  # Flips x and y axes
     labs(title = "Mean Clone Distances by Site",
          subtitle = "Points colored by mean depth of clone group",
          x = "",  
          y = "Mean Distance Between Clones (meters)") +
     theme_light() +
     theme(
         axis.text.y = element_text(size = 13, color = "black"),
         axis.text.x = element_text(size = 13, color = "black"),
         plot.title = element_text(size = 16, face = "bold"),
         plot.subtitle = element_text(size = 12, color = "gray40"),
         panel.grid.major.y = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = "right"
     )

# Save the plot
ggsave("clone_distances_depth.pdf", width = 9, height = 8)

# Print summary statistics
cat("\nSummary statistics by site (distances in meters):\n")
clone_spatial_summary_stats <- plot_data %>%
  group_by(sites) %>%
  summarise(
    n_clones = n(),
    mean_distance = round(mean(distances_mean_distance, na.rm = TRUE), 2),
    max_distance = round(max(distances_max_distance, na.rm = TRUE), 2),
    mean_depth = round(mean(mean_depth, na.rm = TRUE), 1)
  )
print(clone_spatial_summary_stats)
```

# Plot simple histogram of mean distances within clone groups
```{r}
pdf("clone_distances_histo.pdf", width = 9, height = 8) 
hist(plot_data$distances_mean_distance, 
     breaks = 10,  
     col = "gray",  
     border = "black",  
     xlab = "Mean Distance Within Clone Groups (m)", 
     main = ""  # No title
)
dev.off()
```

# Export Clone group information for QGIS
```{r}
# Create table with samples as rows
sample_table <- clone_analysis_updated %>%
  # Start with sample info
  select(sample_id, SITE, LONG, LAT, DEPTH) %>%
  # Join with clone group information
  left_join(
    final_clone_details %>%
      mutate(sample_id = strsplit(samples, ", ")) %>%
      unnest(sample_id) %>%
      select(clone_group, sample_id),
    by = "sample_id"
  ) %>%
  # Arrange by site and clone group
  arrange(SITE, clone_group, sample_id)

# Write as tab-delimited text file
write.table(sample_table, 
            "samples_with_clones.txt", 
            sep = "\t",            
            row.names = FALSE,     
            quote = FALSE)         

# Preview the table
head(sample_table)
```

# Table all singletons included in snps genid object before filtering
This can be used for mapping all sample points 
```{r}
# Print initial counts
cat("Initial samples in snps genind:", nInd(snps), "\n")
cat("Sample names in snps:", length(indNames(snps)), "\n")

# Create table with samples that are in snps
sample_table <- clone_analysis_updated %>%
  # Filter for samples that are in snps
  filter(sample_id %in% indNames(snps)) %>%
  # Start with sample info
  select(sample_id, SITE, LONG, LAT, DEPTH) %>%
  # Join with clone group information
  left_join(
    final_clone_details %>%
      mutate(sample_id = strsplit(samples, ", ")) %>%
      unnest(sample_id) %>%
      select(clone_group, sample_id),
    by = "sample_id"
  ) %>%
  # Arrange by site and clone group
  arrange(SITE, clone_group, sample_id)

# Print sample counts at different stages
cat("\nSamples in clone_analysis_updated:", nrow(clone_analysis_updated), "\n")
cat("Samples after filtering for snps:", nrow(sample_table), "\n")
cat("Samples with clone groups:", sum(!is.na(sample_table$clone_group)), "\n")

# Write as tab-delimited text file
write.table(sample_table, 
            "Swiftia_popgen_all_samples_postfilter.txt", 
            sep = "\t",            
            row.names = FALSE,     
            quote = FALSE)         

# Preview the table
head(sample_table)
dim(sample_table)
```

# Create table of just clone groups
With clone group categories 1-8 for mapping
```{r}
# Create table with numeric clone ID
clone_sample_table <- clone_analysis_updated %>%
  select(sample_id, SITE, LONG, LAT, DEPTH) %>%
  # Join with clone group information
  left_join(
    final_clone_details %>%
      mutate(sample_id = strsplit(samples, ", ")) %>%
      unnest(sample_id) %>%
      select(clone_group, sample_id),
    by = "sample_id"
  ) %>%
  # Remove rows where clone_group is NA
  filter(!is.na(clone_group)) %>%
  # Add new cloneID column by extracting number after underscore
  mutate(cloneID = sub(".*_", "", clone_group)) %>%
  # Arrange by site and clone group
  arrange(SITE, clone_group, sample_id)

# Write as tab-delimited text file
write.table(clone_sample_table, 
            "clone_samples_with_ID.txt", 
            sep = "\t",            
            row.names = FALSE,     
            quote = FALSE)         

# Preview the table
head(clone_sample_table)
```


# Checking Edisto distances to make sure we aren't missing any clones
```{r}
# Get Edisto samples
edisto_snps <- snps[pop='Edisto']

# Calculate distances
distgenDISS_ED <- diss.dist(edisto_snps, percent = FALSE, mat = FALSE)

# Create histogram
hist(distgenDISS_ED, 
     breaks = 30,
     main = "Distribution of Genetic Distances - Edisto",
     xlab = "Genetic Distance",
     ylab = "Frequency")

# Add vertical line at threshold
abline(v = 2000, col = "red", lwd = 2)
text(2000, max(hist(distgenDISS_ED, breaks = 30, plot = FALSE)$counts), 
     "Threshold (2000)", pos = 4, col = "red")

# Print summary statistics
cat("\nSummary statistics for Edisto distances:\n")
print(summary(distgenDISS_ED))

# Get number of pairs below threshold
n_below_threshold <- sum(distgenDISS_ED < 2000)
cat("\nNumber of pairs below threshold (2000):", n_below_threshold, "\n")

# Show closest pairs
if(n_below_threshold > 0) {
  close_pairs <- which(as.matrix(distgenDISS_ED) < 2000 & as.matrix(distgenDISS_ED) > 0, arr.ind = TRUE)
  if(nrow(close_pairs) > 0) {
    cat("\nPairs below threshold:\n")
    dist_matrix <- as.matrix(distgenDISS_ED)
    close_pairs_df <- data.frame(
      Sample1 = rownames(dist_matrix)[close_pairs[,1]],
      Sample2 = colnames(dist_matrix)[close_pairs[,2]],
      Distance = dist_matrix[close_pairs]
    )
    print(close_pairs_df[order(close_pairs_df$Distance),])
  }
}
```

#Euclidean distance b/w allele frequencies
```{r}
distgenEUCL <- dist(snps_mlg, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
hist(distgenEUCL)
save(distgenEUCL, file = './distgenEUCL_se.Rdata')
```

#Pairwise number of allelic differences between two individuals (diss.dist {poppr})
```{r}
distgenDISS <- diss.dist(snps_mlg, percent = FALSE, mat = FALSE)
hist(distgenDISS)
save(distgenDISS, file = './distgenDISS_se.Rdata')
#Visualize data via histograms and a network graph
qgraph(1/distgenDISS, layout='spring', vsize=5)
```

================================================================================
------------------------Further Data Filtering:He- HWequil.---------------------
================================================================================

#Convert genetic data to a hierfstat format and calcuate basic stats
```{r}
# Check the distribution of Het
snps_mlg_hfs <- genind2hierfstat(snps_mlg) # Create hierfstat object

# Perform basic statistics
basicstat_all <- basic.stats(snps_mlg_hfs, diploid = TRUE, digits = 2)
save(basicstat_all, file='basicstat_all.Rdata')

# Extract heterozygosity data
Het_all <- as.data.frame(basicstat_all$perloc)

# Create histogram and check the distribution of observed heterozygosity (Ho)
ggplot(Het_all, aes(x=Ho))+
  geom_histogram(binwidth = 0.01,color="black", fill="white")+
  theme_classic()

# Check the markers with high heterozygosity Ho > 0.5
Het_tab_all <- Het_all[which(Het_all$Ho  > 0.5), ]
dim(Het_tab_all)
```

#Filter out markers with high Ho
```{r}
#If any are found they can be removed by creating a vector of loci names and using adegenet to filter them out of dataset
toRemove_all <- which(Het_all$Ho  > 0.5)
snps_mlg_he = snps_mlg[loc=-toRemove_all]
snps_mlg_he
save(snps_mlg_he, file='snps_mlg_he.Rdata')
```

# Hardy-Weinberg Equilibrium (HWE)
Package: Pegas
Test for Hardy-Weinberg Equilibrium and remove loci that deviate significantly from HWE (p<0.01)
identifies loci that deviate from neutral expectations due to various factors including selection, but also non-selective forces like assortive mating, population structure, or genotyping errors
```{r}
#library(pegas)
snps_mlg_he.hw <- hw.test(snps_mlg_he, B = 1000)
save(snps_mlg_he.hw, file='snps_mlg_he.hw.Rdata')
# View histogram of p-values
hist(snps_mlg_he.hw[,4], breaks=50, 
     main="Distribution of Hardy-Weinberg P-values",
     xlab="P-value", ylab="Frequency")
```

#Identify and remove loci not in HWE
```{r}
loci_toRemove_hw <- which(snps_mlg_he.hw[,4] < 0.01)
length(loci_toRemove_hw)

# Create new neutral dataset
snps_mlg.neutral <- snps_mlg_he[loc=-loci_toRemove_hw]
snps_mlg.neutral
```

# Verify new HWE neutral dataset
```{r}
# Check object types
class(snps_mlg_he)
class(snps_mlg.neutral)

# Verify genind object integrity
str(snps_mlg_he)

#Object's structure should confirm:
#Individuals: 192
#Loci: 23,605
#Populations: 13 (unique pop levels)
#Data Type:
#Codominant genetic data
#Diploid (ploidy = 2)

# Summarize p-value distribution
summary(snps_mlg_he.hw[,4])
```

===================================================
-----------------------BayeScan--------------------
===================================================
This section needs more refinement before publication
# Prepare data for BayeScan
```{r}
# More efficient BayeScan conversion
#library(adegenet)
#library(dplyr)

genind2bayescan_efficient <- function(genind_obj, outfile) {
  # Get populations and loci info
  pops <- levels(genind_obj$pop)
  n_pops <- length(pops)
  n_loci <- nLoc(genind_obj)
  
  # Open file connection
  con <- file(outfile, "w")
  
  # Write header
  writeLines(paste("[loci]=", n_loci), con)
  writeLines(paste("[populations]=", n_pops), con)
  
  # Efficient processing
  for(loc in 1:n_loci) {
    # Get locus name
    loc_name <- locNames(genind_obj)[loc]
    
    # Extract allele data for this locus
    loc_data <- genind_obj@tab[, grep(paste0("^", loc_name), colnames(genind_obj@tab))]
    
    # Count total genes and allele frequencies
    total_genes <- nrow(loc_data) * 2
    allele_counts <- colSums(loc_data, na.rm = TRUE)
    
    # Ensure we have two allele counts
    if(length(allele_counts) >= 2) {
      # Write locus information
      writeLines(paste(loc, total_genes, "2", 
                       allele_counts[1], allele_counts[2]), 
                 con)
    }
  }
  
  # Close file connection
  close(con)
}

# Use the function
genind2bayescan_efficient(snps_mlg.neutral, "swiftia_snps_bayescan_neutral.txt")

# Verify file creation
file.info("swiftia_snps_bayescan_neutral.txt")
```

# Verify BayeScan file
```{r}
# Read first few lines
readLines("swiftia_snps_bayescan_neutral.txt", n=10)

# Count total lines
total_lines <- length(readLines("swiftia_snps_bayescan_neutral.txt"))

# Verify header and first few data lines
validate_bayescan_file <- function(filename) {
  lines <- readLines(filename)
  
  # Check header
  header_loci <- lines[1]
  header_pops <- lines[2]
  
  cat("Loci header:", header_loci, "\n")
  cat("Populations header:", header_pops, "\n")
  
  # Extract number of loci and populations
  n_loci <- as.numeric(gsub("\\[loci\\]=", "", header_loci))
  n_pops <- as.numeric(gsub("\\[populations\\]=", "", header_pops))
  
  cat("Number of loci:", n_loci, "\n")
  cat("Number of populations:", n_pops, "\n")
  cat("Total lines in file:", total_lines, "\n")
  
  # Check first few data lines
  data_lines <- lines[3:min(10, length(lines))]
  cat("\nFirst few data lines:\n")
  print(data_lines)
}

# Run validation
validate_bayescan_file("swiftia_snps_bayescan_neutral.txt")
```

#Further checks
```{r}
# Check for negative allele counts
if (any(gen_data < 0, na.rm = TRUE)) {
  stop("Error: Negative allele frequencies detected.")
}

# Check for missing values
if (any(is.na(gen_data))) {
  warning("Warning: Missing values detected in dataset.")
}

# Check for allele count consistency
allele_sums <- apply(snps_mlg.neutral@tab, 2, sum, na.rm = TRUE)
if (any(allele_sums > (2 * nrow(snps_mlg.neutral@tab)))) {
  warning("Warning: Some loci exceed maximum diploid allele count.")
}
```

# Generate population assignment file in R
```{r}
# Create a numeric population assignment
pops <- levels(snps_mlg.neutral$pop)
pop_numeric <- as.numeric(factor(pop(snps_mlg.neutral), levels=pops))

# Write to file
write(pop_numeric, file = "swiftia_population_assignments.txt", ncolumns = 1)

# Verify contents
cat("First few lines of population assignments:\n")
head(readLines("swiftia_population_assignments.txt"))

# Verify total length and population distribution
length(pop_numeric)  # Should be 192
table(pop_numeric, pop(snps_mlg.neutral))
```

# If above population assignment does not work from R, use bash to create population assignment file:
This was ultimately what I ended up using
```{bash}
# Populations and their sizes
# Edisto: 16, AAR: 12, RTR: 26, Alderdice: 13, Geyer: 5, EFGB: 20, WFGB: 5, 
# Bouma: 13, Elvers: 18, Bright: 8, Diaphus: 22, McGrail: 16, FTR: 18

awk 'BEGIN {
    populations[1] = 16   # Edisto
    populations[2] = 12   # AAR
    populations[3] = 26   # RTR
    populations[4] = 13   # Alderdice
    populations[5] = 5    # Geyer
    populations[6] = 20   # EFGB
    populations[7] = 5    # WFGB
    populations[8] = 13   # Bouma
    populations[9] = 18   # Elvers
    populations[10] = 8   # Bright
    populations[11] = 22  # Diaphus
    populations[12] = 16  # McGrail
    populations[13] = 18  # FTR
}
{
    for(pop in populations) {
        for(i=1; i<=populations[pop]; i++) {
            print pop
        }
    }
}' /dev/null > swiftia_population_assignments_precise.txt
```


# Format the swiftia_snps_bayescan_neutral.txt file further for Bayescan
Creates a population-specific with complete locus coverage
This will take a few minutes on server
```{bash}
# Log into server

# Move swiftia_snps_bayescan_neutral.txt to server 
# Navigate to correct directory


# Format file
(
echo "[loci]=16764"
echo ""
echo "[populations]=13"
echo ""
pop_sizes=(16 12 26 13 5 20 5 13 18 8 22 16 18)
total_loci=16764
current_line=1

for pop in "${!pop_sizes[@]}"; do
    pop_num=$((pop+1))
    pop_size=${pop_sizes[$pop]}
    
    echo "[pop]=$pop_num"
    
    # Create an associative array to store locus data
    declare -A locus_data
    
    # First pass: read all data for this population into memory
    while IFS=' ' read -r locus_num _ _ allele1 allele2; do
        if [[ $locus_num =~ ^[0-9]+$ ]]; then
            locus_data[$locus_num]="$allele1 $allele2"
        fi
    done < swiftia_snps_bayescan_neutral.txt
    
    # Second pass: output data for each locus
    for ((locus=1; locus<=total_loci; locus++)); do
        # Progress indicator
        if ((locus % 1000 == 0)); then
            echo "Processing population $pop_num, locus $locus / $total_loci" >&2
        fi
        
        # Use the stored data if available, otherwise use zeros
        if [[ -n "${locus_data[$locus]}" ]]; then
            read allele1 allele2 <<< "${locus_data[$locus]}"
        else
            allele1=0
            allele2=0
        fi
        
        printf "  %d  384  2  %d  %d\n" $locus $allele1 $allele2
    done
    
    echo ""
    unset locus_data
done
) > swiftia_snps_bayescan_neutral_comprehensive.txt
```

# Run Bayescan
```{bash}
#Find exact path to bayscan on server for following command
ls /opt/mabaforge/envs/popgen/bin/bayescan*

# Run using nohup, using 100 threads takes about 8 hours to run on server
nohup /opt/mambaforge/envs/popgen/bin/bayescan_2.1 swiftia_snps_bayescan_neutral_comprehensive.txt -population_prior swiftia_population_assignments_precise.txt -threads 100 -od . -o swiftia_bayescan_output

# Open additional terminal windo to server and view verify file
head -n 20 swiftia_bayescan_output_Verif.txt
```


# Analyze BayeScan results
- removes outliers and creates a neutral dataset
-The burnin (-burn) and number of iterations (-n) might need adjustment based on your dataset
-The prior odds (-pr_odds) can be adjusted based on how many outliers you expect
-Monitor convergence using the chain files produced by BayeScan
-Consider running multiple times with different seeds to ensure consistency
-The q-value threshold (0.05 in this example) can be adjusted based on how conservative you want to be

# Analyze BayeScan Outliers
```{r}
# Load required libraries
library(ggplot2)

# Read BayeScan results
results <- read.table("swiftia_bayescan_output_fst.txt", header=TRUE)
colnames(results) <- c("prob", "log10.PO", "qval", "alpha", "fst")

# Identify outliers
outliers <- which(results$qval < 0.05)  # Using q-value threshold of 0.05
print(paste("Number of outliers:", length(outliers)))

# Create more detailed plots
# Plot 1: FST vs log10(PO)
ggplot(results, aes(x = fst, y = log10.PO)) +
  geom_point(alpha = 0.6) +
  geom_point(data = results[outliers,], color = "red", size = 3) +
  theme_bw() +
  labs(x = "FST", y = "log10(PO)", 
       title = "BayeScan Results: FST vs log10(PO)")

# Plot 2: Q-value distribution
ggplot(results, aes(x = qval)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  theme_bw() +
  labs(x = "q-value", y = "Count", 
       title = "Distribution of q-values")

# Create summary table of outlier loci
if(length(outliers) > 0) {
  outlier_summary <- results[outliers,]
  outlier_summary$locus <- outliers
  write.csv(outlier_summary, "outlier_summary.csv", row.names = FALSE)
}

# Calculate summary statistics
summary_stats <- data.frame(
  mean_fst = mean(results$fst),
  sd_fst = sd(results$fst),
  median_fst = median(results$fst),
  mean_alpha = mean(results$alpha),
  sd_alpha = sd(results$alpha),
  num_outliers = length(outliers),
  total_loci = nrow(results)
)
write.csv(summary_stats, "bayescan_summary_stats.csv", row.names = FALSE)
#   mean_fst     sd_fst median_fst mean_alpha  sd_alpha num_outliers total_loci
#1 0.1122051 0.02377327    0.11551 -0.0872178 0.3953113          776      16764

# 5. Plot alpha distribution to examine selection patterns
ggplot(results, aes(x = alpha)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  labs(x = "Alpha", y = "Count", 
       title = "Distribution of alpha values")

# Save neutral dataset
if(exists("snps_mlg.neutral")) {
  neutral_dataset_byscn <- snps_mlg.neutral[loc=-outliers]
  save(neutral_dataset_byscn, file="neutral_dataset_byscn.Rdata")
}
```

# Process Bayescan output to analyze and visualize outlier loci across populations:
1. Creates a frequency table from a genetic dataset
2. Filters outlier indices to ensure they're within valid range
3. Extracts the relevant columns from the frequency table
4. Creates a matrix of population-specific allele frequencies for the outlier loci
5. Calculates mean allele frequencies for each population
6. Generates a heatmap to visualize how these outlier loci vary across populations
7. Calculates and prints summary statistics (mean frequencies and most variable loci)
8. Saves the outlier summary data to a CSV file
```{r}
# Get the frequency table
freq_table <- tab(neutral_dataset_byscn, freq=TRUE)

# Convert outlier indices to column pairs and filter out any that are too large
max_locus <- ncol(freq_table)/2  # Maximum valid locus number
valid_outliers <- outliers[outliers <= max_locus]

# Convert valid outliers to column pairs
outlier_cols <- sort(c(valid_outliers * 2 - 1, valid_outliers * 2))

# Create matrix with outlier data
outlier_matrix <- freq_table[, outlier_cols]

# Verify our dimensions
print("Number of valid outliers:")
print(length(valid_outliers))
print("Dimensions of outlier matrix:")
print(dim(outlier_matrix))

# Reshape to get one row per locus
n_outliers <- length(valid_outliers)
pop_freqs <- matrix(nrow=n_outliers, ncol=nPop(neutral_dataset_byscn))
colnames(pop_freqs) <- levels(pop(neutral_dataset_byscn))

# Calculate population means
for(i in 1:nPop(neutral_dataset_byscn)) {
    pop_indices <- which(pop(neutral_dataset_byscn) == levels(pop(neutral_dataset_byscn))[i])
    # Take only first allele frequency for each locus
    pop_freqs[,i] <- colMeans(outlier_matrix[pop_indices, seq(1, ncol(outlier_matrix), 2)], na.rm=TRUE)
}

# Create heatmap
library(pheatmap)
pheatmap(pop_freqs,
         main = "Outlier Loci Allele Frequencies Across Populations",
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# Calculate and print summary statistics
print("Mean outlier allele frequencies per population:")
print(colMeans(pop_freqs))

outlier_var <- apply(pop_freqs, 1, var)
print("\nTop 10 most variable outlier loci:")
print(head(sort(outlier_var, decreasing=TRUE), 10))

# Save summary
outlier_summary <- data.frame(
    locus_id = valid_outliers,
    mean_freq = rowMeans(pop_freqs),
    var_freq = outlier_var
)
write.csv(outlier_summary, "outlier_summary.csv", row.names=FALSE)
```

# Identify loci under strong selection using a conservative threshold
Extract loci that have strong evidence of sesletion by filtering for log10Po values grater than 2 from Bayescan results 
log10PO>2 is a standard threshold in Bayescan analysis
Then extract neutral loci (loci not under selection) where log10PO<2
Print summary statistics showing how many loci fall into each category
```{r}
# Identify loci under strong selection using a conservative threshold
outliers <- which(bayescan_results$log10PO > 2)  # Common threshold is 2
neutral_loci <- which(bayescan_results$log10PO <= 2)

# Print summary
print(paste("Number of loci under selection:", length(outliers)))
print(paste("Number of neutral loci:", length(neutral_loci)))
```

#Verifying outliers
```{r}
# Let's look at the distribution of log10(PO) values
hist(bayescan_results$log10PO, 
     breaks=50, 
     main="Distribution of log10(PO) values",
     xlab="log10(PO)")

# Print summary statistics
summary(bayescan_results$log10PO)

# Let's see how many loci we get with different thresholds
thresholds <- c(1, 1.5, 2, 2.5, 3)
for(thresh in thresholds) {
    n_outliers <- sum(bayescan_results$log10PO > thresh)
    print(paste("Threshold", thresh, ":", n_outliers, "outliers"))
}
```

# Plot Bayescan output in R
```{r}
# 	 This file is used to plot figures for the software Bayescan in R.

#    This program, BayeScan, aims at detecting genetics markers under selection,
#	 based on allele frequency differences between population. 
#    Copyright (C) 2010  Matthieu Foll
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Arguments:
# - file is the name of your file ex: "output_fst.txt"
# - the q-value threshold corresponding to the target False Discovery Rate (FDR)
# - size is the size of the points and text labels for outliers
# - pos is the distance between the points and the labels 
# - highlight is a optional list of marker indices to display in red.
# - name_highlighted alows to write the indices of highlighted markers instead of using a point like the other markers
# - add_text adds the indices of the outlier markers

# Output:
# This function returns different paremeters in a list
# - outliers: the list of outliers
# - nb_outliers: the number of outliers

# Typical usage: 
# - load this file into R (file/source R code)
# - in R, go to the directory where "output_fst.txt" is (file/change current dir)
# - at the R prompt, type 
# > plot_bayescan("output_fst.txt",0,FDR=0.05)
# if you save the output in a variable, you can recall the different results:
# results<-plot_bayescan("output_fst.txt",0,FDR=0.05)
# results$outliers
# results$nb_outliers

#
# plotting posterior distribution is very easy in R with the output of BayeScan:
# first load the output file *.sel produced by BayeScan
# > mydata=read.table("bi.sel",colClasses="numeric")
# choose the parameter you want to plot by setting for example:
# parameter="Fst1"
# then this line will make the plot for:
# > plot(density(mydata[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))
# you can plot population specific Fst coefficient by setting
# parameter="Fst1"
# if you have non-codominant data you can plot posterior for Fis coefficients in each population:
# parameter="Fis1"
# if you test for selection, you can plot the posterior for alpha coefficient for selection:
# parameter="alpha1"
# you also have access to the likelihood with:
# parameter="logL"
# if you have the package "boa" installed, you can very easily obtain Highest Probability 
# Density Interval (HPDI) for your parameter of interest (example for the 95% interval):
# > boa.hpd(mydata[[parameter]],0.05)


plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,name_highlighted=F,add_text=T)
{
if (is.character(res))
  res=read.table(res)

colfstat=5
colq=colfstat-2

highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
non_highlight_rows=setdiff(1:nrow(res),highlight_rows)

outliers=as.integer(row.names(res[res[,colq]<=FDR,]))

ok_outliers=TRUE
if (sum(res[,colq]<=FDR)==0)
	ok_outliers=FALSE;

res[res[,colq]<=0.0001,colq]=0.0001

# plot
plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)

if (name_highlighted) {
 	if (length(highlight_rows)>0) {
 		text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
 	}
}
else {
	points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
	# add names of loci over p and vertical line
	if (ok_outliers & add_text) {
		text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
	}
}
lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)

return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}

```



```{r}
# First, source the plot_bayescan function (if you have it in a separate file)
# source("plot_bayescan.R")

# Or you can copy and paste the entire function definition into R

# Then use the function with your results file
results <- plot_bayescan("swiftia_bayescan_output_fst.txt", FDR=0.05)

# The function will:
# 1. Create a plot of FST vs log10(q-values)
# 2. Add a vertical line at your chosen FDR threshold (0.05)
# 3. Return a list containing:
#    - outliers: the list of outlier loci
#    - nb_outliers: the number of outliers

# You can access the results:
print(paste("Number of outliers:", results$nb_outliers))
print("Outlier loci:")
print(results$outliers)

# If you want to modify the plot appearance, you can adjust parameters:
results <- plot_bayescan("swiftia_bayescan_output_fst.txt", 
                      FDR=0.05,        # False Discovery Rate threshold
                     size=1,          # Point size
                         pos=0.35,        # Label position
                     add_text=TRUE)   # Add locus labels for outliers
```


```{r}
# Identify loci under selection using log(q-value) > 2 threshold
outlier_loci <- which(-log10(bayescan_results$qval) > 2)
outlier_loci
```

# Identify outliers with a more stringent approach
1. Create a composite score combining FST, log10(PO), and q-values
2. Filer for loci with strong evidence of selection (log10(PO)>2)
3. Select the top 30 outliers based on the composite score
4. Generate a visualization highlighting these outliers
5. Provide summary statistics for the selected loci
```{r}
# Read BayeScan results
results <- read.table("swiftia_bayescan_output_fst.txt", header=TRUE)
colnames(results) <- c("prob", "log10.PO", "qval", "alpha", "fst")

# Add row numbers as locus IDs
results$locus <- 1:nrow(results)

# Calculate -log10(q-value) for ranking
results$neg_log_qval <- -log10(results$qval)

# Create a composite score for ranking outliers
# Normalize each metric to 0-1 scale
results$fst_norm <- (results$fst - min(results$fst)) / (max(results$fst) - min(results$fst))
results$log10PO_norm <- (results$log10.PO - min(results$log10.PO)) / (max(results$log10.PO) - min(results$log10.PO))
results$qval_norm <- (results$neg_log_qval - min(results$neg_log_qval)) / (max(results$neg_log_qval) - min(results$neg_log_qval))

# Create composite score (equal weights)
results$composite_score <- (results$fst_norm + results$log10PO_norm + results$qval_norm) / 3

# Filter for strong candidates
strong_candidates <- results[results$log10.PO > 2 & results$alpha > 0, ]

# Sort by composite score
strong_candidates <- strong_candidates[order(-strong_candidates$composite_score), ]

# Select top 30 candidates
top_outliers <- head(strong_candidates, 30)

# Write results
write.csv(top_outliers, "top_30_outliers.csv", row.names=FALSE)

# Create summary plot
library(ggplot2)

ggplot(results, aes(x = neg_log_qval, y = fst)) +
  geom_point(aes(color = alpha > 0), alpha = 0.6) +
  scale_color_manual(values = c("blue", "red"),
                    labels = c("Balancing", "Diversifying"),
                    name = "Selection Type") +
  theme_bw() +
  labs(x = "-log10(q-value)",
       y = "FST",
       title = "BayeScan Results") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

# Print summary statistics for top outliers
cat("\nSummary of top 30 outliers:\n")
print(summary(top_outliers[c("fst", "log10.PO", "qval")]))
```

# Get summary statistics of outliers
```{r}
# Read BayeScan results
results <- read.table("swiftia_bayescan_output_fst.txt", header=TRUE)
colnames(results) <- c("prob", "log10.PO", "qval", "alpha", "fst")

# Add locus numbers
results$locus <- 1:nrow(results)

# Filter for high FST loci (> 0.3)
high_fst_loci <- results[results$fst > 0.3, ]

# Sort by FST value
high_fst_loci <- high_fst_loci[order(-high_fst_loci$fst), ]

# Create summary statistics
summary_stats <- data.frame(
  total_loci = nrow(high_fst_loci),
  mean_fst = mean(high_fst_loci$fst),
  min_fst = min(high_fst_loci$fst),
  max_fst = max(high_fst_loci$fst),
  mean_log10PO = mean(high_fst_loci$log10.PO),
  mean_qval = mean(high_fst_loci$qval)
)

# Write results to file
write.csv(high_fst_loci, "high_fst_loci.csv", row.names=FALSE)
write.csv(summary_stats, "high_fst_summary.csv", row.names=FALSE)

# Create visualization
library(ggplot2)

# Plot FST distribution for high FST loci
p1 <- ggplot(high_fst_loci, aes(x = fst)) +
  geom_histogram(bins = 20, fill = "darkred", alpha = 0.7) +
  theme_bw() +
  labs(x = "FST", y = "Count",
       title = "Distribution of FST Values > 0.3")

# Plot FST vs log10(PO) for high FST loci
p2 <- ggplot(high_fst_loci, aes(x = log10.PO, y = fst)) +
  geom_point(color = "darkred", size = 2) +
  theme_bw() +
  labs(x = "log10(PO)", y = "FST",
       title = "FST vs log10(PO) for High FST Loci")

# Print summary information
cat("\nAnalysis of High FST Loci (FST > 0.3):\n")
cat("\nNumber of loci:", nrow(high_fst_loci))
cat("\nFST range:", round(min(high_fst_loci$fst), 3), "to", round(max(high_fst_loci$fst), 3))
cat("\nMean FST:", round(mean(high_fst_loci$fst), 3))
cat("\nMean log10(PO):", round(mean(high_fst_loci$log10.PO), 3))

# Create sorted list of locus numbers
cat("\n\nLocus numbers for high FST variants (sorted by FST):\n")
print(high_fst_loci[order(-high_fst_loci$fst), c("locus", "fst", "log10.PO", "qval")])
```

# Plot with highlighting 18 outliers
```{r}
# Calculate -log10(q-value)
results$neg_log_qval <- -log10(results$qval)

# Create a new column to identify points in different regions
results$point_type <- case_when(
  results$fst > 0.35 ~ "high_outlier",  # Points in yellow box
  results$alpha > 0 ~ "diversifying",    # Other diversifying selection
  TRUE ~ "balancing"                     # Balancing selection
)

# Create the plot
plot <- ggplot(results, aes(x = neg_log_qval, y = fst)) +
  geom_point(aes(color = point_type), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("blue", "red", "orange"),
                    labels = c("Balancing", "Diversifying", "High Outliers"),
                    name = "Selection Type") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "-log10(q-value)",
       y = "FST",
       title = "BayeScan Results") +
  geom_vline(xintercept = 1.3, linetype = "dashed")

# Save
ggsave("BayeScan_Results_Outliers.pdf", plot = plot, device = "pdf", 
       path = "/File/Path/Feb18Output",
       scale = 1, width = 11, height = 10, units = "in", dpi = 300)
plot
```


# Create final neutral dataset
```{r}
# Create final neutral dataset
neutral_dataset_final <- snps_mlg.neutral[loc=-outliers_to_remove]
save(neutral_dataset_final, file="neutral_dataset_final.Rdata")

# Create neutral dataset metadata
# Extract individual names from the genind object
genind_samples_neu <- rownames(neutral_dataset_final@tab)

# Subset se_metadata to include only matching samples
se_metadata_neutral_all <- se_metadata[rownames(se_metadata) %in% genind_samples_neu, ]

# Create population map for downstream analyses
# Define all populations
pop_order_all <- c("Edisto", "AAR", "RTR", "Alderdice", "Geyer", "EFGB", "WFGB", "Bouma", "Elvers", 
               "Bright", "Diaphus", "McGrail", "FTR")

# Get the population for each sample
all_pops <- as.character(neutral_dataset_final@pop)

# Create subset_samples_all by selecting samples belonging to all populations
subset_samples_all <- genind_samples_neu[all_pops %in% pop_order_all]

# Check 
length(subset_samples_all)

# Create pop_code_all using subset_samples_all
pop_code_all <- data.frame(
  SAMPLE = subset_samples_all,
  STRATA = as.character(neutral_dataset_final@pop[match(subset_samples_all, rownames(neutral_dataset_final@tab))]),
  stringsAsFactors = FALSE
)

# Set row names for easier subsetting
rownames(pop_code_all) <- pop_code_all$SAMPLE

# Check the result
head(pop_code_all) 
#                     SAMPLE STRATA
# 14-v-18-1-001 14-v-18-1-001 Edisto
# 14-v-18-1-002 14-v-18-1-002 Edisto
# 14-v-18-1-004 14-v-18-1-004 Edisto
dim(pop_code_all)

# Create west to east order 
WE_order_all <- c("WFGB", "EFGB", "Bright", "Geyer", "Elvers", "McGrail", "Bouma", "Alderdice", "Diaphus", "FTR", "AAR", "RTR", "Edisto")

# Clean up
rm(subset_samples_all, all_pops, genind_samples_neu, pop_names_all)
```

# -----------Subset data for downstream analyses------------------
```{r}
# Subset with Edisto and all GoM populations with > 10 individuals (_sub)
snps_mlg_sub_n <- popsub(neutral_dataset_final, exclude=c('WFGB','Bright','Geyer'))
# Subset metadata (_sub)
se_metadata_neutral_sub <- subset(se_metadata_neutral, 
                                    !(SITE %in% c("WFGB", "Bright", "Geyer")))
save(se_metadata_neutral_sub, file='se_metadata_neutral_sub.Rdata')


# Subset without Edisto, all GoM populations (_GoMa)
snps_mlg_gom_all_n <- popsub(neutral_dataset_final, exclude = c('Edisto'))
save(snps_mlg_gom_all_n, file = 'snps_mlg_gom_all_n.Rdata')
# Subset metadata (_GoMa)
se_metadata_neutral_GoMa <- subset(se_metadata_neutral, SITE != "Edisto")
save(se_metadata_neutral_GoMa, file='se_metadata_neutral_GoMa.Rdata')


#Change all downstream snps_mlg_GoM10_n to snps_mlg_GoM10_n
# Subset without Edisto and populations with less than 10 individuals (_GoM10)
snps_mlg_GoM10_n <- popsub(neutral_dataset_final, exclude=c('WFGB','Bright','Geyer','Edisto'))
save(snps_mlg_GoM10_n, file='snps_mlg_GoM10_n.Rdata')
# Subset metadata (_GoM10)
se_metadata_neutral_GoM10 <- subset(se_metadata_neutral, 
                                    !(SITE %in% c("WFGB", "Bright", "Geyer", "Edisto")))
save(se_metadata_neutral_GoM10, file='se_metadata_neutral_GoM10.Rdata')


# Create population maps for downstream analyses
# Function to create pop_code objects from a genind object
create_pop_code <- function(genind_obj, output_name) {
  # Extract sample names from the genind object
  samples <- rownames(genind_obj@tab)
  
  # Get the population for each sample
  pops <- as.character(genind_obj@pop)
  
  # Create the pop_code data frame
  pop_code <- data.frame(
    SAMPLE = samples,
    STRATA = pops,
    stringsAsFactors = FALSE
  )
  
  # Set row names for easier subsetting
  rownames(pop_code) <- pop_code$SAMPLE
  
  # Print summary statistics
  cat(paste("\nCreated", output_name, "with", nrow(pop_code), "samples\n"))
  cat("Population breakdown:\n")
  print(table(pop_code$STRATA))
  
  # Return the pop_code object
  return(pop_code)
}

# 1 For snps_mlg_sub_n (Subset with Edisto and all GoM populations with > 10 individuals)
pop_code_sub <- create_pop_code(snps_mlg_sub_n, "pop_code_sub")
save(pop_code_sub, file="pop_code_sub.Rdata")
pop_order_sub <- c("Edisto", "AAR", "RTR", "Alderdice", "EFGB", "Bouma", "Elvers", "Diaphus", "McGrail", "FTR")
WE_order_sub <- c("EFGB", "Elvers", "McGrail", "Bouma", "Alderdice", "Diaphus", "FTR", "AAR", "RTR", "Edisto")

# 2 For snps_mlg_gom_all_n (Subset without Edisto, all GoM populations)
pop_code_GoMa <- create_pop_code(snps_mlg_gom_all_n, "pop_code_GoMa")
save(pop_code_GoMa, file="pop_code_GoMa.Rdata")
pop_order_GoMa <- c("AAR", "RTR", "Alderdice", "Geyer", "EFGB", "WFGB", "Bouma", "Elvers", "Bright", "Diaphus", "McGrail", "FTR")
WE_order_GoMa <- c("WFGB", "EFGB", "Bright", "Geyer", "Elvers", "McGrail", "Bouma", "Alderdice", "Diaphus", "FTR", "AAR", "RTR")

# 3 For snps_mlg_GoM10 (Subset without Edisto and populations with less than 10 individuals)
pop_code_GoM10 <- create_pop_code(snps_mlg_GoM10_n, "pop_code_GoM10")
save(pop_code_GoM10, file="pop_code_GoM10.Rdata")
pop_order_GoM10 <- c("AAR", "RTR", "Alderdice", "EFGB", "Bouma", "Elvers", "Diaphus", "McGrail", "FTR")
WE_order_GoM10 <- c("EFGB", "Elvers", "McGrail", "Bouma", "Alderdice", "Diaphus", "FTR", "AAR", "RTR")

# Print comparison of sample sizes
cat("\nComparison of sample sizes:\n")
cat("pop_code_sub: ", nrow(pop_code_sub), "samples\n")
cat("pop_code_GoMa: ", nrow(pop_code_GoMa), "samples\n")
cat("pop_code_GoM10: ", nrow(pop_code_GoM10), "samples\n")

#Comparison of sample sizes:
#pop_code_sub:  174 samples
#pop_code_GoMa:  176 samples
#pop_code_GoM10:  158 samples

# Save all objects into a single .Rdata file
save(
  # Dataset objects
  snps_mlg_sub_n,
  snps_mlg_gom_all_n,
  snps_mlg_GoM10,
  
  # Metadata objects
  se_metadata_neutral_sub,
  se_metadata_neutral_GoMa,
  se_metadata_neutral_GoM10,
  
  # Population code objects
  pop_code_sub,
  pop_code_GoMa,
  pop_code_GoM10,
  
  # Population order vectors
  pop_order_sub,
  WE_order_sub,
  pop_order_GoMa,
  WE_order_GoMa,
  pop_order_GoM10,
  WE_order_GoM10,
  
  # Function (optional, but can be useful)
  create_pop_code,
  
  file = "nonclone_subsetted_data.Rdata"
)
```

# Clean-up variables
```{r}
# Remove BayeScan results objects
rm(results)
rm(outliers_to_remove)
rm(outlier_freq_matrix)
rm(high_fst_loci)
rm(jittered_matrix)
rm(summary_stats)

# Remove visualization objects
rm(p1)
rm(p2)
rm(pca_result)
rm(pca_data)
rm(pop_info)
rm(pop_colors)

# Remove temporary data objects
rm(freq_mat)
rm(X)
rm(locus_map)
rm(outlier_prefixes)
rm(col_prefixes)
rm(outlier_genind)
```

#-------------------Summary of subsetted objects------------------
All populations (_all) 
192 Individuals
genid object: neutral_dataset_final
population mapping: pop_code_all
metadata: se_metadata_neutral_all
Pop order: pop_order_all ("Edisto", "AAR", "RTR", "Alderdice", "Geyer", "EFGB", "WFGB", "Bouma", "Elvers", 
               "Bright", "Diaphus", "McGrail", "FTR")
W2E order object: WE_order_all ("WFGB", "EFGB", "Bright", "Geyer", "Elvers", "McGrail", "Bouma", "Alderdice", "Diaphus", "FTR", "AAR", "RTR", "Edisto")

Edisto and GoM populations > 10 individuals (_sub)
174 samples
genid object: snps_mlg_sub_n 
population map: pop_code_sub
metadata: se_metadata_neutral_sub
pop order: pop_order_sub ("Edisto", "AAR", "RTR", "Alderdice", "EFGB", "Bouma", "Elvers", "Diaphus", "McGrail", "FTR")
W2E order object: WE_order_sub 

All GoM samples (_GoMa)
176 samples
genid object: snps_mlg_gom_all_n 
population map: pop_code_GoMa
metadata: se_metadata_neutral_GoMa
pop order: pop_order_GoMa ("AAR", "RTR", "Alderdice", "Geyer", "EFGB", "WFGB", "Bouma", "Elvers", "Bright", "Diaphus", "McGrail", "FTR")
W2E order object: WE_order_GoMa ("WFGB", "EFGB", "Bright", "Geyer", "Elvers", "McGrail", "Bouma", "Alderdice", "Diaphus", "FTR", "AAR", "RTR")

GoM populations greater than 10 individuals (_GoM10)
158 samples
snps_mlg_GoM10 → pop_code_GoM10
se_metadata_neutral_GoM10
pop order: pop_order_GoM10 ("AAR", "RTR", "Alderdice" "EFGB", "Bouma", "Elvers", "Diaphus", "McGrail", "FTR")
W2E order object: WE_order_GoM10 ("EFGB", "Elvers", "McGrail", "Bouma", "Alderdice", "Diaphus", "FTR", "AAR", "RTR")
