#!/bin/sh

# Create directory for admixture results
mkdir -p admixture_results_GoM10

# Loop through K values 1-10 (using seq instead of {1..10})
for K in $(seq 1 10)
do
    # Create directory for this K value
    mkdir -p admixture_results_Gom10/K$K
    
    # Run 10 replicates for each K
    for rep in $(seq 1 10)
    do
        mkdir -p admixture_results_Gom10/K$K/rep$rep
        
        # Run ADMIXTURE with 100 threads and cross-validation
        # Using a fixed seed for reproducibility
        SEED=$(expr $K \* 100 + $rep)
        
        admixture -j100 --cv sswiftia_GoM10.bed $K \
          -s $SEED > admixture_results_GoM10/K$K/rep$rep/log${K}_${rep}.out
    done
done

echo "ADMIXTURE GoM10 analysis complete!"
