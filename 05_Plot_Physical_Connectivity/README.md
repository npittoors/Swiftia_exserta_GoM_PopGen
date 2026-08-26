# Physical Connectivity Plotting — 05_Plot_Physical_Connectivity.Rmd

This script plots the biophysical larval dispersal model output (Ichthyop/CROCO)
for *Swiftia exserta* in the northern Gulf of Mexico. It is
standalone — it does not depend on any object from the population genetics
scripts (`01`-`04`), only on the raw connectivity CSVs.

## What this script does

1. **Import connectivity data** — reads the release-by-recruitment ensemble
   matrix and its corresponding year-of-connection matrix, and parses both
   into a consistent release-order x recruitment-order layout.
2. **GoM10 subset and fixes** — converts `"<0.1"` floor labels to a numeric
   value, fixes RTR's row label (misread as `"NA"` from the source CSV),
   adds an RTR recruitment column of zeros if absent, and subsets to the 9
   GoM10 populations used elsewhere in the pipeline (drops WFGB, Bright,
   Geyer).
3. **Figure 1 — Connectivity matrix** — full 12-population ensemble matrix,
   magma-filled by connectivity value, with point shape indicating release
   year(s) (2014 / 2015 / Both).
4. **Figure 2 prep** — builds the directed edge list (release -> recruit,
   nonzero values only), splits it into main (>=0.1% transport) and trace
   (<0.1%, nonzero) edge sets, assigns arrow curvature and magma coloring,
   and pulls the Gulf of Mexico bathymetry/coastline basemap.
5. **Figure 2 — Connectivity map** — GoM10 spatial network map: main edges
   drawn as magma-colored solid arrows sized/colored by exchange strength;
   trace edges drawn as thin dashed grey arrows.

## Requirements

R packages used across this script:

tidyverse (dplyr, tidyr, readr, purrr, forcats), viridis, scales, marmap, maps

## Required input files

- `Se_popgen_Connectivity_matrix.csv` — larval settlement proportions (%),
  release zones (rows) x recruitment zones (columns)
- `Se_popgen_Connectivity_matrix_years.csv` — same dimensions; release
  year(s) ("2014", "2015", "Both") per cell, NA where no modeled connection
  exists

## Outputs

**Figures (PDF/PNG):**
- `figures/phys_connectivity_matrix.pdf/png` — full 12-population ensemble
  connectivity matrix
- `figures/swiftia_GoM10_connectivity_map.pdf/png` — 9-population GoM10
  directed network map

## Notes

- **RTR fix:** the source CSV has RTR's row label read in as `NA` rather
  than `"RTR"`. This is fixed programmatically (not by editing the source
  CSV) so the fix is visible and reproducible on rerun.
- **`<0.1` values:** treated numerically as `0.05` for fill color in the
  matrix figure. For the map figure, `<0.1` (but nonzero) edges are kept as
  a separate "trace exchange" layer rather than dropped entirely, so a
  reader can distinguish "no modeled connection" from "a very weak modeled
  connection."
- **GoM10 subset:** the map uses the same 9 populations as the GoM10
  grouping used throughout the population genetics pipeline, for direct
  comparison with the genomic co-ancestry results in `04_fineRADstructure_p.Rmd`.
- **Manual tuning still needed:** arrow curvature for the trace-exchange
  layer and per-site label offsets (`label_df`) are set to reasonable
  uniform defaults but were originally tuned by eye in earlier drafts —
  check that no arrow or label overlaps another element once you regenerate
  the figure, and adjust the `curve`/`label_hjust`/`label_vjust` values in
  Section 4 if so.
- This script only plots the larval dispersal model output on its own — it
  does not overlay genomic co-ancestry (FST or fineRADstructure) data.
