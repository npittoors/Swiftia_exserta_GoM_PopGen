# =============================================================================
# 02_process_results.R
# Swiftia exserta population genomics — ADMIXTURE results processing
#
# Inputs:
#   - pop_data_admixture.Rdata              (from 01_data_prep.R)
#   - se_metadata_neutral (loaded from your main metadata object)
#   - ADMIXTURE log files in admixture_results_<dataset>/K{k}/rep{rep}/
#   - Best-run .Q files in ADMIXTURE_<dataset>/best_run_files/
#
# Outputs:
#   - CV error plots (mean ± range across replicates) per dataset
#   - Faceted structure bar plots (west → east) per dataset and K
#   - Saved PDFs in Feb18Output/
#
# Datasets: Sub (Gulf-Atlantic), GoM10, GoMa
# =============================================================================

library(ggplot2)
library(tidyr)
library(adegenet)

# ---- Paths ------------------------------------------------------------------
base_path    <- "update/base/path"
out_path     <- file.path(base_path, "update/out/path")
dir.create(out_path, showWarnings = FALSE)

sub_results  <- file.path(base_path, "admixture_results")           # Sub log files
gom10_results <- file.path(base_path, "ADMIXTURE_GoM10/admixture_results_GoM10")
goma_results  <- file.path(base_path, "admixture_results_GoMa")

sub_best     <- file.path(base_path, "admixture_results/best_run_files")
gom10_best   <- file.path(base_path, "ADMIXTURE_GoM10/best_run_files")
goma_best    <- file.path(base_path, "admixture_results_GoMa/best_run_files")

# ---- Load pop_data ----------------------------------------------------------
load(file.path(base_path, "pop_data_admixture.Rdata"))
# Loads: pop_data_sub, pop_data_GoM10, pop_data_GoMa

# ---- Load metadata ----------------------------------------------------------
# se_metadata_neutral must be loaded from your main workspace
# Rows are individual IDs; must contain column LONG (longitude)
# Subset to each dataset after confirming individual order:
se_metadata_neutral_GoM10 <- se_metadata_neutral[pop_data_GoM10$individual, ]
se_metadata_neutral_GoMa  <- se_metadata_neutral[pop_data_GoMa$individual, ]
se_metadata_neutral_sub   <- se_metadata_neutral[pop_data_sub$individual, ]

# Sanity checks — all must return TRUE before proceeding
stopifnot(all(pop_data_GoM10$individual == rownames(se_metadata_neutral_GoM10)))
stopifnot(all(pop_data_GoMa$individual  == rownames(se_metadata_neutral_GoMa)))
stopifnot(all(pop_data_sub$individual   == rownames(se_metadata_neutral_sub)))


# =============================================================================
# FUNCTION: Parse CV errors from ADMIXTURE log files
# Reads all K × rep log files from a results directory and returns a summary
# dataframe with mean, min, and max CV error per K.
# =============================================================================
parse_cv_errors <- function(results_dir, k_max = 10, n_reps = 10) {
  cv_data <- data.frame(K = integer(), CV_Error = numeric())
  
  for (k in 1:k_max) {
    for (rep in 1:n_reps) {
      log_file <- file.path(results_dir,
                            paste0("K", k),
                            paste0("rep", rep),
                            paste0("log", k, "_", rep, ".out"))
      if (file.exists(log_file)) {
        log_content <- readLines(log_file)
        cv_line     <- grep("CV error", log_content, value = TRUE)
        if (length(cv_line) > 0) {
          cv_error <- as.numeric(gsub(".*: ", "", cv_line))
          cv_data  <- rbind(cv_data, data.frame(K = k, CV_Error = cv_error))
        }
      }
    }
  }
  
  if (nrow(cv_data) == 0) {
    stop("No CV error data found. Confirm ADMIXTURE was run with --cv flag and log files exist.")
  }
  
  cv_summary <- aggregate(CV_Error ~ K, data = cv_data,
                          FUN = function(x) c(mean = mean(x), min = min(x), max = max(x)))
  cv_summary <- do.call(data.frame, cv_summary)
  names(cv_summary) <- c("K", "Mean", "Min", "Max")
  cv_summary <- cv_summary[order(cv_summary$K), ]
  
  return(cv_summary)
}


# =============================================================================
# FUNCTION: Plot CV errors
# Plots mean CV error ± range (shaded ribbon) across replicates by K.
# =============================================================================
plot_cv_errors <- function(cv_summary, dataset_label = "", save_plot = FALSE,
                           out_file = NULL) {
  min_k <- cv_summary$K[which.min(cv_summary$Mean)]
  
  p <- ggplot(cv_summary, aes(x = K, y = Mean)) +
    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = 0.2, fill = "#3491A8FF") +
    geom_line(color = "grey30") +
    geom_point(size = 3, color = "grey20") +
    geom_point(data = cv_summary[cv_summary$K == min_k, ],
               aes(x = K, y = Mean), color = "#3491A8FF", size = 4) +
    theme_minimal() +
    labs(title    = paste0("Cross-validation error by K", 
                           if (nchar(dataset_label) > 0) paste0(" (", dataset_label, ")") else ""),
         subtitle = paste0("Optimal K = ", min_k, " (lowest mean CV error)"),
         x        = "Number of clusters (K)",
         y        = "Cross-validation error") +
    scale_x_continuous(breaks = cv_summary$K) +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  cat("\nCV Error summary —", dataset_label, "\n")
  print(cv_summary)
  cat("Optimal K =", min_k, "\n\n")
  print(p)
  
  if (save_plot && !is.null(out_file)) {
    ggsave(out_file, p, width = 6, height = 4, dpi = 300)
    cat("Saved:", out_file, "\n")
  }
  
  return(invisible(list(plot = p, data = cv_summary, optimal_k = min_k)))
}


# =============================================================================
# FUNCTION: Plot ADMIXTURE structure bar plot
# Faceted west → east by longitude, one bar per individual.
#
# IMPORTANT — label switching:
#   ADMIXTURE column assignments (K1, K2...) are arbitrary across runs.
#   Before finalizing figures, confirm which column corresponds to which
#   population cluster by inspecting mean ancestry in Edisto individuals
#   (Sub dataset) or EFGB individuals (GoM10/GoMa).
# =============================================================================
plot_admixture_structure <- function(q_file, pop_data, metadata, k,
                                     pop_order, dataset_label = "",
                                     colors = NULL, facet = TRUE,
                                     save_plot = FALSE, out_file = NULL) {
  # Default color palette (teal, light green, deep purple)
  if (is.null(colors)) {
    colors <- c("#3491A8FF", "#C3E9CEFF", "#35264CFF",
                "#FDE725FF", "#453781FF", "#20A387FF")
  }
  colors <- colors[1:k]
  
  q_matrix <- read.table(q_file)
  
  if (nrow(q_matrix) != nrow(pop_data)) {
    stop("Q matrix rows (", nrow(q_matrix), ") do not match pop_data rows (",
         nrow(pop_data), "). Check that pop_data matches the dataset used for this run.")
  }
  
  adm_plot_data <- data.frame(
    Sample     = pop_data$individual,
    Population = pop_data$population,
    Longitude  = metadata$LONG,
    q_matrix
  )
  names(adm_plot_data)[names(adm_plot_data) %in% paste0("V", 1:k)] <- paste0("K", 1:k)
  
  # Filter to populations present in pop_order; warn about any missing
  adm_plot_data <- adm_plot_data[adm_plot_data$Population %in% pop_order, ]
  missing_pops  <- setdiff(pop_order, unique(adm_plot_data$Population))
  if (length(missing_pops) > 0) {
    warning("Populations in pop_order not found in data: ",
            paste(missing_pops, collapse = ", "))
  }
  
  adm_plot_data$Population <- factor(adm_plot_data$Population, levels = pop_order)
  adm_plot_data <- adm_plot_data[order(adm_plot_data$Population, adm_plot_data$Longitude), ]
  adm_plot_data$Position <- 1:nrow(adm_plot_data)
  
  plot_data_long <- gather(adm_plot_data, key = "Ancestry", value = "Proportion",
                           starts_with("K"))
  plot_data_long$Ancestry <- factor(plot_data_long$Ancestry, levels = paste0("K", 1:k))
  
  p <- ggplot(plot_data_long, aes(x = Position, y = Proportion, fill = Ancestry)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = colors,
                      labels = paste0("Cluster ", 1:k)) +
    theme_minimal() +
    labs(title  = if (nchar(dataset_label) > 0) paste0(dataset_label, " — K = ", k) else paste0("K = ", k),
         y      = "Ancestry proportion",
         x      = NULL,
         fill   = "Genetic cluster") +
    scale_y_continuous(expand = c(0, 0))
  
  if (facet) {
    p <- p + facet_grid(~Population, scales = "free_x", space = "free_x")
  }
  
  p <- p + theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    panel.spacing   = unit(0.1, "lines"),
    strip.text.x    = element_text(size = 10, angle = 90),
    plot.margin     = margin(t = 40, r = 20, b = 20, l = 20, unit = "pt"),
    axis.title.y    = element_text(size = 12),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  )
  
  print(p)
  
  if (save_plot && !is.null(out_file)) {
    ggsave(out_file, p, width = 12, height = 5, dpi = 300)
    cat("Saved:", out_file, "\n")
  }
  
  return(invisible(list(
    plot = p,
    data = adm_plot_data[, c("Sample", "Population", "Longitude", "Position")]
  )))
}


# =============================================================================
# POPULATION ORDER VECTORS
# West → east for each dataset
# Note: Edisto is an Atlantic outgroup and must not appear on Gulf maps;
# it is included in Sub only, with a separate notation in figure captions.
# =============================================================================
pop_order_sub   <- c('WFGB', 'EFGB', 'Bright', 'Geyer', 'Elvers', 'McGrail',
                     'Bouma', 'Alderdice', 'Diaphus', 'FTR', 'AAR', 'RTR', 'Edisto')
pop_order_GoM10 <- c('WFGB', 'EFGB', 'Bright', 'Geyer', 'Elvers', 'McGrail',
                     'Bouma', 'Alderdice', 'Diaphus', 'FTR', 'AAR', 'RTR')
pop_order_GoMa  <- pop_order_GoM10


# =============================================================================
# RUN: Sub (Gulf-Atlantic dataset)
# Optimal K = 2 (Gulf vs. Atlantic)
# =============================================================================

cv_sub <- parse_cv_errors(sub_results)

cv_results_sub <- plot_cv_errors(
  cv_sub,
  dataset_label = "Sub (Gulf-Atlantic)",
  save_plot     = TRUE,
  out_file      = file.path(out_path, "cv_error_Sub.pdf")
)

# ---- Structure plot K=2 ----
# NOTE: Confirm which column (V1 or V2) has higher mean ancestry in Edisto
# before accepting the color-to-cluster mapping below.
adm_sub_K2 <- plot_admixture_structure(
  q_file        = file.path(sub_best, "swiftia_neutral.2.Q"),
  pop_data      = pop_data_sub,
  metadata      = se_metadata_neutral_sub,
  k             = 2,
  pop_order     = pop_order_sub,
  dataset_label = "Sub",
  save_plot     = TRUE,
  out_file      = file.path(out_path, "admixture_Sub_K2_faceted.pdf")
)


# =============================================================================
# RUN: GoM10 (Gulf of Mexico only, 158 ind)
# Optimal K = 1 (monotonically increasing CV; consistent with IBD)
# K=2 retained for supplementary
# =============================================================================

cv_GoM10 <- parse_cv_errors(gom10_results)

cv_results_GoM10 <- plot_cv_errors(
  cv_GoM10,
  dataset_label = "GoM10",
  save_plot     = TRUE,
  out_file      = file.path(out_path, "cv_error_GoM10.pdf")
)

# ---- Structure plot K=2 (supplementary) ----
adm_GoM10_K2 <- plot_admixture_structure(
  q_file        = file.path(gom10_best, "swiftia_GoM10.2.Q"),
  pop_data      = pop_data_GoM10,
  metadata      = se_metadata_neutral_GoM10,
  k             = 2,
  pop_order     = pop_order_GoM10,
  dataset_label = "GoM10",
  save_plot     = TRUE,
  out_file      = file.path(out_path, "admixture_GoM10_K2_faceted.pdf")
)


# =============================================================================
# RUN: GoMa (GoM with non-clonal singletons re-added, 176 ind)
# =============================================================================

cv_GoMa <- parse_cv_errors(goma_results)

cv_results_GoMa <- plot_cv_errors(
  cv_GoMa,
  dataset_label = "GoMa",
  save_plot     = TRUE,
  out_file      = file.path(out_path, "cv_error_GoMa.pdf")
)

# ---- Structure plot K=2 ----
adm_GoMa_K2 <- plot_admixture_structure(
  q_file        = file.path(goma_best, "swiftia_GoMa.2.Q"),
  pop_data      = pop_data_GoMa,
  metadata      = se_metadata_neutral_GoMa,
  k             = 2,
  pop_order     = pop_order_GoMa,
  dataset_label = "GoMa",
  save_plot     = TRUE,
  out_file      = file.path(out_path, "admixture_GoMa_K2_faceted.pdf")
)


# =============================================================================
# OPTIONAL: Check individual plotting order for a given dataset
# Use to verify west-east arrangement before finalizing figures
# =============================================================================
check_plot_order <- function(q_file, pop_data, metadata, pop_order) {
  q_matrix  <- read.table(q_file)
  plot_data <- data.frame(
    Sample     = pop_data$individual,
    Population = pop_data$population,
    Longitude  = metadata$LONG,
    q_matrix
  )
  plot_data <- plot_data[plot_data$Population %in% pop_order, ]
  plot_data$Population <- factor(plot_data$Population, levels = pop_order)
  plot_data <- plot_data[order(plot_data$Population, plot_data$Longitude), ]
  
  cat("Samples ordered by population (W→E) and longitude:\n")
  print(data.frame(Order      = 1:nrow(plot_data),
                   Sample     = plot_data$Sample,
                   Population = plot_data$Population,
                   Longitude  = round(plot_data$Longitude, 3)))
  
  cat("\nLongitude ranges by population:\n")
  long_ranges <- tapply(plot_data$Longitude, plot_data$Population,
                        function(x) paste(round(range(x), 3), collapse = " to "))
  print(long_ranges)
  
  cat("\nSamples per population:\n")
  print(table(plot_data$Population))
  
  return(invisible(plot_data))
}

# Example usage:
# order_check <- check_plot_order(
#   q_file    = file.path(gom10_best, "swiftia_GoM10.2.Q"),
#   pop_data  = pop_data_GoM10,
#   metadata  = se_metadata_neutral_GoM10,
#   pop_order = pop_order_GoM10
# )