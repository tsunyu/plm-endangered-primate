#!/usr/bin/env Rscript
# Phase 3a Summary Plot Generation Script

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

analysis_dir <- Sys.getenv("PLM_ANALYSIS_DIR")
base_dir <- Sys.getenv("PLM_BASE_DIR")
if (!nzchar(analysis_dir)) {
  file_arg <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
  analysis_dir <- normalizePath(file.path(dirname(file_arg), ".."))
}
if (!nzchar(base_dir)) {
  stop("Set PLM_BASE_DIR or run: bash configure_base_dir.sh /path/to/analysis_root")
}

# Load utility functions
source(file.path(analysis_dir, "utils.R"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  output_dir <- args[1]
} else {
  output_dir <- file.path(base_dir, "output/phase3a_population_genomics/plots")
}

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Generating Phase 3a summary plots...\n")
cat("Output directory:", output_dir, "\n")

# ============================================================================
# 1. ROH Analysis Summary
# ============================================================================

roh_summary_file <- file.path(base_dir, "output/phase3a_population_genomics/roh_analysis/roh_summary_per_individual.csv")
if (file.exists(roh_summary_file)) {
  cat("Loading ROH data...\n")
  roh_stats <- read.csv(roh_summary_file)
  
  # F_ROH distribution with better formatting
  p_froh <- ggplot(roh_stats, aes(x = F_ROH)) +
    geom_histogram(bins = 20, fill = "coral", color = "black", alpha = 0.7) +
    geom_vline(aes(xintercept = mean(F_ROH)), color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = "Distribution of Genomic Inbreeding Coefficient",
      x = expression(F[ROH]),
      y = "Count"
    ) +
    theme_publication()
  
  ggsave(file.path(output_dir, "summary_froh_distribution.png"), p_froh, 
         width = 8, height = 6, dpi = 300)
  
  # Number of ROH per individual with better x-axis
  p_roh_count <- ggplot(roh_stats, aes(x = reorder(IID, Num_ROH), y = Num_ROH)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    labs(
      title = "Number of ROH per Individual",
      x = "Individual",
      y = "Number of ROH"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  ggsave(file.path(output_dir, "summary_roh_count_per_individual.png"), p_roh_count, 
         width = 12, height = 6, dpi = 300)
  
  cat("ROH summary plots generated\n")
}

# ============================================================================
# 2. Diversity Metrics Summary
# ============================================================================

het_file <- file.path(base_dir, "output/phase3a_population_genomics/diversity_metrics/individual_heterozygosity.txt")
if (file.exists(het_file)) {
  cat("Loading heterozygosity data...\n")
  het_data <- read.table(het_file, header = TRUE)
  
  # Heterozygosity distribution
  p_het <- ggplot(het_data, aes(x = OBS_HET)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(aes(xintercept = mean(OBS_HET)), color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = "Distribution of Observed Heterozygosity",
      x = "Observed Heterozygosity",
      y = "Count"
    ) +
    theme_publication()
  
  ggsave(file.path(output_dir, "summary_heterozygosity_distribution.png"), p_het, 
         width = 8, height = 6, dpi = 300)
  
  cat("Heterozygosity summary plot generated\n")
}

# ============================================================================
# 3. Population Structure Summary
# ============================================================================

pca_file <- file.path(base_dir, "output/phase3a_population_genomics/population_structure/pca_results.txt")
if (file.exists(pca_file)) {
  cat("Loading PCA data...\n")
  pca_data <- read.table(pca_file, header = TRUE)
  
  # PCA plot
  p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    labs(
      title = "Principal Component Analysis",
      x = paste0("PC1 (", round(var(pca_data$PC1) / sum(apply(pca_data[,1:min(10, ncol(pca_data))], 2, var)) * 100, 1), "%)"),
      y = paste0("PC2 (", round(var(pca_data$PC2) / sum(apply(pca_data[,1:min(10, ncol(pca_data))], 2, var)) * 100, 1), "%)")
    ) +
    theme_publication()
  
  ggsave(file.path(output_dir, "summary_pca_plot.png"), p_pca, 
         width = 8, height = 6, dpi = 300)
  
  cat("PCA summary plot generated\n")
}

# ============================================================================
# 4. Effective Population Size Summary
# ============================================================================

ne_file <- file.path(base_dir, "output/phase3a_population_genomics/ne_estimation/ld_based/ne_estimates.txt")
if (file.exists(ne_file)) {
  cat("Loading Ne estimation data...\n")
  ne_data <- read.table(ne_file, header = TRUE)
  
  # Ne over time plot
  p_ne <- ggplot(ne_data, aes(x = Generation, y = Ne)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "steelblue", size = 2) +
    scale_y_log10(labels = scales::comma) +
    labs(
      title = "Effective Population Size Over Time",
      x = "Generation",
      y = "Effective Population Size (log scale)"
    ) +
    theme_publication()
  
  ggsave(file.path(output_dir, "summary_ne_over_time.png"), p_ne, 
         width = 10, height = 6, dpi = 300)
  
  cat("Ne estimation summary plot generated\n")
}

# ============================================================================
# 5. Combined Summary Dashboard
# ============================================================================

cat("Generating combined summary dashboard...\n")

# Create a combined plot if we have multiple datasets
plots_list <- list()

if (exists("p_froh")) plots_list[["F_ROH"]] <- p_froh
if (exists("p_het")) plots_list[["Heterozygosity"]] <- p_het
if (exists("p_pca")) plots_list[["PCA"]] <- p_pca
if (exists("p_ne")) plots_list[["Ne"]] <- p_ne

if (length(plots_list) >= 2) {
  # Combine plots using patchwork
  combined_plot <- wrap_plots(plots_list, ncol = 2)
  
  ggsave(file.path(output_dir, "phase3a_summary_dashboard.png"), combined_plot, 
         width = 16, height = 12, dpi = 300)
  
  cat("Combined summary dashboard generated\n")
}

cat("Phase 3a summary plots generation completed!\n")
cat("Plots saved to:", output_dir, "\n")
