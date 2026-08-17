# Utility functions for monkey inbreeding analysis (R version)
# Common functions used across multiple phases

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(scales)  # For scales::comma in Manhattan plot
})

# ============================================================================
# CONFIGURATION LOADING
# ============================================================================

#' Load configuration from YAML file
#' 
#' @param config_path Path to config.yaml file
#' @return List containing configuration
load_config <- function(config_path = "config.yaml") {
  if (!file.exists(config_path)) {
    # Try in scripts directory
    # Use multiple methods to find script directory
    script_dir <- tryCatch({
      # Method 1: sys.frame
      if (!is.null(sys.frame(1)$ofile)) {
        dirname(sys.frame(1)$ofile)
      } else {
        # Method 2: current working directory
        getwd()
      }
    }, error = function(e) {
      # Fallback: use working directory
      getwd()
    })
    
    config_path <- file.path(script_dir, "config.yaml")
  }
  
  if (!file.exists(config_path)) {
    stop(sprintf("Configuration file not found: %s", config_path))
  }
  
  config <- yaml::read_yaml(config_path)
  return(config)
}

# ============================================================================
# LOGGING AND OUTPUT
# ============================================================================

#' Print log message with timestamp
#' 
#' @param message Message to log
#' @param level Log level (INFO, WARNING, ERROR)
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
}

#' Setup output directory
#' 
#' @param path Directory path
setup_output_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    log_message(sprintf("Created output directory: %s", path))
  }
}

# ============================================================================
# DATA LOADING UTILITIES
# ============================================================================

#' Read PLINK .fam file
#' 
#' @param fam_path Path to .fam file
#' @return Data frame with family information
read_plink_fam <- function(fam_path) {
  df <- fread(fam_path, header = FALSE, col.names = c("FID", "IID", "Father", "Mother", "Sex", "Phenotype"))
  return(df)
}

#' Read PLINK .bim file
#' 
#' @param bim_path Path to .bim file
#' @return Data frame with variant information
read_plink_bim <- function(bim_path) {
  df <- fread(bim_path, header = FALSE, col.names = c("CHR", "SNP", "cM", "POS", "A1", "A2"))
  return(df)
}

#' Read PLINK heterozygosity output
#' 
#' @param het_path Path to .het file
#' @return Data frame with heterozygosity values
read_plink_het <- function(het_path) {
  df <- fread(het_path)
  df$HET <- (df$`N(NM)` - df$`O(HOM)`) / df$`N(NM)`
  return(df)
}

#' Read PLINK ROH output
#' 
#' @param roh_path Path to .hom file
#' @return Data frame with ROH information
read_plink_roh <- function(roh_path) {
  df <- fread(roh_path)
  
  # Add ROH length category
  df$ROH_Category <- cut(
    df$KB * 1000,
    breaks = c(0, 1e6, 3e6, Inf),
    labels = c("short", "medium", "long"),
    include.lowest = TRUE
  )
  
  return(df)
}

# ============================================================================
# ROH ANALYSIS FUNCTIONS
# ============================================================================

#' Calculate F_ROH per individual
#' 
#' @param roh_df Data frame from PLINK .hom file
#' @param genome_length Autosomal genome length in bp
#' @return Data frame with F_ROH per individual
calculate_froh <- function(roh_df, genome_length = 2.8e9) {
  froh <- roh_df %>%
    group_by(IID) %>%
    summarise(
      Total_ROH_Length = sum(KB * 1000),
      Num_ROH = n(),
      F_ROH = Total_ROH_Length / genome_length,
      Mean_ROH_Length = mean(KB * 1000),
      Max_ROH_Length = max(KB * 1000)
    )
  
  return(froh)
}

#' Summarize ROH by length category
#' 
#' @param roh_df Data frame from PLINK .hom file
#' @return Summary statistics by ROH category
summarize_roh_categories <- function(roh_df) {
  summary <- roh_df %>%
    group_by(IID, ROH_Category) %>%
    summarise(
      Count = n(),
      Total_Length = sum(KB * 1000),
      Mean_Length = mean(KB * 1000),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = ROH_Category,
      values_from = c(Count, Total_Length, Mean_Length),
      values_fill = 0
    )
  
  return(summary)
}

# ============================================================================
# DIVERSITY METRICS FUNCTIONS
# ============================================================================

#' Calculate summary statistics for heterozygosity
#' 
#' @param het_df Data frame with heterozygosity values
#' @return Summary statistics
summarize_heterozygosity <- function(het_df) {
  summary <- data.frame(
    Mean_Het = mean(het_df$HET, na.rm = TRUE),
    SD_Het = sd(het_df$HET, na.rm = TRUE),
    Min_Het = min(het_df$HET, na.rm = TRUE),
    Max_Het = max(het_df$HET, na.rm = TRUE),
    Median_Het = median(het_df$HET, na.rm = TRUE),
    Q1_Het = quantile(het_df$HET, 0.25, na.rm = TRUE),
    Q3_Het = quantile(het_df$HET, 0.75, na.rm = TRUE)
  )
  
  return(summary)
}

#' Calculate summary statistics for nucleotide diversity
#' 
#' @param pi_df Data frame with pi values
#' @return Summary statistics
summarize_pi <- function(pi_df) {
  summary <- data.frame(
    Mean_PI = mean(pi_df$PI, na.rm = TRUE),
    SD_PI = sd(pi_df$PI, na.rm = TRUE),
    Min_PI = min(pi_df$PI, na.rm = TRUE),
    Max_PI = max(pi_df$PI, na.rm = TRUE),
    Median_PI = median(pi_df$PI, na.rm = TRUE)
  )
  
  return(summary)
}

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

#' Create custom ggplot2 theme for publication
#' 
#' @return ggplot2 theme object
theme_publication <- function() {
  theme_bw(base_size = 12, base_family = "Arial") +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title = element_text(face = "bold", size = 12),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key = element_blank(),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text = element_text(size = 10, face = "bold")
    )
}

#' Plot heterozygosity distribution
#' 
#' @param het_df Data frame with heterozygosity values
#' @param output_file Path to save plot
plot_heterozygosity <- function(het_df, output_file = NULL) {
  p <- ggplot(het_df, aes(x = HET)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(aes(xintercept = mean(HET)), color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = "Distribution of Individual Heterozygosity",
      x = "Observed Heterozygosity",
      y = "Count"
    ) +
    theme_publication()
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  }
  
  return(p)
}

#' Plot F_ROH distribution
#' 
#' @param froh_df Data frame with F_ROH values
#' @param output_file Path to save plot
plot_froh_distribution <- function(froh_df, output_file = NULL) {
  p <- ggplot(froh_df, aes(x = F_ROH)) +
    geom_histogram(bins = 20, fill = "coral", color = "black", alpha = 0.7) +
    geom_vline(aes(xintercept = mean(F_ROH)), color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = "Distribution of Genomic Inbreeding Coefficient",
      x = expression(F[ROH]),
      y = "Count"
    ) +
    theme_publication()
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  }
  
  return(p)
}

#' Plot ROH length distribution
#' 
#' @param roh_df Data frame from PLINK .hom file
#' @param output_file Path to save plot
plot_roh_length_distribution <- function(roh_df, output_file = NULL) {
  p <- ggplot(roh_df, aes(x = KB * 1000 / 1e6, fill = ROH_Category)) +
    geom_histogram(bins = 50, alpha = 0.7, color = "black") +
    scale_x_log10(labels = scales::comma) +
    scale_fill_manual(
      values = c("short" = "lightblue", "medium" = "orange", "long" = "darkred"),
      labels = c("Short (100kb-1Mb)", "Medium (1-3Mb)", "Long (>3Mb)")
    ) +
    labs(
      title = "Distribution of ROH Lengths",
      x = "ROH Length (Mb, log scale)",
      y = "Count",
      fill = "ROH Category"
    ) +
    theme_publication()
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  }
  
  return(p)
}

#' Manhattan plot for genomic statistics
#' 
#' @param data Data frame with CHR or CHROM, POS, and statistic columns
#' @param stat_col Name of statistic column to plot
#' @param output_file Path to save plot
plot_manhattan <- function(data, stat_col = "PI", output_file = NULL) {
  # Handle column name variations (CHROM vs CHR)
  if ("CHROM" %in% names(data) && !"CHR" %in% names(data)) {
    data$CHR <- data$CHROM
  }
  
  # Convert CHR to numeric for proper sorting, then back to character
  data$CHR_numeric <- as.numeric(as.character(data$CHR))
  
  # Sort chromosomes numerically (1, 2, 3, ..., 22)
  chr_levels <- sort(unique(data$CHR_numeric))
  
  # Prepare data with proper chromosome ordering
  data$CHR <- factor(data$CHR_numeric, levels = chr_levels, ordered = TRUE)
  data$CHR_num <- as.numeric(data$CHR)
  
  # Calculate chromosome positions for x-axis
  # Use equal spacing for each chromosome to avoid label overlap
  data <- data %>%
    arrange(CHR_numeric) %>%
    group_by(CHR) %>%
    mutate(
      n_windows_chr = n(),
      within_chr_pos = row_number()
    ) %>%
    ungroup()
  
  # Assign fixed width to each chromosome for equal spacing
  chr_info <- data %>%
    group_by(CHR) %>%
    summarize(
      n_windows_chr = first(n_windows_chr),
      .groups = "drop"
    ) %>%
    arrange(CHR) %>%
    mutate(
      # Use fixed width per chromosome (equal spacing)
      chr_width = max(n_windows_chr),  # Equal width based on largest chromosome
      chr_start = lag(cumsum(chr_width), default = 0),
      chr_end = chr_start + chr_width,
      chr_mid = chr_start + chr_width / 2
    )
  
  # Join back and calculate positions with equal chromosome spacing
  data <- data %>%
    left_join(chr_info %>% select(CHR, chr_start, chr_width), by = "CHR") %>%
    mutate(
      # Scale within-chromosome position to fit the fixed width
      cumulative_pos = chr_start + (within_chr_pos - 1) / n_windows_chr * chr_width
    )
  
  # Chromosome midpoints for labels
  chr_mids <- chr_info %>%
    select(CHR, mid = chr_mid)
  
  # Create labels vector - show all chromosomes, but use smaller font for many chromosomes
  n_chr <- length(chr_levels)
  chr_labels <- as.character(chr_mids$CHR)
  
  # Determine font size and plot width based on number of chromosomes
  if (n_chr > 20) {
    font_size <- 8
    plot_width <- 18
  } else if (n_chr > 15) {
    font_size <- 8.5
    plot_width <- 16
  } else {
    font_size <- 9
    plot_width <- 14
  }
  
  # Format labels for display (convert PI to π symbol)
  y_label <- if (stat_col == "PI") {
    expression(pi)
  } else {
    stat_col
  }
  
  title_text <- if (stat_col == "PI") {
    expression(paste("Genome-wide ", pi))
  } else {
    paste("Genome-wide", stat_col)
  }
  
  p <- ggplot(data, aes(x = cumulative_pos, y = .data[[stat_col]], color = CHR_num %% 2 == 0)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("steelblue", "coral"), guide = "none") +
    scale_x_continuous(
      labels = chr_labels,
      breaks = chr_mids$mid
    ) +
    labs(
      title = title_text,
      x = "Chromosome",
      y = y_label
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = font_size),
      plot.margin = margin(t = 10, r = 10, b = 60, l = 10, unit = "pt")
    )
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = plot_width, height = 6, dpi = 300)
  }
  
  return(p)
}

#' PCA biplot
#' 
#' @param pca_data Data frame with PC1, PC2, and sample IDs
#' @param output_file Path to save plot
plot_pca <- function(pca_data, output_file = NULL) {
  p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    labs(
      title = "Principal Component Analysis",
      x = sprintf("PC1 (%.1f%% variance)", attr(pca_data, "var_explained")[1]),
      y = sprintf("PC2 (%.1f%% variance)", attr(pca_data, "var_explained")[2])
    ) +
    theme_publication()
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 8, height = 8, dpi = 300)
  }
  
  return(p)
}

# ============================================================================
# STATISTICAL FUNCTIONS
# ============================================================================

#' Calculate confidence intervals
#' 
#' @param x Numeric vector
#' @param conf_level Confidence level (default 0.95)
#' @return Vector with mean, lower CI, upper CI
calculate_ci <- function(x, conf_level = 0.95) {
  n <- length(x)
  mean_x <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(n)
  margin <- qt((1 + conf_level) / 2, df = n - 1) * se
  
  return(c(
    mean = mean_x,
    lower = mean_x - margin,
    upper = mean_x + margin
  ))
}

#' Perform bootstrap resampling
#' 
#' @param data Data vector or data frame
#' @param statistic Function to calculate statistic
#' @param n_bootstrap Number of bootstrap replicates
#' @return Vector of bootstrap statistics
bootstrap <- function(data, statistic, n_bootstrap = 1000) {
  boot_stats <- replicate(n_bootstrap, {
    resample <- sample(1:nrow(data), replace = TRUE)
    statistic(data[resample, ])
  })
  
  return(boot_stats)
}

# ============================================================================
# DATA VALIDATION
# ============================================================================

#' Validate file exists
#' 
#' @param filepath Path to file
#' @param description File description
validate_file_exists <- function(filepath, description = "File") {
  if (!file.exists(filepath)) {
    stop(sprintf("%s not found: %s", description, filepath))
  }
}

#' Validate data frame columns
#' 
#' @param df Data frame
#' @param required_cols Required column names
validate_columns <- function(df, required_cols) {
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
  }
}

# ============================================================================
# OUTPUT FORMATTING
# ============================================================================

#' Format p-value for display
#' 
#' @param pval P-value
#' @param threshold Threshold for scientific notation
#' @return Formatted string
format_pvalue <- function(pval, threshold = 0.001) {
  if (pval < threshold) {
    return(sprintf("%.2e", pval))
  } else {
    return(sprintf("%.4f", pval))
  }
}

#' Write summary table
#' 
#' @param df Data frame
#' @param output_file Output file path
#' @param row_names Include row names
write_summary_table <- function(df, output_file, row_names = FALSE) {
  write.csv(df, output_file, row.names = row_names, quote = FALSE)
  log_message(sprintf("Summary table written to: %s", output_file))
}

# ============================================================================
# MAIN (for testing)
# ============================================================================

if (sys.nframe() == 0) {
  # Test configuration loading
  tryCatch({
    config <- load_config()
    log_message("Configuration loaded successfully")
    log_message(sprintf("Project directory: %s", config$project$base_dir))
  }, error = function(e) {
    log_message(sprintf("Error loading configuration: %s", e$message), "ERROR")
  })
}


