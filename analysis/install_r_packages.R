#!/usr/bin/env Rscript
# Install CRAN packages used by phase 3a plotting and analysis/utils.R

cran_packages <- c(
  "yaml",
  "dplyr",
  "tidyr",
  "data.table",
  "ggplot2",
  "scales",
  "patchwork"
)

cat("Installing R packages for phase 3a plotting...\n\n")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

cat("\nVerifying installations:\n")
all_installed <- TRUE
for (pkg in cran_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  OK  %s\n", pkg))
  } else {
    cat(sprintf("  FAIL  %s\n", pkg))
    all_installed <- FALSE
  }
}

if (all_installed) {
  cat("\nAll packages installed successfully.\n")
} else {
  cat("\nSome packages failed to install. Check the errors above.\n")
  quit(status = 1)
}
