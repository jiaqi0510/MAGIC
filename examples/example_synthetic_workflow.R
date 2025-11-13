#!/usr/bin/env Rscript
#
# MAGIC: Methylation Analysis with Genomic Inferred Contexts
#
# Example: Complete workflow with synthetic data
# Description: Demonstrates the entire MAGIC workflow using synthetic data:
#              1. Generate synthetic data with known parameters
#              2. Optimize mixture model
#              3. Compare fitted vs true parameters
#              4. Evaluate parameter recovery
#
#              This is useful for:
# - Testing the software installation
# - Understanding the complete workflow
# - Validating the method
# - Benchmarking performance

# ==============================================================================
# Configuration
# ==============================================================================

# Set parameters for synthetic data
K <- 3                    # Number of components
n_cpgs <- 50000          # Number of CpG sites (reduced for speed)
n_samples <- 50          # Number of samples
seed <- 12345            # Random seed for reproducibility

# Create output directory
output_dir <- "synthetic_workflow_test"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("MAGIC Complete Workflow Example\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("This example will:\n")
cat("  1. Generate synthetic data (K=3, 50K CpGs, 50 samples)\n")
cat("  2. Fit mixture model\n")
cat("  3. Compare fitted vs true parameters\n")
cat(sprintf("\nOutput directory: %s/\n\n", output_dir))

# ==============================================================================
# Step 1: Generate Synthetic Data
# ==============================================================================

cat("Step 1: Generating synthetic data\n")
cat(paste(rep("-", 70), collapse = ""), "\n\n")

synthetic_file <- file.path(output_dir, sprintf("synthetic_K%d.qs", K))

system2("Rscript", c(
  "scripts/generate_synthetic_data.R",
  "-K", as.character(K),
  "-n", as.character(n_cpgs),
  "-s", as.character(n_samples),
  "-o", synthetic_file,
  "--seed", as.character(seed)
))

cat("\n")

# Check that files were created
if (!file.exists(synthetic_file)) {
  stop("Failed to generate synthetic data")
}

true_model_file <- file.path(output_dir, sprintf("synthetic_K%d_trueModel.csv", K))
if (!file.exists(true_model_file)) {
  stop("True model file not created")
}

cat("✓ Synthetic data generated successfully\n")
cat(sprintf("  Data: %s\n", synthetic_file))
cat(sprintf("  True model: %s\n\n", true_model_file))

# ==============================================================================
# Step 2: Fit Mixture Model
# ==============================================================================

cat("\nStep 2: Fitting mixture model\n")
cat(paste(rep("-", 70), collapse = ""), "\n\n")

cat(sprintf("Optimizing K=%d mixture model with 5 convergence runs\n", K))
cat("This may take a few minutes...\n\n")

system2("Rscript", c(
  "R/magicFit.R",
  "--input", synthetic_file,
  "--k_components", as.character(K),
  "--max_iters", "500",
  "--tol", "1e-6",
  "--threads", "4",
  "--n_conv_runs", "5"
))

cat("\n")

# Find the optimization output directory
# The optimization creates the directory in the current working directory
opt_dirs <- list.dirs(".", recursive = FALSE)
opt_dirs <- opt_dirs[grepl(sprintf("_K%d_", K), opt_dirs)]

if (length(opt_dirs) == 0) {
  stop("Optimization output directory not found")
}

# Use most recent
opt_dir <- opt_dirs[length(opt_dirs)]

cat("✓ Model optimization complete\n")
cat(sprintf("  Results: %s\n\n", opt_dir))

# ==============================================================================
# Step 3: Compare Parameters
# ==============================================================================

cat("\nStep 3: Comparing fitted vs true parameters\n")
cat(paste(rep("-", 70), collapse = ""), "\n\n")

comparison_csv <- file.path(output_dir, "parameter_comparison.csv")

system2("Rscript", c(
  "scripts/model_comparison.R",
  "--true-dir", output_dir,
  "--fitted-dir", ".",  # Optimization output is in current directory
  "--output", comparison_csv
))

cat("\n")

cat("✓ Parameter comparison complete\n")
cat(sprintf("  Summary: %s\n\n", comparison_csv))

# ==============================================================================
# Step 4: Display Results
# ==============================================================================

cat("\nStep 4: Results Summary\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Read true parameters
true_model <- read.csv(true_model_file)
cat("True Model Parameters:\n")
cat(sprintf("%-12s %-10s %-10s %-10s %-10s\n", 
            "Component", "Alpha", "Beta", "Pi", "Mean Meth"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in 1:nrow(true_model)) {
  mean_meth <- true_model$alpha[i] / (true_model$alpha[i] + true_model$beta[i])
  cat(sprintf("%-12d %-10.3f %-10.3f %-10.3f %-10.3f\n",
              true_model$component[i],
              true_model$alpha[i],
              true_model$beta[i],
              true_model$pi[i],
              mean_meth))
}

cat("\n")

# Read fitted parameters
fitted_model_file <- file.path(opt_dir, "optModel.csv")
if (file.exists(fitted_model_file)) {
  fitted_model <- read.csv(fitted_model_file)
  
  # Sort by mean methylation for alignment
  fitted_model$mean_meth <- fitted_model$alpha / (fitted_model$alpha + fitted_model$beta)
  fitted_model <- fitted_model[order(fitted_model$mean_meth), ]
  
  cat("Fitted Model Parameters:\n")
  cat(sprintf("%-12s %-10s %-10s %-10s %-10s\n",
              "Component", "Alpha", "Beta", "Pi", "Mean Meth"))
  cat(paste(rep("-", 60), collapse = ""), "\n")
  for (i in 1:nrow(fitted_model)) {
    cat(sprintf("%-12d %-10.3f %-10.3f %-10.3f %-10.3f\n",
                fitted_model$component[i],
                fitted_model$alpha[i],
                fitted_model$beta[i],
                fitted_model$pi[i],
                fitted_model$mean_meth[i]))
  }
  
  cat("\n")
  
  # Calculate errors
  cat("Parameter Recovery Errors:\n")
  cat(sprintf("%-12s %-12s %-12s %-12s\n",
              "Component", "Alpha", "Beta", "Pi"))
  cat(paste(rep("-", 55), collapse = ""), "\n")
  
  for (i in 1:nrow(true_model)) {
    alpha_err <- abs(fitted_model$alpha[i] - true_model$alpha[i]) / true_model$alpha[i] * 100
    beta_err <- abs(fitted_model$beta[i] - true_model$beta[i]) / true_model$beta[i] * 100
    pi_err <- abs(fitted_model$pi[i] - true_model$pi[i]) / true_model$pi[i] * 100
    
    cat(sprintf("%-12d %-12.1f%% %-12.1f%% %-12.1f%%\n",
                i, alpha_err, beta_err, pi_err))
  }
}

cat("\n")

# Read comparison summary if available
if (file.exists(comparison_csv)) {
  comparison <- read.csv(comparison_csv)
  
  cat("Overall Recovery Metrics (RMSE):\n")
  cat(sprintf("  Alpha RMSE:          %.4f\n", comparison$alpha_rmse))
  cat(sprintf("  Beta RMSE:           %.4f\n", comparison$beta_rmse))
  cat(sprintf("  Pi RMSE:             %.4f\n", comparison$pi_rmse))
  cat(sprintf("  Mean Meth RMSE:      %.4f\n", comparison$mean_meth_rmse))
  
  if ("alpha_cv" %in% names(comparison)) {
    cat("\nConvergence Stability (CV):\n")
    cat(sprintf("  Alpha CV:            %.6f\n", comparison$alpha_cv))
    cat(sprintf("  Beta CV:             %.6f\n", comparison$beta_cv))
    cat(sprintf("  Pi CV:               %.6f\n", comparison$pi_cv))
    cat(sprintf("  Convergence rate:    %.1f%%\n", comparison$convergence_rate * 100))
  }
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Workflow Complete\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Generated files:\n")
cat(sprintf("  %s/\n", output_dir))
cat(sprintf("    ├── synthetic_K%d.qs              (synthetic data)\n", K))
cat(sprintf("    ├── synthetic_K%d_trueModel.csv   (true parameters)\n", K))
cat(sprintf("    ├── %s/               (fitted model)\n", basename(opt_dir)))
cat("    └── parameter_comparison.csv     (comparison metrics)\n\n")

cat("Interpretation:\n")
cat("  - RMSE: Lower values indicate better parameter recovery\n")
cat("  - Relative error < 10%: Good recovery\n")
cat("  - CV < 0.01: Excellent convergence stability\n")
cat("  - Convergence rate should be 100%\n\n")

cat("Next steps:\n")
cat("  1. Try different K values to test model selection\n")
cat("  2. Test with larger datasets: increase -n parameter\n")
cat("  3. Test sparse coverage: add --sparse flag to generator\n")
cat("  4. Apply to your own methylation data\n\n")

cat("To clean up this test:\n")
cat(sprintf("  rm -rf %s/\n\n", output_dir))
