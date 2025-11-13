#!/usr/bin/env Rscript
#
# MAGIC: Methylation Analysis with Genomic Inferred Contexts
#
# Example: Basic mixture model optimization
# Description: Demonstrates how to optimize a beta-binomial mixture model using
#              MAGIC. Covers basic parameter estimation, model selection across
#              different K values, convergence diagnostics, and output interpretation.
#
library(bsseq)
library(qs)

# ==============================================================================
# Setup
# ==============================================================================

# This example assumes you have a BSseq object saved as a .qs file
# If you need to create one, see the bsseq package documentation

# Input data (modify path to your data)
input_file <- "data/BSseq_example.qs"

# Check if file exists
if (!file.exists(input_file)) {
  cat("Example data file not found:", input_file, "\n")
  cat("Please prepare a BSseq object and save it using:\n")
  cat("  qs::qsave(bsseq_object, 'data/BSseq_example.qs')\n\n")
  cat("The BSseq object should contain:\n")
  cat("  - Methylation counts: getCoverage(bs, type='M')\n")
  cat("  - Total coverage: getCoverage(bs, type='Cov')\n")
  cat("  - Sample metadata: colData(bs)\n")
  quit(save = "no", status = 1)
}

# ==============================================================================
# Example 1: Single Model Optimization
# ==============================================================================

cat("Example 1: Optimizing K=3 mixture model\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Basic optimization with K=3 components
# This will create an output directory with results

system2("Rscript", c(
  "R/magicFit.R",
  "--input", input_file,
  "--k_components", "3",
  "--max_iters", "500",
  "--tol", "1e-6",
  "--threads", "4",
  "--n_conv_runs", "3"  # Run 3 times with different initializations
))

cat("\n")
cat("Results saved to: [basename]_K3_S10_e06_it500/\n")
cat("Key output files:\n")
cat("  - optModel.csv: Fitted parameters (alpha, beta, pi)\n")
cat("  - methComps.csv: Interpretable parameters (mixing weights, mean, dispersion)\n")
cat("  - runSum.csv: BIC, AIC, convergence info\n")
cat("\n")

# ==============================================================================
# Example 2: Model Selection (Compare K=2,3,4,5)
# ==============================================================================

cat("\nExample 2: Model selection across K values\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Try different numbers of components
k_values <- 2:5

cat("Fitting models with K =", paste(k_values, collapse = ", "), "\n")
cat("This may take several minutes...\n\n")

for (k in k_values) {
  cat(sprintf("Fitting K=%d...\n", k))
  
  system2("Rscript", c(
    "R/magicFit.R",
    "--input", input_file,
    "--k_components", as.character(k),
    "--max_iters", "500",
    "--tol", "1e-6",
    "--threads", "4",
    "--n_conv_runs", "5",
    "--quiet"  # Suppress detailed output
  ))
}

# Compare BIC values
cat("\nComparing models by BIC:\n")
cat("(Lower BIC indicates better model)\n\n")

# Find all runSum.csv files
result_dirs <- list.dirs(".", recursive = FALSE)
result_dirs <- result_dirs[grepl("_K[0-9]+_", result_dirs)]

if (length(result_dirs) > 0) {
  bic_results <- data.frame()
  
  for (dir in result_dirs) {
    runsum_file <- file.path(dir, "runSum.csv")
    if (file.exists(runsum_file)) {
      runsum <- read.csv(runsum_file)
      bic_results <- rbind(bic_results, data.frame(
        K = runsum$k,
        BIC = runsum$bic,
        AIC = runsum$aic,
        LogLik = runsum$train_ll,
        Converged = runsum$conv,
        Directory = basename(dir)
      ))
    }
  }
  
  if (nrow(bic_results) > 0) {
    bic_results <- bic_results[order(bic_results$K), ]
    print(bic_results, row.names = FALSE)
    
    best_k <- bic_results$K[which.min(bic_results$BIC)]
    cat(sprintf("\nBest model: K=%d (lowest BIC)\n", best_k))
  }
}

# ==============================================================================
# Example 3: Chromosome-Specific Analysis
# ==============================================================================

cat("\n\nExample 3: Chromosome-specific analysis\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Analyze a single chromosome (useful for testing or memory constraints)
chr_to_analyze <- "chr21"

cat(sprintf("Analyzing only %s\n", chr_to_analyze))

system2("Rscript", c(
  "R/magicFit.R",
  "--input", input_file,
  "--k_components", "3",
  "--chr", chr_to_analyze,
  "--threads", "4",
  "--n_conv_runs", "3"
))

cat(sprintf("\nResults for %s saved with '_chr21' suffix\n", chr_to_analyze))

# ==============================================================================
# Example 4: With Holdout Validation
# ==============================================================================

cat("\n\nExample 4: Optimization with holdout validation\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Use holdout set to assess overfitting
# 20% of sites reserved for testing

system2("Rscript", c(
  "R/magicFit.R",
  "--input", input_file,
  "--k_components", "3",
  "--threads", "4",
  "--holdout",
  "--holdout_frac", "0.2",
  "--n_conv_runs", "5"
))

cat("\nCheck runSum.csv for:\n")
cat("  - train_ll: Log-likelihood on training set\n")
cat("  - test_ll: Log-likelihood on holdout set\n")
cat("  - overfitting: Ratio test_ll / train_ll\n")
cat("\nIf overfitting < 0.95, model may be overfitting\n")

# ==============================================================================
# Example 5: Reading and Interpreting Results
# ==============================================================================

cat("\n\nExample 5: Reading and interpreting results\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Find the most recent K=3 result directory
result_dirs <- list.dirs(".", recursive = FALSE)
k3_dirs <- result_dirs[grepl("_K3_", result_dirs)]

if (length(k3_dirs) > 0) {
  # Use most recent
  result_dir <- k3_dirs[length(k3_dirs)]
  
  cat("Reading results from:", result_dir, "\n\n")
  
  # Read interpretable parameters
  methcomps_file <- file.path(result_dir, "methComps.csv")
  if (file.exists(methcomps_file)) {
    methcomps <- read.csv(methcomps_file)
    
    cat("Mixture Components:\n")
    cat(sprintf("%-10s %-12s %-12s %-12s\n", 
                "Component", "Mix Weight", "Mean Meth", "Dispersion"))
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    for (i in 1:nrow(methcomps)) {
      cat(sprintf("%-10d %-12.3f %-12.3f %-12.5f\n",
                  methcomps$component[i],
                  methcomps$mixWeight[i],
                  methcomps$meanMeth[i],
                  methcomps$methDisp[i]))
    }
    
    cat("\nInterpretation:\n")
    cat("  - Mix Weight: Proportion of sites in this component\n")
    cat("  - Mean Meth: Average methylation level (0-1)\n")
    cat("  - Dispersion: Variability (lower = more consistent)\n")
    cat("\nTypical patterns:\n")
    cat("  - Component 1: Low methylation (unmethylated)\n")
    cat("  - Component 2: Intermediate methylation (partially methylated)\n")
    cat("  - Component 3: High methylation (highly methylated)\n")
  }
  
  # Read convergence info
  cat("\n")
  runsum_file <- file.path(result_dir, "runSum.csv")
  if (file.exists(runsum_file)) {
    runsum <- read.csv(runsum_file)
    
    cat("Convergence Summary:\n")
    cat(sprintf("  Converged: %s\n", runsum$conv))
    cat(sprintf("  Iterations: %d\n", runsum$nIters))
    cat(sprintf("  BIC: %.1f\n", runsum$bic))
    cat(sprintf("  AIC: %.1f\n", runsum$aic))
    cat(sprintf("  Runs: %d converged / %d total\n", runsum$n_conv, runsum$n_runs))
    cat(sprintf("  Runtime: %.1f seconds\n", runsum$runtime))
  }
  
  # Check convergence stability
  allparams_file <- file.path(result_dir, "allParamEst.csv")
  if (file.exists(allparams_file)) {
    allparams <- read.csv(allparams_file)
    
    # Calculate CV for each component
    cat("\nConvergence Stability (across runs):\n")
    for (comp in unique(allparams$component)) {
      comp_data <- allparams[allparams$component == comp, ]
      mean_meth_cv <- sd(comp_data$mean_meth) / mean(comp_data$mean_meth)
      
      cat(sprintf("  Component %d: CV = %.4f", comp, mean_meth_cv))
      if (mean_meth_cv < 0.01) {
        cat(" (excellent)")
      } else if (mean_meth_cv < 0.05) {
        cat(" (good)")
      } else {
        cat(" (check convergence)")
      }
      cat("\n")
    }
    cat("\nLower CV indicates more consistent convergence\n")
  }
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Example Complete\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("What you learned:\n")
cat("  1. Basic model optimization\n")
cat("  2. Model selection using BIC\n")
cat("  3. Chromosome-specific analysis\n")
cat("  4. Holdout validation\n")
cat("  5. Interpreting results\n\n")

cat("Next steps:\n")
cat("  1. Apply to your own data\n")
cat("  2. See example_testing.R for differential methylation\n")
cat("  3. See docs/usage.md for complete parameter reference\n\n")

cat("Output directories contain:\n")
cat("  - optModel.csv: Parameters for downstream testing\n")
cat("  - methComps.csv: Human-readable component summaries\n")
cat("  - runSum.csv: Model fit statistics\n")
cat("  - allParamEst.csv: All convergence runs (for diagnostics)\n\n")
