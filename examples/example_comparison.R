#!/usr/bin/env Rscript
#
# MAGIC: Methylation Analysis with Genomic Inferred Contexts
#
# Example: Model comparison for simulations
# Description: Demonstrates how to evaluate parameter recovery using the
#              model_comparison.R utility. Useful for validating the method on
#              simulated data, assessing parameter recovery accuracy, comparing
#              performance across different K values, and publishing method
#              validation results.
#
library(ggplot2)

# ==============================================================================
# Setup
# ==============================================================================

# This example assumes you have:
# 1. Simulated data with known true parameters
# 2. Fitted models from running magicFit.R on simulated data

# Directory structure should be:
# simulations/
#   ├── sim1_K3_trueModel.csv
#   ├── sim2_K4_trueModel.csv
#   └── ...
# results/
#   ├── sim1_K3_S10_e06_it500/
#   │   └── optModel.csv
#   ├── sim2_K4_S10_e06_it500/
#   │   └── optModel.csv
#   └── ...

true_model_dir <- "simulations"
fitted_model_dir <- "results"

# Check if directories exist
if (!dir.exists(true_model_dir)) {
  cat("Error: True model directory not found:", true_model_dir, "\n")
  cat("\nThis example requires simulated data with known parameters.\n")
  cat("True model files should be named: *_K[n]_trueModel.csv\n")
  cat("Format: CSV with columns 'alpha', 'beta', 'pi'\n\n")
  quit(save = "no", status = 1)
}

if (!dir.exists(fitted_model_dir)) {
  cat("Error: Fitted model directory not found:", fitted_model_dir, "\n")
  cat("Please run magicFit.R on your simulated data first\n")
  quit(save = "no", status = 1)
}

# ==============================================================================
# Example 1: Basic Model Comparison
# ==============================================================================

cat("Example 1: Comparing fitted vs true parameters\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Run comparison script
system2("Rscript", c(
  "scripts/model_comparison.R",
  "--true-dir", true_model_dir,
  "--fitted-dir", fitted_model_dir,
  "--output", "comparison_summary.csv",
  "--plots", "comparison_plots.pdf"
))

cat("\nOutputs created:\n")
cat("  - comparison_summary.csv: RMSE metrics for each model\n")
cat("  - comparison_plots.pdf: Visualization of parameter recovery\n\n")

# ==============================================================================
# Example 2: Reading Comparison Results
# ==============================================================================

cat("Example 2: Interpreting comparison metrics\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

if (file.exists("comparison_summary.csv")) {
  comparison <- read.csv("comparison_summary.csv")
  
  cat("Comparison summary contains:\n")
  cat("  K: Number of mixture components\n")
  cat("  alpha_rmse: Root mean squared error for alpha parameters\n")
  cat("  beta_rmse: RMSE for beta parameters\n")
  cat("  pi_rmse: RMSE for mixing proportions\n")
  cat("  mean_meth_rmse: RMSE for mean methylation\n")
  cat("  alpha_cv: Coefficient of variation across convergence runs\n")
  cat("  convergence_rate: Fraction of runs that converged\n\n")
  
  cat("Results by K:\n")
  print(comparison, row.names = FALSE)
  
  cat("\n\nInterpretation:\n")
  cat("  Lower RMSE = better parameter recovery\n")
  cat("  Lower CV = more stable convergence\n")
  cat("  convergence_rate should be close to 1.0\n\n")
  
  # Identify best K
  if ("mean_meth_rmse" %in% names(comparison)) {
    best_k <- comparison$K[which.min(comparison$mean_meth_rmse)]
    cat(sprintf("Best parameter recovery: K=%d (lowest mean methylation RMSE)\n", best_k))
  }
}

# ==============================================================================
# Example 3: Creating Custom Visualizations
# ==============================================================================

cat("\n\nExample 3: Custom visualization of results\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

if (file.exists("comparison_summary.csv")) {
  comparison <- read.csv("comparison_summary.csv")
  
  # RMSE by parameter type
  cat("Creating parameter recovery plot...\n")
  
  # Reshape data for plotting
  plot_data <- data.frame(
    K = rep(comparison$K, 3),
    RMSE = c(comparison$alpha_rmse, 
             comparison$beta_rmse, 
             comparison$pi_rmse),
    Parameter = rep(c("Alpha", "Beta", "Pi"), 
                    each = nrow(comparison))
  )
  
  p1 <- ggplot(plot_data, aes(x = K, y = RMSE, color = Parameter)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "Parameter Recovery by K",
         x = "Number of Components (K)",
         y = "Root Mean Squared Error",
         color = "Parameter") +
    theme(legend.position = "bottom")
  
  ggsave("parameter_recovery.png", p1, width = 8, height = 6)
  cat("  Saved: parameter_recovery.png\n")
  
  # Convergence stability
  if ("alpha_cv" %in% names(comparison)) {
    cat("Creating convergence stability plot...\n")
    
    stability_data <- data.frame(
      K = rep(comparison$K, 3),
      CV = c(comparison$alpha_cv,
             comparison$beta_cv,
             comparison$pi_cv),
      Parameter = rep(c("Alpha", "Beta", "Pi"),
                      each = nrow(comparison))
    )
    
    p2 <- ggplot(stability_data, aes(x = K, y = CV, color = Parameter)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      theme_minimal() +
      labs(title = "Convergence Stability by K",
           x = "Number of Components (K)",
           y = "Coefficient of Variation",
           subtitle = "Lower values indicate more consistent convergence",
           color = "Parameter") +
      theme(legend.position = "bottom")
    
    ggsave("convergence_stability.png", p2, width = 8, height = 6)
    cat("  Saved: convergence_stability.png\n")
  }
}

# ==============================================================================
# Example 4: Component-Wise Comparison
# ==============================================================================

cat("\n\nExample 4: Detailed component-wise comparison\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# For a specific K value, compare each component
k_to_examine <- 3

true_file <- file.path(true_model_dir, 
                       sprintf("sim1_K%d_trueModel.csv", k_to_examine))
fitted_dir <- list.dirs(fitted_model_dir, recursive = FALSE)
fitted_dir <- fitted_dir[grepl(sprintf("_K%d_", k_to_examine), fitted_dir)][1]
fitted_file <- file.path(fitted_dir, "optModel.csv")

if (file.exists(true_file) && file.exists(fitted_file)) {
  true_model <- read.csv(true_file)
  fitted_model <- read.csv(fitted_file)
  
  cat(sprintf("Comparing K=%d model components:\n\n", k_to_examine))
  
  # Sort by mean methylation for alignment
  true_model$mean_meth <- true_model$alpha / (true_model$alpha + true_model$beta)
  fitted_model$mean_meth <- fitted_model$alpha / (fitted_model$alpha + fitted_model$beta)
  
  true_model <- true_model[order(true_model$mean_meth), ]
  fitted_model <- fitted_model[order(fitted_model$mean_meth), ]
  
  cat(sprintf("%-10s %-12s %-12s %-12s %-12s\n",
              "Component", "Parameter", "True", "Fitted", "Rel Error"))
  cat(paste(rep("-", 65), collapse = ""), "\n")
  
  for (comp in 1:nrow(true_model)) {
    # Alpha
    alpha_err <- (fitted_model$alpha[comp] - true_model$alpha[comp]) / true_model$alpha[comp]
    cat(sprintf("%-10d %-12s %-12.3f %-12.3f %-12.1f%%\n",
                comp, "Alpha",
                true_model$alpha[comp],
                fitted_model$alpha[comp],
                100 * alpha_err))
    
    # Beta
    beta_err <- (fitted_model$beta[comp] - true_model$beta[comp]) / true_model$beta[comp]
    cat(sprintf("%-10s %-12s %-12.3f %-12.3f %-12.1f%%\n",
                "", "Beta",
                true_model$beta[comp],
                fitted_model$beta[comp],
                100 * beta_err))
    
    # Pi
    pi_err <- (fitted_model$pi[comp] - true_model$pi[comp]) / true_model$pi[comp]
    cat(sprintf("%-10s %-12s %-12.3f %-12.3f %-12.1f%%\n",
                "", "Pi",
                true_model$pi[comp],
                fitted_model$pi[comp],
                100 * pi_err))
    
    # Mean methylation
    mean_err <- (fitted_model$mean_meth[comp] - true_model$mean_meth[comp]) / 
                true_model$mean_meth[comp]
    cat(sprintf("%-10s %-12s %-12.3f %-12.3f %-12.1f%%\n\n",
                "", "Mean Meth",
                true_model$mean_meth[comp],
                fitted_model$mean_meth[comp],
                100 * mean_err))
  }
  
  cat("Interpretation:\n")
  cat("  Relative error < 5%: Excellent recovery\n")
  cat("  Relative error < 10%: Good recovery\n")
  cat("  Relative error > 20%: Check convergence or model fit\n")
  
} else {
  cat(sprintf("Could not find files for K=%d comparison\n", k_to_examine))
  cat("  True file:", true_file, "\n")
  cat("  Fitted file:", fitted_file, "\n")
}

# ==============================================================================
# Example 5: Multiple Simulation Replicates
# ==============================================================================

cat("\n\nExample 5: Comparing across simulation replicates\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# If you have multiple simulations with the same K
k_value <- 3
sim_pattern <- sprintf("*_K%d_trueModel.csv", k_value)

true_files <- list.files(true_model_dir, pattern = sprintf("K%d_trueModel\\.csv$", k_value), 
                         full.names = TRUE)

if (length(true_files) > 1) {
  cat(sprintf("Found %d simulation replicates for K=%d\n\n", length(true_files), k_value))
  
  # Collect RMSE across replicates
  rmse_results <- data.frame()
  
  for (true_file in true_files) {
    # Find corresponding fitted model
    base_name <- sub("_trueModel\\.csv$", "", basename(true_file))
    fitted_dirs <- list.dirs(fitted_model_dir, recursive = FALSE)
    fitted_dir <- fitted_dirs[grepl(base_name, fitted_dirs)][1]
    
    if (!is.na(fitted_dir) && dir.exists(fitted_dir)) {
      fitted_file <- file.path(fitted_dir, "optModel.csv")
      
      if (file.exists(fitted_file)) {
        true_model <- read.csv(true_file)
        fitted_model <- read.csv(fitted_file)
        
        # Calculate RMSE
        true_model$mean_meth <- true_model$alpha / (true_model$alpha + true_model$beta)
        fitted_model$mean_meth <- fitted_model$alpha / (fitted_model$alpha + fitted_model$beta)
        
        true_model <- true_model[order(true_model$mean_meth), ]
        fitted_model <- fitted_model[order(fitted_model$mean_meth), ]
        
        mean_meth_rmse <- sqrt(mean((fitted_model$mean_meth - true_model$mean_meth)^2))
        
        rmse_results <- rbind(rmse_results, data.frame(
          simulation = basename(true_file),
          rmse = mean_meth_rmse
        ))
      }
    }
  }
  
  if (nrow(rmse_results) > 0) {
    cat("Mean methylation RMSE by simulation:\n")
    print(rmse_results, row.names = FALSE)
    
    cat(sprintf("\nSummary across %d replicates:\n", nrow(rmse_results)))
    cat(sprintf("  Mean RMSE: %.4f\n", mean(rmse_results$rmse)))
    cat(sprintf("  SD RMSE: %.4f\n", sd(rmse_results$rmse)))
    cat(sprintf("  Min RMSE: %.4f\n", min(rmse_results$rmse)))
    cat(sprintf("  Max RMSE: %.4f\n", max(rmse_results$rmse)))
  }
} else {
  cat(sprintf("Only found %d simulation for K=%d\n", length(true_files), k_value))
  cat("This example requires multiple simulation replicates\n")
}

# ==============================================================================
# Example 6: Creating Publication-Quality Figures
# ==============================================================================

cat("\n\nExample 6: Publication-quality comparison figures\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

if (file.exists("comparison_summary.csv")) {
  comparison <- read.csv("comparison_summary.csv")
  
  cat("Creating publication figures...\n")
  
  # Combined RMSE plot
  combined_data <- data.frame(
    K = comparison$K,
    Mean_Methylation = comparison$mean_meth_rmse,
    Alpha = comparison$alpha_rmse,
    Beta = comparison$beta_rmse,
    Pi = comparison$pi_rmse
  )
  
  # Normalize to make comparable
  combined_data$Mean_Methylation_norm <- combined_data$Mean_Methylation * 1000
  
  plot_data <- data.frame(
    K = rep(comparison$K, 4),
    RMSE = c(combined_data$Alpha,
             combined_data$Beta,
             combined_data$Pi,
             combined_data$Mean_Methylation_norm),
    Parameter = factor(rep(c("α", "β", "π", "Mean Meth (×1000)"),
                           each = nrow(comparison)),
                       levels = c("α", "β", "π", "Mean Meth (×1000)"))
  )
  
  p <- ggplot(plot_data, aes(x = K, y = RMSE, color = Parameter, shape = Parameter)) +
    geom_line(size = 1.2) +
    geom_point(size = 4) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")) +
    scale_shape_manual(values = c(16, 17, 15, 18)) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Parameter Recovery Performance",
         x = "Number of Mixture Components (K)",
         y = "Root Mean Squared Error",
         color = "Parameter",
         shape = "Parameter")
  
  ggsave("publication_figure.png", p, width = 10, height = 7, dpi = 300)
  ggsave("publication_figure.pdf", p, width = 10, height = 7)
  
  cat("  Saved: publication_figure.png (300 DPI)\n")
  cat("  Saved: publication_figure.pdf (vector)\n")
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Example Complete\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("What you learned:\n")
cat("  1. Running model comparison utility\n")
cat("  2. Interpreting RMSE and CV metrics\n")
cat("  3. Component-wise parameter evaluation\n")
cat("  4. Analyzing multiple simulation replicates\n")
cat("  5. Creating publication-quality figures\n\n")

cat("Key metrics:\n")
cat("  - RMSE: Lower = better parameter recovery\n")
cat("  - CV: Lower = more stable convergence\n")
cat("  - Convergence rate: Should be close to 1.0\n")
cat("  - Relative error < 10%: Good recovery\n\n")

cat("Use cases:\n")
cat("  - Method validation for manuscripts\n")
cat("  - Comparing optimization strategies\n")
cat("  - Assessing identifiability\n")
cat("  - Choosing optimal K\n\n")

cat("Next steps:\n")
cat("  1. Run on your simulated datasets\n")
cat("  2. Report RMSE metrics in methods validation\n")
cat("  3. Use plots in supplementary materials\n")
cat("  4. See docs/usage.md for more options\n\n")
