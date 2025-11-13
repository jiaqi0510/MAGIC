#!/usr/bin/env Rscript
#
# MAGIC: Methylation Analysis with Genomic Inferred Contexts
#
# Module: Model comparison and evaluation
# Description: Compares fitted mixture models against true models for validation.
#              Computes parameter recovery metrics (RMSE, relative errors), analyzes
#              convergence stability across multiple runs, and generates performance
#              visualization plots.
#
# Authors: Jiaqi Han, Michael Thompson, Matteo Pellegrini
# Contact: mjthompson69@gmail.com
# License: MIT
# Date: 2025

# No additional libraries needed

compare_models <- function(true_model_file, fitted_model_file) {
  true_model <- read.csv(true_model_file)
  fitted_model <- read.csv(fitted_model_file)
  
  K <- nrow(true_model)
  
  # Sort both models by mean methylation to align components
  true_model$mean_meth <- true_model$alpha / (true_model$alpha + true_model$beta)
  fitted_model$mean_meth <- fitted_model$alpha / (fitted_model$alpha + fitted_model$beta)
  
  true_model <- true_model[order(true_model$mean_meth), ]
  fitted_model <- fitted_model[order(fitted_model$mean_meth), ]
  
  # Calculate parameter differences
  alpha_diff <- fitted_model$alpha - true_model$alpha
  beta_diff <- fitted_model$beta - true_model$beta
  pi_diff <- fitted_model$pi - true_model$pi
  mean_meth_diff <- fitted_model$mean_meth - true_model$mean_meth
  
  # Calculate relative errors with bounds checking
  alpha_rel_err <- ifelse(true_model$alpha > 1e-10, 
                          alpha_diff / true_model$alpha, 
                          NA_real_)
  beta_rel_err <- ifelse(true_model$beta > 1e-10, 
                         beta_diff / true_model$beta, 
                         NA_real_)
  pi_rel_err <- ifelse(true_model$pi > 1e-10, 
                       pi_diff / true_model$pi, 
                       NA_real_)
  mean_meth_rel_err <- ifelse(true_model$mean_meth > 1e-10, 
                              mean_meth_diff / true_model$mean_meth, 
                              NA_real_)
  
  if (any(is.na(alpha_rel_err))) {
    warning(sprintf("[VALIDATION] %d component(s) have near-zero true alpha - relative error set to NA", 
                   sum(is.na(alpha_rel_err))))
  }
  if (any(is.na(beta_rel_err))) {
    warning(sprintf("[VALIDATION] %d component(s) have near-zero true beta - relative error set to NA",
                   sum(is.na(beta_rel_err))))
  }
  if (any(is.na(pi_rel_err))) {
    warning(sprintf("[VALIDATION] %d component(s) have near-zero true pi - relative error set to NA",
                   sum(is.na(pi_rel_err))))
  }
  if (any(is.na(mean_meth_rel_err))) {
    warning(sprintf("[VALIDATION] %d component(s) have near-zero true mean methylation - relative error set to NA",
                   sum(is.na(mean_meth_rel_err))))
  }
  
  comparison <- data.frame(
    component = 1:K,
    true_alpha = true_model$alpha,
    fitted_alpha = fitted_model$alpha,
    alpha_diff = alpha_diff,
    alpha_rel_err = alpha_rel_err,
    true_beta = true_model$beta,
    fitted_beta = fitted_model$beta,
    beta_diff = beta_diff,
    beta_rel_err = beta_rel_err,
    true_pi = true_model$pi,
    fitted_pi = fitted_model$pi,
    pi_diff = pi_diff,
    pi_rel_err = pi_rel_err,
    true_mean_meth = true_model$mean_meth,
    fitted_mean_meth = fitted_model$mean_meth,
    mean_meth_diff = mean_meth_diff,
    mean_meth_rel_err = mean_meth_rel_err
  )
  
  # Print summary
  cat("Parameter Recovery:\n")
  cat(sprintf("  Alpha RMSE: %.4f\n", sqrt(mean(alpha_diff^2))))
  cat(sprintf("  Beta RMSE:  %.4f\n", sqrt(mean(beta_diff^2))))
  cat(sprintf("  Pi RMSE:    %.4f\n", sqrt(mean(pi_diff^2))))
  cat(sprintf("  Mean methylation RMSE: %.4f\n\n", sqrt(mean(mean_meth_diff^2))))
  
  cat("Mean Absolute Relative Errors:\n")
  cat(sprintf("  Alpha: %.1f%%\n", 100 * mean(abs(alpha_rel_err), na.rm = TRUE)))
  cat(sprintf("  Beta:  %.1f%%\n", 100 * mean(abs(beta_rel_err), na.rm = TRUE)))
  cat(sprintf("  Pi:    %.1f%%\n", 100 * mean(abs(pi_rel_err), na.rm = TRUE)))
  cat(sprintf("  Mean methylation: %.1f%%\n\n", 100 * mean(abs(mean_meth_rel_err), na.rm = TRUE)))
  
  cat("Component-wise comparison:\n")
  for (i in 1:K) {
    cat(sprintf("Component %d:\n", i))
    cat(sprintf("  True:   α=%.3f, β=%.3f, π=%.3f, mean=%.3f\n", 
                true_model$alpha[i], true_model$beta[i], true_model$pi[i], true_model$mean_meth[i]))
    cat(sprintf("  Fitted: α=%.3f, β=%.3f, π=%.3f, mean=%.3f\n", 
                fitted_model$alpha[i], fitted_model$beta[i], fitted_model$pi[i], fitted_model$mean_meth[i]))
    cat(sprintf("  Rel err: α=%.1f%%, β=%.1f%%, π=%.1f%%, mean=%.1f%%\n\n",
                100*alpha_rel_err[i], 100*beta_rel_err[i], 100*pi_rel_err[i], 100*mean_meth_rel_err[i]))
  }
  
  return(comparison)
}

analyze_convergence_stability <- function(all_param_file, true_model_file) {
  all_params <- read.csv(all_param_file)
  true_model <- read.csv(true_model_file)
  
  n_components <- nrow(true_model)
  n_runs <- nrow(all_params) / n_components
  
  # Align components by mean methylation for each run
  aligned_params <- data.frame()
  for (run in unique(all_params$run_id)) {
    run_data <- all_params[all_params$run_id == run, ]
    run_data <- run_data[order(run_data$mean_meth), ]
    run_data$aligned_component <- 1:nrow(run_data)
    aligned_params <- rbind(aligned_params, run_data)
  }
  
  # Calculate CV for each parameter within each component
  component_cvs <- data.frame()
  for (comp in 1:n_components) {
    comp_data <- aligned_params[aligned_params$aligned_component == comp, ]
    
    alpha_cv <- sd(comp_data$alpha) / mean(comp_data$alpha)
    beta_cv <- sd(comp_data$beta) / mean(comp_data$beta)
    pi_cv <- sd(comp_data$pi) / mean(comp_data$pi)
    mean_meth_cv <- sd(comp_data$mean_meth) / mean(comp_data$mean_meth)
    
    component_cvs <- rbind(component_cvs, data.frame(
      component = comp, alpha_cv = alpha_cv, beta_cv = beta_cv, 
      pi_cv = pi_cv, mean_meth_cv = mean_meth_cv
    ))
  }
  
  # Overall convergence statistics
  return(list(
    n_runs = n_runs,
    alpha_cv = mean(component_cvs$alpha_cv),
    beta_cv = mean(component_cvs$beta_cv),
    pi_cv = mean(component_cvs$pi_cv),
    mean_meth_cv = mean(component_cvs$mean_meth_cv),
    convergence_rate = length(unique(all_params$run_id)) / n_runs
  ))
}

# Command line interface
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  cat("Compare true and fitted beta-binomial mixture models\n\n")
  cat("Usage: Rscript model_comparison.R -t <true_dir> -f <fitted_dir> [options]\n\n")
  cat("Required:\n")
  cat("  -t, --true-dir <dir>     Directory containing true model files (*_trueModel.csv)\n")
  cat("  -f, --fitted-dir <dir>   Directory containing optimization subdirs (*_K<number>_*)\n\n")
  cat("Options:\n")
  cat("  -o, --output <file>      Save summary comparison to CSV\n")
  cat("  -h, --help               Show this help\n\n")
  quit(status = 0)
}

# Parse flags
true_dir <- NULL
fitted_dir <- NULL
output_file <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "-t" || args[i] == "--true-dir") {
    true_dir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "-f" || args[i] == "--fitted-dir") {
    fitted_dir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "-o" || args[i] == "--output") {
    output_file <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

# Find true model files and optimization directories
true_files <- list.files(true_dir, pattern = "_K[0-9]+.*_trueModel\\.csv$", full.names = TRUE)
opt_dirs <- list.dirs(fitted_dir, recursive = FALSE)
opt_dirs <- opt_dirs[grepl("_K[0-9]+_", opt_dirs)]

cat("Found", length(true_files), "true model files and", length(opt_dirs), "optimization directories\n\n")

# Extract K values
extract_k_from_file <- function(filename) {
  match <- regexpr("_K([0-9]+)", basename(filename))
  as.integer(sub("_K([0-9]+)", "\\1", regmatches(basename(filename), match)))
}

extract_k_from_dir <- function(dirname) {
  match <- regexpr("_K([0-9]+)_", basename(dirname))
  as.integer(sub("_K([0-9]+)_", "\\1", regmatches(basename(dirname), match)))
}

true_k_values <- sapply(true_files, extract_k_from_file)
names(true_k_values) <- true_files
opt_k_values <- sapply(opt_dirs, extract_k_from_dir)
names(opt_k_values) <- opt_dirs

# Find matching pairs and compare
all_results <- list()

for (k in sort(unique(c(true_k_values, opt_k_values)))) {
  true_file_for_k <- names(true_k_values)[true_k_values == k]
  opt_dir_for_k <- names(opt_k_values)[opt_k_values == k]
  
  if (length(true_file_for_k) == 0 || length(opt_dir_for_k) == 0) next
  
  true_file <- true_file_for_k[1]
  opt_dir <- opt_dir_for_k[1]
  opt_model_file <- file.path(opt_dir, "optModel.csv")
  
  if (!file.exists(opt_model_file)) next
  
  cat("Evaluating K =", k, ":\n")
  cat(paste(rep("=", 40), collapse=""), "\n")
  
  result <- compare_models(true_file, opt_model_file)
  
  # Analyze convergence stability
  all_param_file <- file.path(opt_dir, "allParamEst.csv")
  convergence_stats <- NULL
  if (file.exists(all_param_file)) {
    convergence_stats <- analyze_convergence_stability(all_param_file, true_file)
  }
  
  # Store summary metrics
  alpha_rmse <- sqrt(mean((result$fitted_alpha - result$true_alpha)^2))
  beta_rmse <- sqrt(mean((result$fitted_beta - result$true_beta)^2))
  pi_rmse <- sqrt(mean((result$fitted_pi - result$true_pi)^2))
  mean_meth_rmse <- sqrt(mean((result$fitted_mean_meth - result$true_mean_meth)^2))
  
  result_row <- data.frame(
    K = k,
    true_file = basename(true_file),
    opt_dir = basename(opt_dir),
    alpha_rmse = alpha_rmse,
    beta_rmse = beta_rmse,
    pi_rmse = pi_rmse,
    mean_meth_rmse = mean_meth_rmse
  )
  
  if (!is.null(convergence_stats)) {
    result_row$n_runs <- convergence_stats$n_runs
    result_row$alpha_cv <- convergence_stats$alpha_cv
    result_row$beta_cv <- convergence_stats$beta_cv
    result_row$pi_cv <- convergence_stats$pi_cv
    result_row$mean_meth_cv <- convergence_stats$mean_meth_cv
    result_row$convergence_rate <- convergence_stats$convergence_rate
  }
  
  all_results[[paste0("K", k)]] <- result_row
  cat("\n")
}

if (length(all_results) > 0) {
  summary_df <- do.call(rbind, all_results)
  
  cat("Summary of Model Comparisons:\n")
  cat("============================\n")
  print(summary_df, row.names = FALSE)
  
  if (!is.null(output_file)) {
    write.csv(summary_df, output_file, row.names = FALSE)
    cat("\nSummary saved to:", output_file, "\n")
  }
}