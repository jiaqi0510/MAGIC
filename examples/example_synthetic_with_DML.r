#!/usr/bin/env Rscript
#
# MAGIC: Methylation Analysis with Genomic Inferred Contexts
#
# Example: Complete workflow with synthetic data (including true DML)
#
# Description: Demonstrates a full MAGIC workflow using synthetic WGBS data
#              with controlled differential methylation (DML):
#
#              1. Generate synthetic methylation data with pre-defined
#                 effect sizes
#              2. Run MAGIC differential methylation analysis on
#                 simulated case vs control groups
#              3. Visualize methylation patterns and distributions 
#              4. Evaluate performance using ROC curves and other metrics
#                 by comparing estimated DMLs to ground truth
#
#              This example is useful for:
#              - Validating MAGIC’s DML-detection performance
#              - Testing installation with a full working example
#              - Benchmarking MAGIC against other DML detection packages


# Load required libraries
suppressPackageStartupMessages({
  library(bsseq)
  library(qs)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(pROC)
})

# ==============================================================================
# STEP 1: Generate Synthetic Data
# ==============================================================================

cat("=== STEP 1: Generating Synthetic Data ===\n")

# Generate synthetic dataset with:
# - 3 mixture components (K=3)
# - 100,000 CpG sites
# - 10 samples (5 Control + 5 Treatment)
# - 10% DML with 0.1 methylation difference
# - Mean coverage of 10×
# - No filtering (exact CpG count)

system_cmd <- paste(
  "Rscript scripts/synthetic_bb_generator_with_DML.r",
  "-K 3",
  "-n 100000",
  "-s 10",
  "--dml_fraction 0.1",
  "--effect_size 0.1",
  "-c 10",
  "-f FALSE",
  "-o synthetic_example/synthetic_bb_100k_s10_d01_cov10.qs"
)

system(system_cmd)

cat("\n✓ Synthetic data generated successfully!\n")
cat("  Output files created in: synthetic_example/\n")
cat("  - synthetic_bb_100k_s10_d01_cov10.qs (BSseq object)\n")
cat("  - synthetic_bb_100k_s10_d01_cov10_truth_DML.csv (ground truth)\n")
cat("  - synthetic_bb_100k_s10_d01_cov10_trueModel.csv (mixture parameters)\n\n")

# ==============================================================================
# STEP 2: Run MAGIC Analysis
# ==============================================================================

cat("=== STEP 2: Running MAGIC Differential Methylation Analysis ===\n")

# Step 2a: Fit mixture model on Control group
cat("Step 2a: Fitting mixture model on Control samples...\n")
fit_cmd <- paste(
  "Rscript R/magicFit.R",
  "--input synthetic_example/synthetic_bb_100k_s10_d01_cov10.qs",
  "--k_components 3",
  "--max_iters 500",
  "--tol 1e-6",
  "--threads 4",
  "--subset_trait group",
  "--subset_trait_val Control"
)

system(fit_cmd)

cat("✓ Mixture model fitting completed!\n\n")

# Step 2b: Test for differential methylation between groups
cat("Step 2b: Testing for differential methylation...\n")
test_cmd <- paste(
  "Rscript R/magicTest.R",
  "--input synthetic_example/synthetic_bb_100k_s10_d01_cov10.qs",
  "--group group",
  "--compare Control,Treatment",
  "--mixture_model synthetic_example/synthetic_bb_100k_s10_d01_cov10_K3_methOnly_trait_Control_S5_e06_it500/",
  "--cores 4"
)

system(test_cmd)

cat("✓ MAGIC analysis completed!\n")
cat("  Results saved to: synthetic_example/magic_results/\n\n")

# ==============================================================================
# STEP 3: Visualize Methylation Distributions by Component
# ==============================================================================

cat("=== STEP 3: Creating Methylation Distribution Plots ===\n")

# Function to plot methylation histogram faceted by mixture component
plot_site_group_hist <- function(qs_file, out_file) {
  # Read BSseq data
  bs_obj <- qread(qs_file)
  
  # Extract methylation matrix (site × sample)
  meth_mat <- getMeth(bs_obj, type = "raw")  
  pdata <- pData(bs_obj)
  
  # Extract true component assignments from metadata 
  if (!is.null(metadata(bs_obj)$true_cpg_components)) {
    comp_vec <- metadata(bs_obj)$true_cpg_components
  } else {
    comp_vec <- rep(1, nrow(bs_obj))
    warning("No true_cpg_components found in metadata; using component = 1 for all sites.")
  }
  
  # Calculate mean methylation level in each group for each site
  group_means <- sapply(levels(pdata$group), function(g) {
    rowMeans(meth_mat[, pdata$group == g, drop = FALSE], na.rm = TRUE)
  })
  
  group_means <- as.data.frame(group_means)
  group_means$site_id <- seq_len(nrow(group_means))
  group_means$component <- comp_vec
  
  # Convert to long format for plotting
  df <- group_means %>%
    pivot_longer(cols = -c(site_id, component),
                 names_to = "Group",
                 values_to = "Methylation")
  
  # Extract delta value from filename
  delta_str <- str_extract(qs_file, "d\\d+")
  delta_val <- switch(delta_str,
                      "d01" = 0.1,
                      "d02" = 0.2,
                      "d05" = 0.05,
                      NA)
  
  # Create histogram plot faceted by component
  p <- ggplot(df, aes(x = Methylation, fill = Group)) +
    geom_histogram(position = "identity", bins = 50, alpha = 0.6) +
    facet_wrap(~ component, ncol = 3, scales = "fixed") +
    coord_cartesian(ylim = c(0, 6000)) +
    theme_bw(base_size = 14) +
    labs(title = paste0("Methylation Histogram by Component (Δ = ", delta_val, ")"),
         x = "Group-mean Methylation",
         y = "Number of CpG Sites") +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "plain"),
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text = element_text(face = "plain", size = 12)
    )
  
  # Save plot
  ggsave(out_file, p, width = 12, height = 6, dpi = 300)
  
  cat("✓ Plot saved to:", out_file, "\n")
}

# Generate methylation distribution plot
plot_site_group_hist(
  qs_file = "synthetic_example/synthetic_bb_100k_s10_d01_cov10.qs",
  out_file = "synthetic_example/methylation_histogram_by_component.png"
)


# ==============================================================================
# STEP 4: Generate ROC Curve for Method Comparison
# ==============================================================================

cat("=== STEP 4: Creating ROC Curves for Method Comparison ===\n")

# Function to calculate ROC curves
get_roc_df <- function(score, label, model_name) {
  roc_obj <- roc(response = label, predictor = score, levels = c(0, 1))
  auc_val <- auc(roc_obj)
  data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    model = paste0(model_name, " (AUC=", round(auc_val, 3), ")")
  )
}

# Function to plot ROC curves comparing MAGIC methods and baseline
plot_roc_magic <- function(magic_file,
                           allsites,
                           truth_labels,
                           methdiff_df,
                           method = c("log", "BH", "quantile", "wald"),
                           out_file = NULL,
                           delta_threshold = 0.1) {
  
  method <- match.arg(method)
  roc_df_all <- data.frame()
  
  # Transform scores for ROC computation
  transform_score <- function(score_vec, method) {
    if (method == "log") {
      score <- -log10(score_vec)
      score[!is.finite(score)] <- 0
    } else if (method == "BH") {
      score <- 1 - score_vec
    } else if (method == "quantile") {
      score <- rank(score_vec, ties.method = "average") / length(score_vec)
    } else if (method == "wald") {
      score <- abs(score_vec)
      score[is.na(score)] <- 0
    }
    return(score)
  }
  
  # Create binary labels (1 = true DML, 0 = non-DML)
  label_binary <- ifelse(allsites %in% truth_labels, 1, 0)
  
  # === Load MAGIC results ===
  if (is.character(magic_file)) {
    magic <- fread(magic_file)
  } else {
    magic <- as.data.table(magic_file)
  }
  
  magic$id <- paste0(magic$chr, "_", magic$start)
  
  # === MAGIC w3 test ===
  if (method == "wald") {
    wald_vals <- magic$w3Wald[match(allsites, magic$id)]
    score <- transform_score(wald_vals, method)
    roc_df_all <- rbind(roc_df_all, get_roc_df(score, label_binary, "MAGIC w3 Wald"))
  } else {
    q_vals <- magic$w3Q[match(allsites, magic$id)]
    q_vals[is.na(q_vals)] <- 1
    score <- transform_score(q_vals, method)
    roc_df_all <- rbind(roc_df_all, get_roc_df(score, label_binary, "MAGIC w3"))
  }
  
  # === MAGIC Bayes Factor ===
  bf_vals <- magic$bfBF[match(allsites, magic$id)]
  bf_vals[is.na(bf_vals)] <- 0
  score_bf <- if (method == "quantile") {
    transform_score(bf_vals, method)
  } else if (method == "log") {
    bf_vals
  } else {
    bf_vals / max(bf_vals, na.rm = TRUE)
  }
  roc_df_all <- rbind(roc_df_all, get_roc_df(score_bf, label_binary, "MAGIC BF"))
  
  # === Baseline: Methylation delta (full curve) ===
  delta_score <- abs(methdiff_df$meth_diff)
  roc_df_all <- rbind(
    roc_df_all,
    get_roc_df(delta_score, label_binary, "Methylation Delta")
  )
  
  # === Baseline: Single threshold point ===
  df_delta <- methdiff_df
  df_delta$label <- label_binary
  pred_single <- ifelse(abs(df_delta$meth_diff) >= delta_threshold, 1, 0)
  tp <- sum(pred_single == 1 & df_delta$label == 1)
  fp <- sum(pred_single == 1 & df_delta$label == 0)
  fn <- sum(pred_single == 0 & df_delta$label == 1)
  tn <- sum(pred_single == 0 & df_delta$label == 0)
  tpr_single <- tp / (tp + fn)
  fpr_single <- fp / (fp + tn)
  roc_point_df <- data.frame(
    fpr = fpr_single,
    tpr = tpr_single,
    model = paste0("Delta ≥ ", delta_threshold)
  )
  
  # === Create ROC plot ===
  p_roc <- ggplot(roc_df_all, aes(x = fpr, y = tpr, color = model)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(data = roc_point_df, aes(x = fpr, y = tpr, color = model),
               size = 3, shape = 17) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_brewer(palette = "Set1") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = paste("ROC Curve: MAGIC vs Baseline (Delta = 0.1, Coverage = 10×)"),
         x = "False Positive Rate",
         y = "True Positive Rate") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain"),
      legend.position = c(0.7, 0.3),
      legend.background = element_rect(fill = "white", color = "gray80"),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(out_file)) {
    ggsave(out_file, p_roc, width = 8, height = 6, dpi = 300)
    cat("✓ ROC plot saved to:", out_file, "\n")
  }
  
  return(list(roc_df = roc_df_all, plot = p_roc))
}

# Load ground truth and prepare data
truth <- fread("synthetic_example/synthetic_bb_100k_s10_d01_cov10_truth_DML.csv")
truth$id <- paste0(truth$chr, "_", truth$pos)
truth_labels <- truth$id

# Load BSseq object and compute methylation differences
bs_data <- qread("synthetic_example/synthetic_bb_100k_s10_d01_cov10.qs")
group1_idx <- 1:5
group2_idx <- 6:10
meth_mat <- getMeth(bs_data, type = "raw")
mean_group1 <- rowMeans(meth_mat[, group1_idx], na.rm = TRUE)
mean_group2 <- rowMeans(meth_mat[, group2_idx], na.rm = TRUE)
methDiff <- mean_group2 - mean_group1
gr <- granges(bs_data)
allsites <- paste0(seqnames(gr), "_", start(gr))
methdiff_df <- data.frame(id = allsites, meth_diff = methDiff)

# Load MAGIC results
magic_results <- fread("synthetic_example/magic_results/magic_results.csv")

# Generate ROC plot
roc_results <- plot_roc_magic(
  magic_file = magic_results,
  allsites = allsites,
  truth_labels = truth_labels,
  methdiff_df = methdiff_df,
  method = "wald",
  out_file = "synthetic_example/roc_curve_comparison.png",
  delta_threshold = 0.1
)

cat("\n")

# ==============================================================================
# Alternative: QQ Plot for Type I Error Validation (Null Simulation)
# ==============================================================================

cat("=== Alternative: QQ Plot for Type I Error Validation ===\n")

# Generate null simulation data (no true DML)
cat("Generating null simulation dataset...\n")
null_cmd <- paste(
  "Rscript scripts/synthetic_bb_generator_with_DML.r",
  "-K 3",
  "-n 100000",
  "-s 10",
  "--dml_fraction 0",
  "--effect_size 0",
  "-c 10",
  "-f FALSE",
  "-o synthetic_example/null_simulation.qs"
)

system(null_cmd)

cat("✓ Null simulation data generated!\n\n")

# Fit mixture model on null data
cat("Fitting mixture model on null data...\n")
null_fit_cmd <- paste(
  "Rscript R/magicFit.R",
  "--input synthetic_example/null_simulation.qs",
  "--k_components 3",
  "--max_iters 500",
  "--tol 1e-6",
  "--threads 4",
  "--subset_trait group",
  "--subset_trait_val Control"
)

system(null_fit_cmd)

# Test on null data
cat("Testing on null data...\n")
null_test_cmd <- paste(
  "Rscript R/magicTest.R",
  "--input synthetic_example/null_simulation.qs",
  "--group group",
  "--compare Control,Treatment",
  "--mixture_model synthetic_example/null_simulation_K3_methOnly_trait_Control_S5_e06_it500/",
  "--cores 4"
)

system(null_test_cmd)

cat("✓ Null simulation analysis completed!\n\n")

# Function to create QQ plot data
qq_data <- function(pvals, method) {
  n <- length(pvals)
  exp <- -log10((1:n) / n)
  obs <- -log10(sort(pvals))
  data.frame(expected = exp, observed = obs, method = method)
}

# Load null simulation results
magic_null <- fread("synthetic_example/null_simulation_results/magic_results.csv")

# Create QQ plot data
qq_magic_w3 <- qq_data(magic_null$w3P, "MAGIC w3P")

# Create QQ plot
p_qq <- ggplot(qq_magic_w3, aes(x = expected, y = observed)) +
  geom_point(size = 0.8, alpha = 0.8, color = "#377eb8") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "QQ Plot: MAGIC w3 Test Under Null Simulation",
       x = "Expected -log10(p)",
       y = "Observed -log10(p)")

ggsave("synthetic_example/qq_plot_null_simulation.png", p_qq, width = 6, height = 6, dpi = 300)

cat("✓ QQ plot saved to: synthetic_example/qq_plot_null_simulation.png\n\n")


# ==============================================================================
# STEP 6: Summary Statistics
# ==============================================================================

cat("=== STEP 6: Summary Statistics ===\n\n")

# Load MAGIC results
magic_results <- fread("synthetic_example/magic_results/magic_results.csv")

# Data generation summary
cat("Synthetic Data Parameters:\n")
cat(sprintf("  - Number of CpG sites: %d\n", nrow(bs_data)))
cat(sprintf("  - Number of samples: %d (5 Control + 5 Treatment)\n", ncol(bs_data)))
cat(sprintf("  - True DML sites: %d (%.1f%%)\n", 
            length(metadata(bs_data)$true_dml_sites),
            100 * length(metadata(bs_data)$true_dml_sites) / nrow(bs_data)))
cat(sprintf("  - Effect size: 0.1\n"))
cat(sprintf("  - Mean coverage: 10×\n\n"))

# MAGIC detection summary (FDR < 0.05)
magic_sig <- magic_results[magic_results$w3Q < 0.05, ]
magic_sig$id <- paste0(magic_sig$chr, "_", magic_sig$start)
n_detected <- sum(magic_sig$id %in% truth_labels)
n_false_pos <- nrow(magic_sig) - n_detected

cat("MAGIC Results (w3 test, FDR < 0.05):\n")
cat(sprintf("  - Significant CpGs detected: %d\n", nrow(magic_sig)))
cat(sprintf("  - True positives: %d\n", n_detected))
cat(sprintf("  - False positives: %d\n", n_false_pos))
cat(sprintf("  - Sensitivity (Recall): %.3f\n", 
            n_detected / length(metadata(bs_data)$true_dml_sites)))
cat(sprintf("  - Precision (PPV): %.3f\n", 
            n_detected / nrow(magic_sig)))
cat(sprintf("  - F1 Score: %.3f\n\n",
            2 * (n_detected / nrow(magic_sig)) * (n_detected / length(metadata(bs_data)$true_dml_sites)) /
            ((n_detected / nrow(magic_sig)) + (n_detected / length(metadata(bs_data)$true_dml_sites)))))

# Component distribution
comp_vec <- metadata(bs_data)$true_cpg_components
cat("Mixture Component Distribution:\n")
comp_table <- table(comp_vec)
for (i in 1:length(comp_table)) {
  cat(sprintf("  - Component %d: %d CpGs (%.1f%%)\n", 
              i, comp_table[i], 100 * comp_table[i] / sum(comp_table)))
}

cat("\n=== Analysis Complete! ===\n")
cat("\nGenerated files:\n")
cat("  1. synthetic_example/synthetic_bb_100k_s10_d01_cov10.qs\n")
cat("  2. synthetic_example/synthetic_bb_100k_s10_d01_cov10_truth_DML.csv\n")
cat("  3. synthetic_example/synthetic_bb_100k_s10_d01_cov10_trueModel.csv\n")
cat("  4. synthetic_example/magic_results/magic_results.csv\n")
cat("  5. synthetic_example/methylation_histogram_by_component.png\n")
cat("  6. synthetic_example/roc_curve_comparison.png\n")
cat("  7. synthetic_example/null_simulation.qs\n")
cat("  8. synthetic_example/qq_plot_null_simulation.png\n")
cat("\nYou can now benchmark MAGIC against other methods using the ground truth!\n")