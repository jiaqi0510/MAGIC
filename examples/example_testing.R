#!/usr/bin/env Rscript
#
# MAGIC: Methylation Analysis with Genomic Inferred Contexts
#
# Example: Differential methylation testing
# Description: Demonstrates how to perform differential methylation analysis using
#              MAGIC. Covers case-control comparisons, interpreting test results,
#              filtering and subsetting, and result visualization.
#
library(bsseq)
library(qs)

# ==============================================================================
# Setup
# ==============================================================================

# Input data
input_file <- "data/BSseq_example.qs"

# Fitted mixture model (from optimization step)
# This should be the output directory from magicFit.R
mixture_model_dir <- "results/data_K3_S10_e06_it500"

# Check if files exist
if (!file.exists(input_file)) {
  cat("Error: Input file not found:", input_file, "\n")
  cat("Please run example_optimization.R first\n")
  quit(save = "no", status = 1)
}

if (!dir.exists(mixture_model_dir)) {
  cat("Error: Mixture model directory not found:", mixture_model_dir, "\n")
  cat("Please run example_optimization.R first to create the model\n")
  cat("Or update 'mixture_model_dir' to point to your fitted model\n")
  quit(save = "no", status = 1)
}

# ==============================================================================
# Example 1: Basic Case-Control Comparison
# ==============================================================================

cat("Example 1: Case-control differential methylation\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# This assumes your BSseq object has a column called "condition" 
# with values "case" and "control"
# Adjust --group and --compare parameters for your data

cat("Testing for differential methylation between case and control\n\n")

system2("Rscript", c(
  "R/magicTest.R",
  "--input", input_file,
  "--group", "condition",                    # Column name in colData(BSseq)
  "--compare", "case,control",               # Groups to compare
  "--mixture_model", mixture_model_dir,      # Fitted model directory
  "--cores", "4"
))

cat("\nResults saved to: [analysis]_magic_results_[timestamp]/\n")
cat("Main output file: results.csv\n\n")

# ==============================================================================
# Example 2: Reading and Interpreting Results
# ==============================================================================

cat("\nExample 2: Interpreting test results\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Find the most recent results directory
result_dirs <- list.dirs(".", recursive = FALSE)
result_dirs <- result_dirs[grepl("magic_results_", result_dirs)]

if (length(result_dirs) > 0) {
  # Use most recent
  result_dir <- result_dirs[length(result_dirs)]
  results_file <- file.path(result_dir, "results.csv")
  
  if (file.exists(results_file)) {
    results <- read.csv(results_file)
    
    cat("Results file contains", nrow(results), "sites\n\n")
    
    cat("Column descriptions:\n")
    cat("  chr, start, end: Genomic coordinates\n")
    cat("  mean1, mean2: Mean methylation in each group\n")
    cat("  difference: mean1 - mean2\n")
    cat("  bfProb: Bayes factor probability (0-1, higher = more evidence)\n")
    cat("  mixtureP: P-value from mixture-based test\n")
    cat("  domComp1, domComp2: Dominant component in each group\n")
    cat("  entropy1, entropy2: Assignment uncertainty\n\n")
    
    # Summary statistics
    cat("Summary of results:\n")
    cat(sprintf("  Mean methylation group 1: %.3f\n", mean(results$mean1, na.rm = TRUE)))
    cat(sprintf("  Mean methylation group 2: %.3f\n", mean(results$mean2, na.rm = TRUE)))
    cat(sprintf("  Mean difference: %.3f\n", mean(results$difference, na.rm = TRUE)))
    
    # Significance thresholds
    cat("\nSignificant sites:\n")
    sig_bf <- sum(results$bfProb > 0.9, na.rm = TRUE)
    sig_p <- sum(results$mixtureP < 0.05, na.rm = TRUE)
    
    cat(sprintf("  Bayes factor > 0.9: %d sites (%.1f%%)\n", 
                sig_bf, 100 * sig_bf / nrow(results)))
    cat(sprintf("  P-value < 0.05: %d sites (%.1f%%)\n",
                sig_p, 100 * sig_p / nrow(results)))
    
    # Effect size distribution
    cat("\nEffect size distribution:\n")
    large_effects <- sum(abs(results$difference) > 0.2, na.rm = TRUE)
    cat(sprintf("  |difference| > 0.2: %d sites (%.1f%%)\n",
                large_effects, 100 * large_effects / nrow(results)))
    
    # Show top hits
    cat("\nTop 10 sites by Bayes factor:\n")
    top_sites <- results[order(-results$bfProb), ][1:min(10, nrow(results)), ]
    cat(sprintf("%-6s %-10s %-10s %-8s %-8s %-10s %-8s\n",
                "Chr", "Start", "End", "Mean1", "Mean2", "Difference", "BF_Prob"))
    cat(paste(rep("-", 70), collapse = ""), "\n")
    
    for (i in 1:nrow(top_sites)) {
      cat(sprintf("%-6s %-10d %-10d %-8.3f %-8.3f %-10.3f %-8.3f\n",
                  top_sites$chr[i],
                  top_sites$start[i],
                  top_sites$end[i],
                  top_sites$mean1[i],
                  top_sites$mean2[i],
                  top_sites$difference[i],
                  top_sites$bfProb[i]))
    }
    
    # Component switching analysis
    cat("\nComponent assignments:\n")
    component_switches <- sum(results$domComp1 != results$domComp2, na.rm = TRUE)
    cat(sprintf("  Sites with component switch: %d (%.1f%%)\n",
                component_switches, 100 * component_switches / nrow(results)))
    
    cat("\nComponent switch interpretation:\n")
    cat("  - Different components suggest qualitative state change\n")
    cat("  - Same component with difference suggests quantitative change\n")
  }
}

# ==============================================================================
# Example 3: Chromosome-Specific Analysis
# ==============================================================================

cat("\n\nExample 3: Analyze specific chromosome\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

chr_to_test <- "chr21"
cat(sprintf("Testing only %s\n\n", chr_to_test))

system2("Rscript", c(
  "R/magicTest.R",
  "--input", input_file,
  "--group", "condition",
  "--compare", "case,control",
  "--mixture_model", mixture_model_dir,
  "--chr", chr_to_test,
  "--cores", "4"
))

# ==============================================================================
# Example 4: Subset Analysis
# ==============================================================================

cat("\n\nExample 4: Subset analysis (tissue-specific)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# This example assumes your BSseq object has a "tissue" column
# and you want to analyze only "liver" samples

cat("Analyzing only samples from specific tissue type\n\n")

system2("Rscript", c(
  "R/magicTest.R",
  "--input", input_file,
  "--group", "condition",
  "--compare", "case,control",
  "--mixture_model", mixture_model_dir,
  "--subset", "tissue",        # Column to filter on
  "--keep", "liver",           # Value to keep
  "--cores", "4"
))

cat("\nOnly samples with tissue='liver' were included in analysis\n")

# ==============================================================================
# Example 5: Multiple Testing Correction
# ==============================================================================

cat("\n\nExample 5: Multiple testing correction\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Find most recent results
result_dirs <- list.dirs(".", recursive = FALSE)
result_dirs <- result_dirs[grepl("magic_results_", result_dirs)]

if (length(result_dirs) > 0) {
  result_dir <- result_dirs[length(result_dirs)]
  results_file <- file.path(result_dir, "results.csv")
  
  if (file.exists(results_file)) {
    results <- read.csv(results_file)
    
    # Apply Benjamini-Hochberg FDR correction
    results$fdr <- p.adjust(results$mixtureP, method = "BH")
    
    cat("FDR correction applied (Benjamini-Hochberg method)\n\n")
    
    # Count significant at different thresholds
    cat("Significant sites at different FDR thresholds:\n")
    for (thresh in c(0.01, 0.05, 0.10)) {
      n_sig <- sum(results$fdr < thresh, na.rm = TRUE)
      cat(sprintf("  FDR < %.2f: %d sites\n", thresh, n_sig))
    }
    
    # Save results with FDR
    output_file <- file.path(result_dir, "results_with_fdr.csv")
    write.csv(results, output_file, row.names = FALSE)
    cat(sprintf("\nResults with FDR saved to: %s\n", output_file))
    
    # Extract significant sites
    sig_file <- file.path(result_dir, "significant_sites_fdr05.csv")
    sig_results <- results[results$fdr < 0.05, ]
    if (nrow(sig_results) > 0) {
      write.csv(sig_results, sig_file, row.names = FALSE)
      cat(sprintf("Significant sites (FDR < 0.05) saved to: %s\n", sig_file))
    }
  }
}

# ==============================================================================
# Example 6: Null Distribution Testing
# ==============================================================================

cat("\n\nExample 6: Null distribution (randomized groups)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Running analysis with randomized group labels\n")
cat("This creates null distribution for calibrating FDR\n\n")

system2("Rscript", c(
  "R/magicTest.R",
  "--input", input_file,
  "--group", "condition",
  "--compare", "case,control",
  "--mixture_model", mixture_model_dir,
  "--randomize",               # Randomize group assignments
  "--cores", "4"
))

cat("\nCompare null p-value distribution to real analysis\n")
cat("Expected: uniform distribution under null\n")

# ==============================================================================
# Example 7: Export Results for Visualization
# ==============================================================================

cat("\n\nExample 7: Preparing results for visualization\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

result_dirs <- list.dirs(".", recursive = FALSE)
result_dirs <- result_dirs[grepl("magic_results_", result_dirs)]

if (length(result_dirs) > 0) {
  result_dir <- result_dirs[length(result_dirs)]
  results_file <- file.path(result_dir, "results.csv")
  
  if (file.exists(results_file)) {
    results <- read.csv(results_file)
    
    # Add FDR if not already present
    if (!"fdr" %in% names(results)) {
      results$fdr <- p.adjust(results$mixtureP, method = "BH")
    }
    
    # Create BED file for genome browser
    bed_file <- file.path(result_dir, "significant_sites.bed")
    sig_results <- results[results$fdr < 0.05 & abs(results$difference) > 0.2, ]
    
    if (nrow(sig_results) > 0) {
      # BED format: chr, start, end, name, score
      bed_data <- data.frame(
        chr = sig_results$chr,
        start = sig_results$start,
        end = sig_results$end,
        name = paste0("site_", 1:nrow(sig_results)),
        score = round(1000 * sig_results$bfProb)  # Scale to 0-1000
      )
      
      write.table(bed_data, bed_file, sep = "\t", 
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
      
      cat(sprintf("BED file created: %s\n", bed_file))
      cat(sprintf("  Contains %d sites (FDR < 0.05, |diff| > 0.2)\n", nrow(sig_results)))
      cat("  Load in IGV or UCSC Genome Browser for visualization\n\n")
    }
    
    # Create summary plot data
    cat("Creating summary statistics for plotting:\n")
    
    # Volcano plot data
    volcano_file <- file.path(result_dir, "volcano_data.csv")
    volcano_data <- data.frame(
      difference = results$difference,
      neg_log10_p = -log10(results$mixtureP),
      significant = results$fdr < 0.05
    )
    write.csv(volcano_data, volcano_file, row.names = FALSE)
    cat(sprintf("  Volcano plot data: %s\n", volcano_file))
    
    # Manhattan plot data
    manhattan_file <- file.path(result_dir, "manhattan_data.csv")
    manhattan_data <- data.frame(
      chr = results$chr,
      position = results$start,
      neg_log10_p = -log10(results$mixtureP)
    )
    write.csv(manhattan_data, manhattan_file, row.names = FALSE)
    cat(sprintf("  Manhattan plot data: %s\n", manhattan_file))
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
cat("  1. Basic differential methylation testing\n")
cat("  2. Interpreting Bayes factors and p-values\n")
cat("  3. Chromosome-specific and subset analyses\n")
cat("  4. Multiple testing correction\n")
cat("  5. Null distribution testing\n")
cat("  6. Exporting results for visualization\n\n")

cat("Key result metrics:\n")
cat("  - bfProb > 0.9: Strong evidence for difference\n")
cat("  - FDR < 0.05: Statistically significant after correction\n")
cat("  - |difference| > 0.2: Moderate effect size\n")
cat("  - Component switch: Qualitative state change\n\n")

cat("Next steps:\n")
cat("  1. Apply to your own data\n")
cat("  2. Visualize results in R or genome browser\n")
cat("  3. Perform downstream enrichment analyses\n")
cat("  4. See docs/usage.md for all options\n\n")
