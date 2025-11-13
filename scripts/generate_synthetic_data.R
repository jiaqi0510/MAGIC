#!/usr/bin/env Rscript
#
# MAGIC: Methylation Analysis with Genomic Inferred Contexts
#
# Module: Synthetic data generation
# Description: Generates synthetic methylation data from known beta-binomial mixture
#              models for testing and validation. Creates realistic data with specified
#              component parameters, coverage distributions, and sample sizes. Saves
#              both data and true model parameters for benchmarking.
#
# Authors: Jiaqi Han, Michael Thompson, Matteo Pellegrini
# Contact: mjthompson69@gmail.com
# License: MIT
# Date: 2025

suppressPackageStartupMessages({
  library(bsseq)
  library(qs)
})

generate_synthetic_bb_mixture <- function(K, n_cpgs = 1000000, n_samples = 100, 
                                        output_file = "synthetic_bb_data.qs",
                                        seed = 12345, apply_filtering = TRUE,
                                        min_reads = 5, min_samples = 6, min_var = 0.001,
                                        sparse_coverage = FALSE, sparsity_level = 0.8,
                                        sex_diff = FALSE, chrx_fraction = 0.05) {
  
  set.seed(seed)
  
  if (apply_filtering) {
    n_cpgs_initial <- as.integer(n_cpgs * 2.5)
    cat("Generating", n_cpgs_initial, "initial CpGs to account for filtering\n")
  } else {
    n_cpgs_initial <- n_cpgs
  }
  
  # Generate realistic beta-binomial parameters for K components
  generate_realistic_params <- function(K) {
    params <- list()
    
    if (K == 1) {
      params[[1]] <- list(alpha = 2, beta = 2)
    } else if (K == 2) {
      params[[1]] <- list(alpha = 1, beta = 10)
      params[[2]] <- list(alpha = 10, beta = 1)
    } else if (K == 3) {
      params[[1]] <- list(alpha = 1, beta = 15)
      params[[2]] <- list(alpha = 3, beta = 3)
      params[[3]] <- list(alpha = 15, beta = 1)
    } else {
      alpha_seq <- seq(1, 20, length.out = K)
      beta_seq <- seq(20, 1, length.out = K)
      for (i in 1:K) {
        params[[i]] <- list(alpha = alpha_seq[i], beta = beta_seq[i])
      }
    }
    
    return(params)
  }
  
  # Generate mixing proportions
  mixing_props <- rep(1/K, K)
  
  # Generate component parameters
  component_params <- generate_realistic_params(K)
  
  # Generate sample sex assignments if sex differences requested
  sample_sex <- NULL
  if (sex_diff) {
    sample_sex <- sample(c("Male", "Female"), n_samples, replace = TRUE, prob = c(0.5, 0.5))
    cat("Generating sex-specific differential methylation\n")
    cat("Sample sex distribution:", table(sample_sex), "\n")
  }
  
  # Assign each CpG to a component
  cpg_components <- sample(1:K, n_cpgs_initial, replace = TRUE, prob = mixing_props)
  
  # Designate chrX sites if sex differences requested
  chrx_sites <- NULL
  if (sex_diff) {
    n_chrx_sites <- round(n_cpgs_initial * chrx_fraction)
    chrx_sites <- sample(n_cpgs_initial, n_chrx_sites)
    cat("Designated", n_chrx_sites, "CpGs as chrX sites (", round(chrx_fraction*100, 1), "%)\n")
  }
  
  # Generate synthetic genomic coordinates
  chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
  chr_weights <- c(rep(1, 22), 0.5, 0.1)
  
  cpg_chrs <- sample(chrs, n_cpgs_initial, replace = TRUE, prob = chr_weights)
  
  # Override chromosome assignment for designated chrX sites
  if (sex_diff && !is.null(chrx_sites)) {
    cpg_chrs[chrx_sites] <- "chrX"
  }
  
  # Generate positions within chromosomes
  cpg_positions <- integer(n_cpgs_initial)
  for (chr in unique(cpg_chrs)) {
    chr_mask <- cpg_chrs == chr
    n_chr_cpgs <- sum(chr_mask)
    if (chr %in% paste0("chr", 1:22)) {
      max_pos <- 200000000
    } else if (chr == "chrX") {
      max_pos <- 150000000
    } else {
      max_pos <- 50000000
    }
    cpg_positions[chr_mask] <- sort(sample(1:max_pos, n_chr_cpgs, replace = FALSE))
  }
  
  # Create GRanges object for CpG positions
  cpg_granges <- GRanges(
    seqnames = cpg_chrs,
    ranges = IRanges(start = cpg_positions, width = 1)
  )
  
  # Generate coverage data (vectorized)
  coverage_matrix <- matrix(0, nrow = n_cpgs_initial, ncol = n_samples)
  
  if (sparse_coverage) {
    cat("Generating sparse WGBS-like coverage with", sparsity_level*100, "% zeros\n")
    
    n_total_entries <- n_cpgs_initial * n_samples
    n_nonzero_entries <- round(n_total_entries * (1 - sparsity_level))
    
    # Vectorized approach: generate all nonzero values at once
    nonzero_positions <- sample(n_total_entries, n_nonzero_entries)
    nonzero_coverage <- pmax(1, round(rgamma(n_nonzero_entries, shape = 2, scale = 10)))
    
    # Use linear indexing to assign coverage values
    coverage_matrix[nonzero_positions] <- nonzero_coverage
    
    cat("Generated coverage: mean =", round(mean(coverage_matrix), 1), 
        ", sparsity =", round(sum(coverage_matrix == 0) / n_total_entries, 3), "\n")
    
  } else {
    # Dense coverage generation (vectorized)
    coverage_matrix[] <- rnbinom(n_cpgs_initial * n_samples, size = 10, mu = 20)
  }
  
  # Generate methylation data (fully vectorized)
  meth_matrix <- matrix(0, nrow = n_cpgs_initial, ncol = n_samples)
  
  if (sex_diff && !is.null(chrx_sites)) {
    # Handle chrX sites with sex-specific patterns (vectorized by sex)
    male_samples <- which(sample_sex == "Male")
    female_samples <- which(sample_sex == "Female")
    
    # For chrX sites in males: all get low methylation
    if (length(male_samples) > 0) {
      chrx_coverage_male <- coverage_matrix[chrx_sites, male_samples, drop = FALSE]
      valid_male <- chrx_coverage_male > 0
      n_valid_male <- sum(valid_male)
      
      if (n_valid_male > 0) {
        p_meth_male <- rbeta(n_valid_male, 1, 15)  # Low methylation
        meth_matrix[chrx_sites, male_samples][valid_male] <- rbinom(n_valid_male, 
                                                                   chrx_coverage_male[valid_male], 
                                                                   p_meth_male)
      }
    }
    
    # For chrX sites in females: mix of low and intermediate
    if (length(female_samples) > 0) {
      chrx_coverage_female <- coverage_matrix[chrx_sites, female_samples, drop = FALSE]
      valid_female <- chrx_coverage_female > 0
      n_valid_female <- sum(valid_female)
      
      if (n_valid_female > 0) {
        # 60% low methylation, 40% intermediate
        x_inact_pattern <- sample(c(TRUE, FALSE), n_valid_female, replace = TRUE, prob = c(0.6, 0.4))
        p_meth_female <- ifelse(x_inact_pattern, 
                               rbeta(n_valid_female, 1, 15),  # Low methylation (active X)
                               rbeta(n_valid_female, 3, 3))   # Intermediate (inactive X)
        
        meth_matrix[chrx_sites, female_samples][valid_female] <- rbinom(n_valid_female,
                                                                       chrx_coverage_female[valid_female],
                                                                       p_meth_female)
      }
    }
    
    # For non-chrX sites: use regular component assignment
    non_chrx_sites <- setdiff(1:n_cpgs_initial, chrx_sites)
  } else {
    # No sex differences: all sites use regular component assignment
    non_chrx_sites <- 1:n_cpgs_initial
  }
  
  # Generate methylation for non-chrX sites (vectorized by component)
  for (k in 1:K) {
    comp_mask <- cpg_components %in% k
    if (sex_diff) {
      comp_mask <- comp_mask & (1:n_cpgs_initial %in% non_chrx_sites)
    }
    
    comp_sites <- which(comp_mask)
    if (length(comp_sites) > 0) {
      alpha <- component_params[[k]]$alpha
      beta <- component_params[[k]]$beta
      
      # Get coverage for this component's sites
      comp_coverage <- coverage_matrix[comp_sites, , drop = FALSE]
      valid_entries <- comp_coverage > 0
      n_valid <- sum(valid_entries)
      
      if (n_valid > 0) {
        # Generate all methylation probabilities for this component at once
        p_meth_comp <- rbeta(n_valid, alpha, beta)
        
        # Generate methylation counts for all valid entries at once
        meth_matrix[comp_sites, ][valid_entries] <- rbinom(n_valid, 
                                                          comp_coverage[valid_entries], 
                                                          p_meth_comp)
      }
    }
  }
  
  # Apply filtering if requested
  if (apply_filtering) {
    cat("Applying magicFit filtering...\n")
    
    # Vectorized filtering operations
    valid_mask <- coverage_matrix >= min_reads
    samples_per_feat <- rowSums(valid_mask)
    valid_feats <- samples_per_feat >= min_samples
    
    cat("After coverage/sample filters:", sum(valid_feats), "sites remain\n")
    
    # Apply coverage filter
    meth_matrix_filt <- meth_matrix[valid_feats, , drop = FALSE]
    coverage_matrix_filt <- coverage_matrix[valid_feats, , drop = FALSE]
    valid_mask_filt <- valid_mask[valid_feats, , drop = FALSE]
    cpg_components_filt <- cpg_components[valid_feats]
    cpg_granges_filt <- cpg_granges[valid_feats]
    
    if (sex_diff && !is.null(chrx_sites)) {
      # Update chrX sites indices after filtering
      old_to_new_idx <- cumsum(valid_feats)
      chrx_sites_filt <- old_to_new_idx[chrx_sites[chrx_sites %in% which(valid_feats)]]
      chrx_sites <- chrx_sites_filt[chrx_sites_filt > 0]
    }
    
    # Variance filter (vectorized)
    props_matrix <- meth_matrix_filt / pmax(coverage_matrix_filt, 1)
    props_matrix[!valid_mask_filt] <- NA
    
    # Vectorized variance calculation
    vars <- apply(props_matrix, 1, function(x) var(x, na.rm = TRUE))
    vars[is.na(vars)] <- 0
    
    valid_vars <- vars[vars > 0 & is.finite(vars)]
    var_thresh <- max(min_var, quantile(valid_vars, 0.25, na.rm = TRUE))
    var_pass <- vars >= var_thresh & is.finite(vars)
    
    cat("After variance filter:", sum(var_pass), "sites remain\n")
    
    if (sum(var_pass) < n_cpgs) {
      n_final <- sum(var_pass)
    } else {
      var_pass_indices <- which(var_pass)[1:n_cpgs]
      var_pass[] <- FALSE
      var_pass[var_pass_indices] <- TRUE
      n_final <- n_cpgs
    }
    
    # Apply final filtering
    meth_matrix <- meth_matrix_filt[var_pass, , drop = FALSE]
    coverage_matrix <- coverage_matrix_filt[var_pass, , drop = FALSE]
    cpg_components <- cpg_components_filt[var_pass]
    cpg_granges <- cpg_granges_filt[var_pass]
    
    if (sex_diff && !is.null(chrx_sites)) {
      # Final update of chrX sites after variance filtering
      var_pass_indices <- which(var_pass)
      chrx_sites <- which(var_pass_indices %in% chrx_sites)
    }
    
    # Recalculate mixing proportions (vectorized)
    filtered_component_counts <- tabulate(cpg_components, nbins = K)
    mixing_props <- filtered_component_counts / sum(filtered_component_counts)
    
    cat("Post-filtering mixing proportions:\n")
    for (i in 1:K) {
      cat(sprintf("  Component %d: %.3f (was %.3f)\n", i, mixing_props[i], 1/K))
    }
    
  } else {
    n_final <- n_cpgs_initial
  }
  
  # Create sample names and metadata
  sample_names <- paste0("Sample_", sprintf("%03d", 1:n_samples))
  
  # Create sample metadata data frame
  sample_metadata <- data.frame(
    sampleNames = sample_names,
    stringsAsFactors = FALSE
  )
  rownames(sample_metadata) <- sample_names
  
  if (sex_diff) {
    sample_metadata$sex <- sample_sex
  }
  
  # Create bsseq object
  bs_obj <- BSseq(
    gr = cpg_granges,
    M = meth_matrix,
    Cov = coverage_matrix,
    sampleNames = sample_names,
    pData = sample_metadata
  )
  
  # Store true parameters as metadata
  metadata(bs_obj)$true_K <- K
  metadata(bs_obj)$true_mixing_props <- mixing_props
  metadata(bs_obj)$true_component_params <- component_params
  metadata(bs_obj)$true_cpg_components <- cpg_components
  
  if (sex_diff) {
    metadata(bs_obj)$sex_diff_enabled <- TRUE
    metadata(bs_obj)$chrx_sites <- chrx_sites
    metadata(bs_obj)$sample_sex <- sample_sex
    metadata(bs_obj)$chrx_fraction <- chrx_fraction
  }
  
  # Save to file
  qsave(bs_obj, output_file)
  
  # Save true model parameters to CSV
  true_model_file <- sub("\\.qs$", "_trueModel.csv", output_file)
  true_model_df <- data.frame(
    component = 1:K,
    alpha = sapply(component_params, function(x) x$alpha),
    beta = sapply(component_params, function(x) x$beta),
    pi = mixing_props
  )
  write.csv(true_model_df, true_model_file, row.names = FALSE)
  
  # Print summary
  cat("Generated synthetic beta-binomial mixture data:\n")
  cat("  K =", K, "components\n")
  cat("  CpGs:", nrow(meth_matrix), "\n")
  cat("  Samples:", n_samples, "\n")
  cat("  Mean coverage:", round(mean(coverage_matrix), 1), "\n")
  cat("  Mean methylation:", round(mean(meth_matrix / pmax(coverage_matrix, 1)), 3), "\n")
  cat("  Saved to:", output_file, "\n")
  cat("  True model saved to:", true_model_file, "\n")
  
  # Print component parameters
  cat("\nTrue component parameters:\n")
  for (i in 1:K) {
    alpha <- component_params[[i]]$alpha
    beta <- component_params[[i]]$beta
    mean_meth <- alpha / (alpha + beta)
    cat(sprintf("  Component %d: alpha=%.1f, beta=%.1f, mean=%.3f\n", 
                i, alpha, beta, mean_meth))
  }
  
  return(bs_obj)
}

# Command line interface
args <- commandArgs(trailingOnly = TRUE)

show_help <- function() {
  cat("Generate synthetic beta-binomial mixture methylation data\n\n")
  cat("Usage: Rscript synthetic_bb_generator.R -K <int> [options]\n")
  cat("       Rscript synthetic_bb_generator.R --help\n\n")
  cat("Required:\n")
  cat("  -K, --components <int>    Number of mixture components\n\n")
  cat("Options:\n")
  cat("  -n, --cpgs <int>         Number of CpG sites (default: 1000000)\n")
  cat("  -s, --samples <int>      Number of samples (default: 100)\n")
  cat("  -o, --output <file>      Output .qs file (default: synthetic_bb_K{K}.qs)\n")
  cat("  --seed <int>             Random seed (default: 12345)\n")
  cat("  --no-filter              Skip magicFit filtering\n")
  cat("  --min-reads <int>        Min reads for filtering (default: 5)\n")
  cat("  --min-samples <int>      Min samples for filtering (default: 6)\n")
  cat("  --min-var <num>          Min variance for filtering (default: 0.001)\n")
  cat("  --sparse                 Generate sparse coverage (WGBS-like)\n")
  cat("  --sparsity <num>         Fraction of zero coverage entries (default: 0.8)\n")
  cat("  --sex-diff               Add sex-specific differential methylation\n")
  cat("  --chrx-fraction <num>    Fraction of CpGs designated as chrX (default: 0.05)\n")
  cat("  -h, --help               Show this help\n\n")
  cat("Examples:\n")
  cat("  Rscript synthetic_bb_generator.R -K 3\n")
  cat("  Rscript synthetic_bb_generator.R -K 3 -n 500000 -s 50 -o my_data.qs --seed 42\n")
  cat("  Rscript synthetic_bb_generator.R -K 3 --sparse --sparsity 0.43\n")
  cat("  Rscript synthetic_bb_generator.R -K 3 --sex-diff --chrx-fraction 0.08\n")
}

if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  show_help()
  quit(status = 0)
}

# Parse flags
K <- NULL
n_cpgs <- 1000000
n_samples <- 100
output_file <- NULL
seed <- 12345
apply_filtering <- TRUE
min_reads <- 5
min_samples <- 6
min_var <- 0.001
sparse_coverage <- FALSE
sparsity_level <- 0.8
sex_diff <- FALSE
chrx_fraction <- 0.05

i <- 1
while (i <= length(args)) {
  if (args[i] == "-K" || args[i] == "--components") {
    K <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "-n" || args[i] == "--cpgs") {
    n_cpgs <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "-s" || args[i] == "--samples") {
    n_samples <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "-o" || args[i] == "--output") {
    output_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--seed") {
    seed <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--no-filter") {
    apply_filtering <- FALSE
    i <- i + 1
  } else if (args[i] == "--min-reads") {
    min_reads <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--min-samples") {
    min_samples <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--min-var") {
    min_var <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--sparse") {
    sparse_coverage <- TRUE
    i <- i + 1
  } else if (args[i] == "--sparsity") {
    sparsity_level <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--sex-diff") {
    sex_diff <- TRUE
    i <- i + 1
  } else if (args[i] == "--chrx-fraction") {
    chrx_fraction <- as.numeric(args[i + 1])
    i <- i + 2
  } else {
    stop(paste("Unknown argument:", args[i]))
  }
}

if (is.null(K)) {
  stop("Required argument -K/--components not provided")
}

if (is.null(output_file)) {
  output_file <- paste0("synthetic_bb_K", K, ".qs")
}

# Generate data
cat("Generating synthetic beta-binomial mixture data...\n")
bs_data <- generate_synthetic_bb_mixture(
  K = K,
  n_cpgs = n_cpgs,
  n_samples = n_samples,
  output_file = output_file,
  seed = seed,
  apply_filtering = apply_filtering,
  min_reads = min_reads,
  min_samples = min_samples,
  min_var = min_var,
  sparse_coverage = sparse_coverage,
  sparsity_level = sparsity_level,
  sex_diff = sex_diff,
  chrx_fraction = chrx_fraction
)
