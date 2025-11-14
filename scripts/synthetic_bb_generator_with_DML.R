#!/usr/bin/env Rscript
#
# MAGIC: Methylation Analysis with Genomic Inferred Contexts
#
# Module: Synthetic data generation
# Description: Generates synthetic methylation data from known beta-binomial mixture
#              models for testing and validation. Creates realistic data with specified
#              component parameters, coverage distributions, and sample sizes. Controlled 
#              injection of Differentially Methylated Loci (DML) between groups. Support
#              for null simulations (no DML) for Type I error rate evaluation. 
#              Saves data, true model parameters and true DMLs for benchmarking.
#
# Authors: Jiaqi Han, Michael Thompson, Matteo Pellegrini
# Contact: mjthompson69@gmail.com
# License: MIT
# Date: 2025


# Load required libraries with suppressed startup messages
suppressPackageStartupMessages({
  library(bsseq)      # Bioconductor package for bisulfite sequencing data
  library(qs)         # Fast serialization for R objects
  library(optparse)   # Command-line argument parsing
})

#' Generate Synthetic WGBS Data with Beta-Binomial Mixture Model and DML Injection
#'
#' This function creates realistic synthetic WGBS data by modeling methylation 
#' patterns as mixtures of Beta distributions, then injecting controlled 
#' differential methylation between experimental groups.
#'
#' @param K Integer. Number of Beta distribution components in the mixture model
#' @param n_cpgs Integer. Target number of CpG sites to include in final dataset
#' @param n_samples Integer. Total number of samples (must be even, split equally between groups)
#' @param output_file Character. Output filename for BSseq object (.qs extension)
#' @param seed Integer. Random seed for reproducible results
#' @param apply_filtering Logical. Whether to oversample CpGs to account for filtering
#' @param min_reads Integer. Minimum read coverage threshold for filtering
#' @param min_samples Integer. Minimum number of samples with adequate coverage
#' @param min_var Numeric. Minimum methylation variance threshold for filtering
#' @param sparse_coverage Logical. Whether to simulate sparse coverage patterns typical of WGBS
#' @param sparsity_level Numeric. Fraction of CpG-sample combinations with zero coverage (0-1)
#' @param sex_diff Logical. Whether to simulate sex-specific methylation differences (future feature)
#' @param chrx_fraction Numeric. Proportion of CpGs to place on chromosome X
#' @param dml_fraction Numeric. Fraction of CpGs to designate as differentially methylated (0-1)
#' @param effect_size Numeric. Magnitude of methylation difference for DML sites (0-1)
#' @param true_delta Numeric. Minimum effect size threshold for DML classification
#' @param dml_direction Character. Direction of methylation change ("increase" or "decrease")
#' @param mean_coverage Numeric. Expected mean sequencing coverage per CpG site
#' @param coverage_dispersion Numeric. Dispersion parameter for negative binomial coverage
#' @param mixing_props Numeric vector. Custom mixture proportions (must sum to 1, length K)
#'
#' @return BSseq object containing synthetic methylation data with metadata
#'
#' @examples
#' # Generate basic 3-component mixture with 10% DML
#' bs_data <- generate_synthetic_bb_mixture_with_DML(K = 3, n_cpgs = 10000, n_samples = 20)
#'
#' # Null simulation for Type I error evaluation  
#' null_data <- generate_synthetic_bb_mixture_with_DML(K = 3, dml_fraction = 0)
#'
generate_synthetic_bb_mixture_with_DML <- function(K = 3, 
                                                   n_cpgs = 1000000, 
                                                   n_samples = 100,
                                                   output_file = "synthetic_bb_dml.qs",
                                                   seed = 12345, 
                                                   apply_filtering = TRUE,
                                                   min_reads = 5, 
                                                   min_samples = 6, 
                                                   min_var = 0.001,
                                                   sparse_coverage = FALSE, 
                                                   sparsity_level = 0.8,
                                                   sex_diff = FALSE, 
                                                   chrx_fraction = 0.05,
                                                   dml_fraction = 0.2, 
                                                   effect_size = 0.2,
                                                   true_delta = 0.1, 
                                                   dml_direction = "increase",
                                                   mean_coverage = 20, 
                                                   coverage_dispersion = 10,
                                                   mixing_props = NULL) {
  
  # Set random seed for reproducible results across all random operations
  set.seed(seed)
  
  # Calculate Initial CpG Count with Oversampling for Filtering
  
  # If filtering is enabled, generate extra CpGs to compensate for sites 
  # that will be removed due to low coverage or low variance
  if (apply_filtering) {
    n_cpgs_initial <- as.integer(n_cpgs * 2.5)  # 2.5x oversampling factor
    cat("Generating", n_cpgs_initial, "initial CpGs to account for post-filtering requirements\n")
  } else {
    n_cpgs_initial <- n_cpgs
  }
  
  # Define Beta Distribution Parameters for Mixture Components
  
  #' Generate Biologically Realistic Beta Distribution Parameters
  #' 
  #' Creates parameter sets that reflect typical methylation patterns:
  #' - Low methylation (near 0): unmethylated CpG islands, active promoters
  #' - Intermediate methylation: partially methylated domains, tissue-specific sites  
  #' - High methylation (near 1): repetitive elements, heterochromatin
  #'
  #' @param K Number of mixture components
  #' @return List of alpha/beta parameter pairs for each component

  generate_realistic_params <- function(K = 5, seed = 123) {
    set.seed(seed)
    
    # Define mean methylation levels (μ) evenly spaced from low to high
    means <- seq(0.1, 0.9, length.out = K)
    
    # Define symmetric dispersion (α + β) pattern:
    # Example for K = 5: [80, 40, 20, 40, 80]
    center_strength <- 20
    edge_strength   <- 100

    # Create symmetric dispersion pattern
    if (K %% 2 == 1) {
      # K=3: [100, 20, 100]
      # K=5: [100, 60, 20, 60, 100]
      mid_idx <- ceiling(K / 2)
      left_half <- seq(edge_strength, center_strength, length.out = mid_idx)
      
      right_half <- rev(left_half[-length(left_half)])
      total_strengths <- c(left_half, right_half)
    } else {
      # symmetric between two center points
      left <- seq(edge_strength, center_strength, length.out = K / 2)
      total_strengths <- c(left, rev(left))
    }
    cat("Before noise:", total_strengths, "\n")        
    
    # Compute Beta distribution parameters (α, β) for each component
    params <- list()
    for (i in seq_len(K)) {
      mean_i <- means[i]
      strength_i <- total_strengths[i]
      
      alpha_i <- mean_i * strength_i
      beta_i  <- (1 - mean_i) * strength_i
      
      params[[i]] <- list(
        mean = mean_i,
        strength = strength_i,
        alpha = alpha_i,
        beta = beta_i
      )
    }
    
    return(params)
  }

  # Set Up Mixture Component Proportions
  
  # Use uniform mixing proportions if not specified by user
  if (is.null(mixing_props)) {
    mixing_props <- rep(1/K, K)  # Equal probability for each component
  } else {
    # Validate user-provided mixing proportions
    if (length(mixing_props) != K) {
      stop("Length of mixing_props must equal K (number of components)")
    }
    # Normalize to ensure proportions sum to 1
    mixing_props <- mixing_props / sum(mixing_props)
  }
  
  cat("Using mixture proportions:", paste(round(mixing_props, 3), collapse = ", "), "\n")

  # Generate component parameters and randomly assign CpGs to components
  component_params <- generate_realistic_params(K)
  cpg_components   <- sample(1:K, n_cpgs_initial, replace = TRUE, prob = mixing_props)
  
  # Generate Genomic Coordinates for CpG Sites
  
  # Define human chromosomes with realistic relative frequencies
  chrs <- c(paste0("chr", 1:22), "chrX", "chrY")  # Autosomes + sex chromosomes
  chr_weights <- c(rep(1, 22), 0.5, 0.1)          # Reduced weight for sex chromosomes
  
  # Randomly assign chromosomes to CpG sites based on weights
  cpg_chrs <- sample(chrs, n_cpgs_initial, replace = TRUE, prob = chr_weights)
  
  # Generate realistic genomic positions for each chromosome
  cpg_positions <- integer(n_cpgs_initial)
  for (chr in unique(cpg_chrs)) {
    chr_mask <- cpg_chrs == chr
    n_chr_cpgs <- sum(chr_mask)
    
    # Set chromosome-specific maximum positions (approximate human genome sizes)
    max_pos <- if (chr %in% paste0("chr", 1:22)) {
      200000000    # ~200 Mb for autosomes
    } else if (chr == "chrX") {
      150000000    # ~150 Mb for chromosome X  
    } else {
      50000000     # ~50 Mb for chromosome Y
    }
    
    # Generate sorted positions to maintain genomic order within chromosomes
    cpg_positions[chr_mask] <- sort(sample(1:max_pos, n_chr_cpgs, replace = FALSE))
  }
  
  # Create GenomicRanges object for CpG coordinates
  cpg_granges <- GRanges(seqnames = cpg_chrs, 
                         ranges = IRanges(start = cpg_positions, width = 1))
  
  # Simulate Sequencing Coverage Patterns
  
  # Initialize coverage matrix: rows = CpG sites, columns = samples
  coverage_matrix <- matrix(0, nrow = n_cpgs_initial, ncol = n_samples)
  
  if (sparse_coverage) {
    # Sparse coverage mode: simulates real WGBS data with many zero-coverage sites
    cat("Generating sparse coverage (WGBS-like) with", sparsity_level*100, "% zero coverage\n")
    
    n_total <- n_cpgs_initial * n_samples
    n_nonzero <- round(n_total * (1 - sparsity_level))
    
    # Randomly select positions to have non-zero coverage
    nonzero_positions <- sample(n_total, n_nonzero)
    
    # Generate non-zero coverage values using Gamma distribution (realistic for NGS)
    nonzero_coverage <- pmax(1, round(rgamma(n_nonzero, shape = 2, scale = mean_coverage / 2)))
    
    # Assign non-zero coverage values to selected positions
    coverage_matrix[nonzero_positions] <- nonzero_coverage
    
  } else {
    # Dense coverage mode: all sites have some coverage (Negative Binomial distribution)
    coverage_matrix[] <- rnbinom(n_cpgs_initial * n_samples, 
                                size = coverage_dispersion, 
                                mu = mean_coverage)
  }
  
  # Simulate Methylation Proportions and Read Counts  
  
  # Initialize matrices for methylation proportions and observed counts
  p_matrix    <- matrix(0, nrow = n_cpgs_initial, ncol = n_samples)  # True methylation proportions
  meth_matrix <- matrix(0, nrow = n_cpgs_initial, ncol = n_samples)  # Observed methylated read counts
  
  # Generate methylation data for each mixture component
  for (k in 1:K) {
    component_sites <- which(cpg_components == k)
    
    if (length(component_sites) > 0) {
      alpha <- component_params[[k]]$alpha
      beta  <- component_params[[k]]$beta
      
      # For each CpG site in this component
      for (site in component_sites) {
        # Draw methylation proportions from Beta distribution for all samples
        p_values <- rbeta(n_samples, alpha, beta)
        p_matrix[site, ] <- p_values
        
        # Generate observed methylated read counts using Binomial distribution
        # Each sample's counts depend on its coverage and true methylation proportion
        coverage_values <- coverage_matrix[site, ]
        meth_matrix[site, ] <- rbinom(n_samples, coverage_values, p_values)
      }
    }
  }
  
  # Inject Differentially Methylated Loci (DML) Between Groups
  
  # Calculate number of DML sites to inject
  n_dml <- round(dml_fraction * n_cpgs_initial)
  
  if (n_dml == 0) {
    # Null simulation: no differential methylation (for Type I error evaluation)
    cat("No DML will be injected (dml_fraction = 0) - null simulation mode\n")
    dml_sites <- integer(0)
    
    # Create empty truth table for consistency
    truth_df <- data.frame(
      CpG_ID = integer(0),
      chr = character(0),
      pos = integer(0),
      Effect = numeric(0),
      Direction = character(0),
      mean_meth_control = numeric(0),
      mean_meth_treatment = numeric(0),
      min_cov_control = numeric(0),
      min_cov_treatment = numeric(0),
      delta = numeric(0)
    )
    
    # Save empty truth table
    write.csv(truth_df, sub("\\.qs$", "_truth_DML.csv", output_file), row.names = FALSE)
    
  } else {
    # Define experimental groups (assuming equal group sizes)
    group_A <- 1:(n_samples/2)              # Control group indices  
    group_B <- (n_samples/2 + 1):n_samples  # Treatment group indices
    
    # Calculate baseline methylation levels in control group
    mean_p_control <- rowMeans(p_matrix[, group_A])
    
    # Select eligible sites based on direction of methylation change
    if (dml_direction == "increase") {
      # For methylation increases, avoid already highly methylated sites
      eligible_sites <- which(mean_p_control < 0.8)
    } else if (dml_direction == "decrease") {  
      # For methylation decreases, avoid already lowly methylated sites
      eligible_sites <- which(mean_p_control > 0.2)
    } else {
      stop("dml_direction must be either 'increase' or 'decrease'")
    }
    
    # Handle case where insufficient eligible sites are available
    if (length(eligible_sites) < n_dml) {
      warning("Number of eligible CpGs (", length(eligible_sites), 
              ") is fewer than requested DMLs (", n_dml, "). Using all eligible sites.")
      dml_sites <- eligible_sites
    } else {
      # Randomly select DML sites from eligible candidates
      dml_sites <- sample(eligible_sites, n_dml)
    }
    
    cat("Injecting DML into", length(dml_sites), "CpG sites with effect size", effect_size, 
        "in", dml_direction, "direction\n")
    
    # Apply methylation shifts to treatment group at selected DML sites
    for (site in dml_sites) {
      # Get original methylation proportions in control group
      p_control <- p_matrix[site, group_A]
      
      # Calculate shifted methylation proportions for treatment group
      if (dml_direction == "increase") {
        p_treatment <- pmin(1, p_control + effect_size)  # Cap at 1.0
      } else {
        p_treatment <- pmax(0, p_control - effect_size)  # Floor at 0.0  
      }
      
      # Update treatment group methylation data
      coverage_treatment <- coverage_matrix[site, group_B]
      meth_matrix[site, group_B] <- rbinom(length(group_B), coverage_treatment, p_treatment)
      p_matrix[site, group_B] <- p_treatment
    }
    
    # Create Ground Truth Table for Injected DML
    
    # Calculate summary statistics for each DML site
    mean_meth_control <- rowMeans(meth_matrix[dml_sites, group_A] / 
                                 pmax(coverage_matrix[dml_sites, group_A], 1))
    mean_meth_treatment <- rowMeans(meth_matrix[dml_sites, group_B] / 
                                   pmax(coverage_matrix[dml_sites, group_B], 1))
    min_cov_control <- apply(coverage_matrix[dml_sites, group_A], 1, min)
    min_cov_treatment <- apply(coverage_matrix[dml_sites, group_B], 1, min)
    
    # Build comprehensive truth table
    truth_df <- data.frame(
      CpG_ID = dml_sites,
      chr = as.character(seqnames(cpg_granges[dml_sites])),
      pos = start(cpg_granges[dml_sites]),
      Effect = effect_size,
      Direction = dml_direction,
      mean_meth_control = mean_meth_control,
      mean_meth_treatment = mean_meth_treatment,
      min_cov_control = min_cov_control,
      min_cov_treatment = min_cov_treatment,
      delta = abs(mean_meth_treatment - mean_meth_control)
    )
    
    # Save ground truth table
    write.csv(truth_df, sub("\\.qs$", "_truth_DML.csv", output_file), row.names = FALSE)
  }
  
  # Construct BSseq Object with Sample Metadata
  
  # Generate systematic sample names
  sample_names <- paste0("Sample_", sprintf("%03d", 1:n_samples))
  
  # Create group assignments (equal numbers in each group)
  group_labels <- c(rep("Control", n_samples / 2), rep("Treatment", n_samples / 2))
  
  # Build sample metadata data frame
  sample_metadata <- data.frame(
    sampleNames = sample_names,
    group = factor(group_labels, levels = c("Control", "Treatment")),
    stringsAsFactors = FALSE
  )
  rownames(sample_metadata) <- sample_names
  
  # Create BSseq object with genomic coordinates, methylation counts, and coverage
  bs_object <- BSseq(gr = cpg_granges,           # GenomicRanges object
                     M = meth_matrix,            # Methylated read counts
                     Cov = coverage_matrix,      # Total coverage
                     sampleNames = sample_names, # Sample identifiers
                     pData = sample_metadata)    # Sample metadata
  
  # Store ground truth information in object metadata for future reference
  metadata(bs_object)$true_K <- K
  metadata(bs_object)$true_component_params <- component_params
  metadata(bs_object)$true_cpg_components <- cpg_components
  metadata(bs_object)$true_dml_sites <- dml_sites
  metadata(bs_object)$simulation_parameters <- list(
    n_cpgs = n_cpgs,
    n_samples = n_samples,
    dml_fraction = dml_fraction,
    effect_size = effect_size,
    mean_coverage = mean_coverage,
    seed = seed
  )
  
  # Save BSseq object in efficient .qs format
  qsave(bs_object, output_file)
  
  # Save True Mixture Model Parameters
  
  # Create summary table of mixture model parameters
  true_model_file <- sub("\\.qs$", "_trueModel.csv", output_file)
  true_model_df <- data.frame(
    component = 1:K,
    alpha = sapply(component_params, function(x) x$alpha),
    beta  = sapply(component_params, function(x) x$beta),
    pi    = mixing_props,  # Mixing proportions
    mean_methylation = sapply(component_params, function(x) x$alpha / (x$alpha + x$beta))
  )
  write.csv(true_model_df, true_model_file, row.names = FALSE)
  
  # Output Summary and Return Object
  
  cat("Synthetic data saved to:", output_file, "\n")
  cat("True DML CpGs saved to:", sub("\\.qs$", "_truth_DML.csv", output_file), "\n")
  
  return(bs_object)
}

# Command-Line Interface Configuration

# Define command-line options with detailed help text
option_list <- list(
  make_option(c("-K", "--components"), 
              type = "integer", default = 3, 
              help = "Number of Beta distribution components in mixture model [default: %default]"),
              
  make_option(c("-n", "--n_cpgs"), 
              type = "integer", default = 100000, 
              help = "Target number of CpG sites in final dataset [default: %default]"),
              
  make_option(c("-s", "--samples"), 
              type = "integer", default = 20, 
              help = "Total number of samples (must be even) [default: %default]"),
              
  make_option("--dml_fraction", 
              type = "double", default = 0.1, 
              help = "Fraction of CpGs designated as DML (0 = null simulation) [default: %default]"),
              
  make_option("--effect_size", 
              type = "double", default = 0.2, 
              help = "Magnitude of methylation difference for DML sites [default: %default]"),
              
  make_option("--true_delta", 
              type = "double", default = 0.1, 
              help = "Minimum effect size threshold for DML classification [default: %default]"),
              
  make_option(c("-f", "--apply_filtering"), 
              type = "character", default = "TRUE", 
              help = "Enable CpG oversampling for post-filtering (TRUE/FALSE) [default: %default]"),
              
  make_option(c("-o", "--output"), 
              type = "character", default = "synthetic_bb_dml.qs", 
              help = "Output filename for BSseq object (.qs extension) [default: %default]"),
              
  make_option(c("-c", "--mean_coverage"), 
              type = "double", default = 20, 
              help = "Expected mean sequencing coverage per CpG site [default: %default]"),
              
  make_option("--mixing_props", 
              type = "character", default = NULL,
              help = "Custom mixture proportions as comma-separated values (e.g., '0.1,0.2,0.7')")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list,
                          description = "Generate synthetic WGBS data with controlled differential methylation",
                          usage = "Usage: %prog [options]\n\nSynthetic WGBS Simulator - Generate realistic bisulfite sequencing data with controlled differential methylation\n\nExample usage:\n  %prog -K 3 -n 100000 -s 20 --dml_fraction 0.1 -o output.qs\n  %prog --components 2 --samples 40 --dml_fraction 0 -o null_sim.qs  # Null simulation")
opt <- parse_args(opt_parser)

# Execute Simulation with Parsed Parameters

cat(">>> Running WGBS simulator with parameters:\n")
apply_filtering_flag <- toupper(opt$apply_filtering) == "TRUE"
mixing_props <- NULL
if (!is.null(opt$mixing_props)) {
  mixing_props <- as.numeric(strsplit(opt$mixing_props, ",")[[1]])
}

generate_synthetic_bb_mixture_with_DML(
  K            = opt$components,
  n_cpgs       = opt$n_cpgs,
  n_samples    = opt$samples,
  dml_fraction = opt$dml_fraction,
  effect_size  = opt$effect_size,
  true_delta   = opt$true_delta,
  apply_filtering = apply_filtering_flag,
  output_file  = opt$output,
  mean_coverage = opt$mean_coverage,
  mixing_props = mixing_props
)

cat(">>> Simulation Completed\n")

