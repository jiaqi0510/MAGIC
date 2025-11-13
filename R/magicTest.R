#!/usr/bin/env Rscript
#
# MAGIC: Methylation Analysis with Genomic Inferred Contexts
#
# Module: Bayesian differential testing interface
# Description: R interface for Bayesian differential methylation testing between
#              conditions. Implements mixture-weighted dispersion estimation and
#              hypothesis testing. Supports case-control comparisons and correlation
#              analysis with multiple testing correction.
#
# Authors: Jiaqi Han, Michael Thompson, Matteo Pellegrini
# Contact: mjthompson69@gmail.com
# License: MIT
# Date: 2025

library(bsseq); library(GenomicRanges); library(optparse); library(matrixStats); library(data.table); library(Rcpp); library(tools); library(qs)

args_list <- commandArgs(trailingOnly = FALSE)
script_arg <- args_list[grep("^--file=", args_list)]

if (length(script_arg) > 0) {
  script_path <- sub("^--file=", "", script_arg[1])
  script_dir <- dirname(normalizePath(script_path))
} else {
  script_dir <- getwd()
}

cpp_file <- file.path(dirname(script_dir), "src", "magicTest.cpp")

if (!file.exists(cpp_file)) {
  stop(sprintf("Cannot find C++ module: %s", cpp_file))
}

cat(sprintf("Loading C++ module: %s\n", cpp_file))
sourceCpp(cpp_file)

print_progress <- function(message, verbose = TRUE) {
  if (verbose) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), message))
}

report_memory <- function(label, verbose = TRUE) {
  if (!verbose) return()
  
  gc_info <- gc()
  used_mb <- sum(gc_info[, "used"]) * (8 * .Machine$sizeof.pointer) / 1024^2
  max_mb <- sum(gc_info[, "max used"]) * (8 * .Machine$sizeof.pointer) / 1024^2
  
  sys_mem <- "N/A"
  if (file.exists("/proc/meminfo")) {
    meminfo <- readLines("/proc/meminfo", n = 3)
    total_kb <- as.numeric(gsub("[^0-9]", "", meminfo[grep("MemTotal", meminfo)]))
    avail_kb <- as.numeric(gsub("[^0-9]", "", meminfo[grep("MemAvailable", meminfo)]))
    if (length(total_kb) > 0 && length(avail_kb) > 0) {
      sys_mem <- sprintf("%.1f/%.1f GB avail", avail_kb/1024^2, total_kb/1024^2)
    }
  }
  
  cat(sprintf("[MEM] %s: R=%.1f MB (max %.1f MB) | System=%s\n", 
              label, used_mb, max_mb, sys_mem))
}

validate_test_data <- function(meth_data, total_data, context = "testing", verbose = TRUE) {
  n_obs_total <- length(meth_data)
  invalid_count <- 0
  
  if (any(!is.finite(meth_data))) {
    n_invalid <- sum(!is.finite(meth_data))
    pct <- 100 * n_invalid / n_obs_total
    if (pct > 1.0 && verbose) {
      warning(sprintf("[VALIDATION] %s: Found %d (%.1f%%) NaN/Inf methylation values",
                     context, n_invalid, pct))
    }
    invalid_count <- invalid_count + n_invalid
  }
  
  if (any(!is.finite(total_data))) {
    n_invalid <- sum(!is.finite(total_data))
    pct <- 100 * n_invalid / n_obs_total
    if (pct > 1.0 && verbose) {
      warning(sprintf("[VALIDATION] %s: Found %d (%.1f%%) NaN/Inf coverage values",
                     context, n_invalid, pct))
    }
    invalid_count <- invalid_count + n_invalid
  }
  
  if (any(meth_data < 0, na.rm = TRUE)) {
    n_invalid <- sum(meth_data < 0, na.rm = TRUE)
    if (n_invalid > 0 && verbose) {
      warning(sprintf("[VALIDATION] %s: Found %d negative methylation counts", context, n_invalid))
    }
    invalid_count <- invalid_count + n_invalid
  }
  
  if (any(total_data <= 0, na.rm = TRUE)) {
    n_invalid <- sum(total_data <= 0, na.rm = TRUE)
    if (n_invalid > 0 && verbose) {
      warning(sprintf("[VALIDATION] %s: Found %d non-positive coverage counts", context, n_invalid))
    }
    invalid_count <- invalid_count + n_invalid
  }
  
  if (any(meth_data > total_data, na.rm = TRUE)) {
    violations <- meth_data > total_data
    violations[is.na(violations)] <- FALSE
    n_violations <- sum(violations)
    if (n_violations > 0 && verbose) {
      warning(sprintf("[VALIDATION] %s: Found %d observations where meth > coverage", 
                     context, n_violations))
    }
    invalid_count <- invalid_count + n_violations
  }
  
  valid_obs <- is.finite(meth_data) & is.finite(total_data) & 
               (meth_data >= 0) & (total_data > 0) & (meth_data <= total_data)
  
  pct_valid <- 100 * sum(valid_obs) / n_obs_total
  
  if (pct_valid < 50) {
    stop(sprintf("[ERROR] %s: Only %.1f%% of observations valid (need >= 50%%)", 
                context, pct_valid))
  }
  
  if (invalid_count > 0 && verbose) {
    print_progress(sprintf("Validation: %.1f%% valid observations", pct_valid), verbose)
  }
  
  return(list(valid = valid_obs, pct_valid = pct_valid))
}


load_mixture_model <- function(mixture_model_dir) {
  model_file <- file.path(mixture_model_dir, "optModel.csv")
  comp_file <- file.path(mixture_model_dir, "methComps.csv")
  
  if (file.exists(model_file)) {
    selected_file <- model_file
    cat("Using optModel.csv (natural space parameters)\n")
  } else if (file.exists(comp_file)) {
    selected_file <- comp_file
    cat("Using methComps.csv (interpretable parameters)\n")
  } else {
    stop("No mixture model file found in ", mixture_model_dir)
  }
  
  model_data <- data.table::fread(selected_file, stringsAsFactors = FALSE)
  
  if (basename(selected_file) == "optModel.csv") {
    params <- list(k_comp = nrow(model_data), pi = model_data$pi, alpha = model_data$alpha, beta = model_data$beta)
  } else {
    total_strength <- (1 / model_data$methDisp) - 1
    params <- list(
      k_comp = nrow(model_data),
      pi = model_data$mixWeight,
      alpha = model_data$meanMeth * total_strength,
      beta = (1 - model_data$meanMeth) * total_strength
    )
    cat("Warning: Converting interpretable parameters to alpha/beta\n")
  }
  
  cat(sprintf("Loaded mixture model: K=%d components\n", params$k_comp))
  return(params)
}

# Filtering functions (simplified from original)
aggressive_sample_filter <- function(bs, args, verbose = TRUE) {
  initial_samps <- ncol(bs)
  samp_data <- as.data.frame(colData(bs))
  if (is.null(rownames(samp_data))) rownames(samp_data) <- sampleNames(bs)
  
  if (!identical(sampleNames(bs), rownames(samp_data))) {
    rownames(samp_data) <- sampleNames(bs)
    colData(bs) <- samp_data
  }
  
  if (!args$group %in% colnames(samp_data)) {
    stop("Group column not found in data: ", args$group)
  }
  
  group_vals <- samp_data[[args$group]]
  valid_group_samps <- which(!is.na(group_vals))
  
  if (length(valid_group_samps) < 2 * args$min_samp) {
    stop("Insufficient samples with valid group values")
  }
  
  bs_filtered <- bs[, valid_group_samps]
  rm(bs); gc()
  
  if (verbose) print_progress(sprintf("Sample filter: %d retained (valid groups)", ncol(bs_filtered)), verbose)
  return(bs_filtered)
}

aggressive_site_filter <- function(bs, args, verbose = TRUE) {
  initial_sites <- nrow(bs)
  cov_data <- getCoverage(bs, type = "Cov")
  samp_data <- as.data.frame(colData(bs))
  group_vals <- samp_data[[args$group]]
  
  unique_groups <- unique(group_vals[!is.na(group_vals)])
  if (length(unique_groups) < 2) stop("Need at least 2 group types for analysis")
  
  group1_type <- unique_groups[1]
  group2_type <- unique_groups[2]
  group1_idx <- which(group_vals == group1_type)
  group2_idx <- which(group_vals == group2_type)
  
  val_matrix <- cov_data >= args$min_reads
  group1_valid_counts <- rowSums(val_matrix[, group1_idx, drop = FALSE])
  group2_valid_counts <- rowSums(val_matrix[, group2_idx, drop = FALSE])
  
  keep_sites <- (group1_valid_counts >= args$min_samp) & (group2_valid_counts >= args$min_samp)
  
  rm(cov_data, val_matrix, group1_valid_counts, group2_valid_counts); gc()
  
  bs_filtered <- bs[keep_sites, ]
  rm(bs); gc()
  
  if (verbose) print_progress(sprintf("Site filter: %d retained (≥%d per group)", sum(keep_sites), args$min_samp), verbose)
  return(bs_filtered)
}

apply_subset_filter <- function(bs, args, verbose = TRUE) {
  if (is.null(args$subset_group) || is.null(args$subset_group_type)) return(bs)
  
  samp_data <- as.data.frame(colData(bs))
  if (!args$subset_group %in% colnames(samp_data)) {
    stop("Subset group not found in data: ", args$subset_group)
  }
  
  subset_vals <- as.character(samp_data[[args$subset_group]])
  match_samps <- which(!is.na(subset_vals) & subset_vals == as.character(args$subset_group_type))
  
  if (length(match_samps) < 2 * args$min_samp) stop("Insufficient samples in subset")
  
  bs <- bs[, match_samps]
  if (verbose) print_progress(sprintf("Subset filter: %d samples retained", ncol(bs)), verbose)
  return(bs)
}

apply_chromosome_filter <- function(bs, args, verbose = TRUE) {
  if (is.null(args$chr)) return(bs)
  
  chr_name <- if (grepl("^chr", args$chr)) args$chr else paste0("chr", args$chr)
  chr_index <- which(as.character(seqnames(bs)) == chr_name)
  
  if (length(chr_index) == 0) stop("No sites found on chromosome: ", chr_name)
  
  bs <- bs[chr_index, ]
  if (verbose) print_progress(sprintf("Chr filter: %d sites on %s", length(chr_index), chr_name), verbose)
  return(bs)
}

apply_data_subset <- function(bs, args, verbose = TRUE) {
  if (is.null(args$data_subset) || args$data_subset == "full") return(bs)
  
  n_feat_initial <- nrow(bs)
  n_feat <- switch(args$data_subset, 
                  "tiny" = min(1000L, n_feat_initial), 
                  "small" = min(10000L, n_feat_initial),
                  "medium" = min(50000L, n_feat_initial),
                  "large" = min(100000L, n_feat_initial),
                  stop("Invalid data subset: ", args$data_subset))
  
  set.seed(42L)
  bs <- bs[sample(n_feat_initial, n_feat), ]
  if (verbose) print_progress(sprintf("Data subset (%s): %d sites retained", args$data_subset, n_feat), verbose)
  return(bs)
}

filter_cpg_comp <- function(bsseq_data, min_frac = 0.8, verbose = TRUE) {
  if (!"cpgs_per_samp" %in% names(mcols(bsseq_data))) {
    if (verbose) print_progress("CpG comp: skipping (no data)", verbose)
    return(bsseq_data)
  }
  
  cpgs_matrix <- mcols(bsseq_data)$cpgs_per_samp
  max_cpgs <- apply(cpgs_matrix, 1, max, na.rm = TRUE)
  min_required <- pmax(1L, as.integer(max_cpgs * min_frac))
  
  val_samp <- matrix(FALSE, nrow = nrow(cpgs_matrix), ncol = ncol(cpgs_matrix))
  for (i in seq_len(nrow(cpgs_matrix))) {
    if (max_cpgs[i] > 0L) val_samp[i, ] <- cpgs_matrix[i, ] >= min_required[i]
  }
  
  keep_tiles <- rowSums(val_samp) > 0L
  filt_data <- bsseq_data[keep_tiles, ]
  mcols(filt_data)$val_samp_mask <- val_samp[keep_tiles, , drop = FALSE]
  mcols(filt_data)$cpgs_per_samp <- NULL
  
  if (verbose) print_progress(sprintf("CpG comp: %d retained (%.1f%% min)", sum(keep_tiles), min_frac * 100), verbose)
  return(filt_data)
}

apply_cpg_filter <- function(bs, args, verbose = TRUE) {
  is_tiled <- "cpgs_per_samp" %in% names(mcols(bs))
  if (!is_tiled) return(list(bs = bs, is_tiled = FALSE))
  
  bs <- filter_cpg_comp(bs, args$min_cpg_fraction, verbose)
  return(list(bs = bs, is_tiled = TRUE))
}

establish_groups <- function(bs, args, verbose = TRUE) {
  samp_data <- as.data.frame(colData(bs))
  if (is.null(rownames(samp_data))) rownames(samp_data) <- sampleNames(bs)
  
  if (!args$group %in% colnames(samp_data)) stop("Group column not found in data")
  
  group_vals <- samp_data[[args$group]]
  vals <- group_vals[!is.na(group_vals)]
  vals_char <- as.character(vals)
  n_unique <- length(unique(vals_char))
  num_vals <- suppressWarnings(as.numeric(vals_char))
  is_numeric <- !all(is.na(num_vals))
  categorical <- !is_numeric || n_unique <= 5 || 
    (all(num_vals == floor(num_vals), na.rm = TRUE) && n_unique <= 10 && 
     min(num_vals, na.rm = TRUE) >= 0 && max(num_vals, na.rm = TRUE) <= 10)
  
  grp_info <- list(categorical = categorical, type = ifelse(categorical, "categorical", "continuous"),
                   is_numeric = is_numeric, values = unique(vals))
  
  if (!categorical) {
    if (verbose) print_progress("Group analysis: continuous (correlation mode)", verbose)
    grp_info$is_correlation <- TRUE
    grp_info$group_values <- num_vals
    grp_info$valid_samps <- sampleNames(bs)
    return(list(bs = bs, grp_info = grp_info))
  }
  
  grp_info$is_correlation <- FALSE
  
  if (!is.null(args$group_types)) {
    group_types <- trimws(strsplit(args$group_types, ",")[[1]])
  } else {
    type_counts <- sort(table(group_vals), decreasing = TRUE)
    group_types <- names(type_counts)[1:2]
  }
  
  grp_info$ref_type <- group_types[1]
  grp_info$test_type <- group_types[2]
  grp_info$is_binary <- length(unique(group_vals)) == 2
  
  current_samps <- sampleNames(bs)
  grp_info$type1_samp <- current_samps[!is.na(group_vals) & group_vals == grp_info$ref_type]
  grp_info$type2_samp <- current_samps[!is.na(group_vals) & group_vals == grp_info$test_type]
  
  samps_to_keep <- c(grp_info$type1_samp, grp_info$type2_samp)
  if (length(samps_to_keep) < length(current_samps)) {
    samp_indices <- match(samps_to_keep, current_samps)
    bs <- bs[, samp_indices]
  }
  
  if (verbose) print_progress(sprintf("Group analysis: %s vs %s (%d + %d samps)", 
                                     grp_info$ref_type, grp_info$test_type, 
                                     length(grp_info$type1_samp), length(grp_info$type2_samp)), verbose)
  
  return(list(bs = bs, grp_info = grp_info))
}

randomize_groups <- function(bs, args, verbose = TRUE) {
  if (!args$randomize) return(bs)
  
  samp_data <- as.data.frame(colData(bs))
  group_vals <- samp_data[[args$group]]
  
  set.seed(as.integer(Sys.time()))
  group_vals <- sample(group_vals)
  samp_data[[args$group]] <- group_vals
  colData(bs) <- samp_data
  
  if (verbose) print_progress("Group randomization applied", verbose)
  return(bs)
}

# Main processing
process_features_chunked <- function(bs, grp_info, args, mixture_model, method_list, global_prior, verbose) {
  chunk_size <- 10000
  n_features <- nrow(bs)
  
  if (verbose) print_progress(sprintf("Processing %d features in chunks of %d", n_features, chunk_size), verbose)
  
  chunk_starts <- seq(1, n_features, chunk_size)
  chunk_ends <- pmin(chunk_starts + chunk_size - 1, n_features)
  n_chunks <- length(chunk_starts)
  
  all_results <- vector("list", n_chunks)
  
  sample_names <- sampleNames(bs)
  samp_data <- as.data.frame(colData(bs))
  group_vals <- samp_data[[args$group]]
  n_samps <- length(sample_names)
  
  if (!grp_info$is_correlation) {
    type1_idx <- which(group_vals == grp_info$ref_type)
    type2_idx <- which(group_vals == grp_info$test_type)
    method_type <- if (length(method_list) == 1 && method_list[1] != "all") method_list[1] else "all"
  } else {
    type1_idx <- integer(0)
    type2_idx <- integer(0)
    method_type <- "correlation"
  }
  
  for (chunk_idx in seq_len(n_chunks)) {
    start_idx <- chunk_starts[chunk_idx]
    end_idx <- chunk_ends[chunk_idx]
    chunk_indices <- start_idx:end_idx
    
    if (verbose && chunk_idx %% 50 == 1) {
      print_progress(sprintf("Processing chunk %d (%d-%d)", chunk_idx, start_idx, end_idx), verbose)
    }
    
    bs_chunk <- bs[chunk_indices, ]
    meth_chunk <- getCoverage(bs_chunk, type = "M")
    total_chunk <- getCoverage(bs_chunk, type = "Cov")
    
    validation <- validate_test_data(meth_chunk, total_chunk, 
                                     context = sprintf("chunk %d", chunk_num),
                                     verbose = (chunk_num == 1))
    
    val_samp_masks <- if ("val_samp_mask" %in% names(mcols(bs_chunk))) {
      mcols(bs_chunk)$val_samp_mask
    } else {
      NULL
    }
    
    n_chunk_features <- nrow(meth_chunk)
    
    val_chunk <- (total_chunk >= args$min_reads)
    if (!is.null(val_samp_masks)) {
      val_chunk <- val_chunk & val_samp_masks
    }
    
    if (grp_info$is_correlation) {
      nSamp1_chunk <- rowSums(val_chunk)
      nSamp2_chunk <- rep(NA_integer_, n_chunk_features)
      avgCov_chunk <- rowMeans(total_chunk, na.rm = FALSE)
    } else {
      nSamp1_chunk <- rowSums(val_chunk[, type1_idx, drop = FALSE])
      nSamp2_chunk <- rowSums(val_chunk[, type2_idx, drop = FALSE])
      avgCov_chunk <- rowMeans(total_chunk, na.rm = FALSE)
    }
    
    if (grp_info$is_correlation) {
      trait_chunk <- matrix(rep(as.numeric(group_vals), n_chunk_features), 
                           nrow = n_chunk_features, ncol = n_samps, byrow = TRUE)
      group_chunk <- matrix(0L, nrow = n_chunk_features, ncol = n_samps)
    } else {
      group_chunk <- matrix(0L, nrow = n_chunk_features, ncol = n_samps)
      trait_chunk <- matrix(0, nrow = n_chunk_features, ncol = n_samps)
      
      group_labels <- integer(n_samps)
      group_labels[type1_idx] <- 0L
      group_labels[type2_idx] <- 1L
      group_chunk <- matrix(rep(group_labels, n_chunk_features), 
                           nrow = n_chunk_features, ncol = n_samps, byrow = TRUE)
    }
    
    chunk_results <- magic_tests_chunked(
      meth_matrix = meth_chunk,
      total_matrix = total_chunk,
      val_matrix = val_chunk,
      group_matrix = group_chunk,
      trait_matrix = trait_chunk,
      pi_universal = mixture_model$pi,
      alpha = mixture_model$alpha,
      beta = mixture_model$beta,
      alpha_0 = rep(1.0, mixture_model$k_comp),
      global_m0 = global_prior$m0,
      global_tau = global_prior$tau,
      test_type = method_type,
      n_threads = args$cores,
      is_correlation = grp_info$is_correlation
    )
    
    chunk_results$nSamp1 <- nSamp1_chunk
    chunk_results$nSamp2 <- nSamp2_chunk
    chunk_results$avgCov <- round(avgCov_chunk, 1)
    
    all_results[[chunk_idx]] <- chunk_results
    
    rm(meth_chunk, total_chunk, val_chunk, group_chunk, trait_chunk, bs_chunk)
    rm(nSamp1_chunk, nSamp2_chunk, avgCov_chunk)
    gc()
  }
  
  if (verbose) print_progress(sprintf("Processed %d chunks, combining results", n_chunks), verbose)
  return(rbindlist(all_results))
}

format_results <- function(cpp_results, bs, grp_info, methods, is_correlation = FALSE) {
  results <- data.table(
    chr = as.character(seqnames(bs)),
    start = start(bs)
  )
  
  n_features <- nrow(results)
  
  if (is_correlation) {
    if ("pearR" %in% names(cpp_results)) {
      results$pearR <- round(cpp_results$pearR, 2)
      results$pearP <- signif(cpp_results$pearP, 3)
      valid_p <- !is.na(results$pearP)
      results$pearQ <- NA_real_
      if (sum(valid_p) > 0) {
        results$pearQ[valid_p] <- signif(p.adjust(results$pearP[valid_p], method = "BH"), 3)
      }
    }
    
    if ("spearRho" %in% names(cpp_results)) {
      results$spearRho <- round(cpp_results$spearRho, 2)
      results$spearP <- signif(cpp_results$spearP, 3)
      valid_p <- !is.na(results$spearP)
      results$spearQ <- NA_real_
      if (sum(valid_p) > 0) {
        results$spearQ[valid_p] <- signif(p.adjust(results$spearP[valid_p], method = "BH"), 3)
      }
    }
    
    if ("magicR" %in% names(cpp_results)) {
      results$magicR <- round(cpp_results$magicR, 2)
      results$magicP <- signif(cpp_results$magicP, 3)
      valid_p <- !is.na(results$magicP)
      results$magicQ <- NA_real_
      if (sum(valid_p) > 0) {
        results$magicQ[valid_p] <- signif(p.adjust(results$magicP[valid_p], method = "BH"), 3)
      }
    }
    
    if ("nSamp1" %in% names(cpp_results)) results$nSamp <- cpp_results$nSamp1
    if ("avgCov" %in% names(cpp_results)) results$avgCov <- cpp_results$avgCov
    
    return(results)
  }
  
  # Two-group analysis formatting
  if (any(c("mixture_mean_group1", "magic_mean_group1") %in% names(cpp_results))) {
    if ("magic_mean_group1" %in% names(cpp_results)) {
      results$meth1 <- round(cpp_results$magic_mean_group1, 3)
      results$meth2 <- round(cpp_results$magic_mean_group2, 3)
    } else if ("mixture_mean_group1" %in% names(cpp_results)) {
      results$meth1 <- round(cpp_results$mixture_mean_group1, 3)
      results$meth2 <- round(cpp_results$mixture_mean_group2, 3)
    }
  }
  
  if ("nSamp1" %in% names(cpp_results)) results$nSamp1 <- cpp_results$nSamp1
  if ("nSamp2" %in% names(cpp_results)) results$nSamp2 <- cpp_results$nSamp2
  if ("avgCov" %in% names(cpp_results)) results$avgCov <- cpp_results$avgCov
  
  if ("mixture_mean_group1" %in% names(cpp_results)) {
    results$mixMeth1 <- round(cpp_results$mixture_mean_group1, 3)
    results$mixMeth2 <- round(cpp_results$mixture_mean_group2, 3)
  }
  
  if ("magic_dom_comp1" %in% names(cpp_results)) {
    results$domComp1 <- as.integer(cpp_results$magic_dom_comp1)
    results$domComp2 <- as.integer(cpp_results$magic_dom_comp2)
    
    ent1_vals <- as.numeric(cpp_results$magic_entropy1)
    ent2_vals <- as.numeric(cpp_results$magic_entropy2)
    
    results$ent1 <- ifelse(is.finite(ent1_vals), round(ent1_vals, 3), NA_real_)
    results$ent2 <- ifelse(is.finite(ent2_vals), round(ent2_vals, 3), NA_real_)
  }
  
  if ("mixture_p_value" %in% names(cpp_results)) {
    results$mixtureDiff <- round(cpp_results$mixture_difference, 3)
    results$mixtureWald <- round(cpp_results$mixture_wald_stat, 2)
    results$mixtureP <- signif(cpp_results$mixture_p_value, 3)
    
    valid_p <- !is.na(results$mixtureP)
    results$mixtureQ <- NA_real_
    if (sum(valid_p) > 0) {
      results$mixtureQ[valid_p] <- signif(p.adjust(results$mixtureP[valid_p], method = "BH"), 3)
    }
  }
  
  if ("magic_bf_10" %in% names(cpp_results)) {
    results$bfBF <- round(cpp_results$magic_bf_10, 2)
    results$bfProb <- signif(cpp_results$magic_prob, 3)
  }
  
  return(results)
}

summarize_results <- function(results, methods, mixture_model, verbose, is_correlation = FALSE) {
  if (!verbose) return()
  
  if (is_correlation) {
    if ("pearP" %in% names(results)) {
      n_sig <- sum(results$pearQ < 0.05, na.rm = TRUE)
      print_progress(sprintf("Pearson significant (FDR < 0.05): %d (%.1f%%)", 
                            n_sig, (n_sig / nrow(results)) * 100), verbose)
    }
    
    if ("spearP" %in% names(results)) {
      n_sig <- sum(results$spearQ < 0.05, na.rm = TRUE)
      print_progress(sprintf("Spearman significant (FDR < 0.05): %d (%.1f%%)", 
                            n_sig, (n_sig / nrow(results)) * 100), verbose)
    }
    
    if ("magicP" %in% names(results)) {
      n_sig <- sum(results$magicQ < 0.05, na.rm = TRUE)
      print_progress(sprintf("MAGIC correlation significant (FDR < 0.05): %d (%.1f%%)", 
                            n_sig, (n_sig / nrow(results)) * 100), verbose)
    }
    return()
  }
  
  # Two-group analysis summaries
  if ("bfProb" %in% names(results)) {
    n_sig <- sum(results$bfProb > 0.5, na.rm = TRUE)
    print_progress(sprintf("BF significant (P(diff) > 0.5): %d (%.1f%%)", 
                          n_sig, (n_sig / nrow(results)) * 100), verbose)
  }
  
  if ("mixtureP" %in% names(results)) {
    n_sig <- sum(results$mixtureQ < 0.05, na.rm = TRUE)
    print_progress(sprintf("Mixture test significant (FDR < 0.05): %d (%.1f%%)", 
                          n_sig, (n_sig / nrow(results)) * 100), verbose)
    
    if ("ent1" %in% names(results)) {
      mean_entropy <- mean(c(results$ent1, results$ent2), na.rm = TRUE)
      max_entropy <- log(mixture_model$k_comp)
      print_progress(sprintf("Average posterior entropy: %.3f (max: %.3f)", mean_entropy, max_entropy), verbose)
    }
  }
}

build_dir_name <- function(grp_info, args, is_tiled, methods, mixture_model) {
  method_abbrev <- c("mixture" = "mix", "magic" = "bf", "all" = "all")
  
  dir_parts <- c("magicTests")
  dir_parts <- c(dir_parts, paste0("K", mixture_model$k_comp))
  
  if (length(methods) == 1 && methods[1] == "all") {
    dir_parts <- c(dir_parts, "all")
  } else {
    abbrevs <- method_abbrev[methods]
    sorted_abbrevs <- sort(abbrevs)
    dir_parts <- c(dir_parts, paste(sorted_abbrevs, collapse = ""))
  }
  
  group_camel <- gsub("_([a-z])", "\\U\\1", args$group, perl = TRUE)
  
  if (grp_info$is_binary) {
    dir_parts <- c(dir_parts, group_camel)
  } else {
    clean_ref <- gsub("[^a-zA-Z0-9]", "", tools::toTitleCase(grp_info$ref_type))
    clean_test <- gsub("[^a-zA-Z0-9]", "", tools::toTitleCase(grp_info$test_type))
    comparison <- paste0(group_camel, clean_ref, "v", clean_test)
    dir_parts <- c(dir_parts, comparison)
  }
  
  if (!is.null(args$subset_group) && !is.null(args$subset_group_type)) {
    clean_type <- gsub("[^a-zA-Z0-9]", "", tools::toTitleCase(args$subset_group_type))
    dir_parts <- c(dir_parts, paste0("sub", clean_type))
  }
  
  dir_parts <- c(dir_parts, sprintf("r%dm%d", args$min_reads, args$min_samp))
  
  if (is_tiled) {
    cpg_frac_pct <- round(args$min_cpg_fraction * 100)
    dir_parts <- c(dir_parts, sprintf("cpgfrac%d", cpg_frac_pct))
  }
  
  if (!is.null(args$data_subset) && args$data_subset != "full") {
    dir_parts <- c(dir_parts, args$data_subset)
  }
  
  if (!is.null(args$chr)) {
    dir_parts <- c(dir_parts, paste0("chr", gsub("[^a-zA-Z0-9]", "", args$chr)))
  }
  
  dir_name <- paste(dir_parts, collapse = "_")
  
  if (dir.exists(dir_name)) {
    dir_name <- paste0(dir_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  return(dir_name)
}

save_results_and_summary <- function(results, mixture_model, analysis_type, args, grp_info, 
                                     is_tiled, data_type, run_time, methods, verbose = TRUE) {
  output_dir_name <- build_dir_name(grp_info, args, is_tiled, methods, mixture_model)
  dir.create(output_dir_name, showWarnings = FALSE, recursive = TRUE)

  qs::qsave(results, file.path(output_dir_name, "magic_results.qs"))
  if (verbose) print_progress(sprintf("Results saved: %s", basename(output_dir_name)), verbose)

  method_descriptions <- list(
    magic = "- BF (Bayes Factor): Bayesian mixture component analysis\n  • Uses Dirichlet-multinomial marginal likelihoods\n  • Evidence quantified via posterior model probabilities", 
    mixture = "- Mixture Test: Law of Total Variance-based differential methylation\n  • Accounts for uncertainty in component assignments\n  • Uses mixture-weighted expectations with within+between component variance"
  )

  summary_lines <- c(
    "MAGIC Analysis Summary",
    paste0(rep("=", 24L), collapse = ""),
    paste("Date:", Sys.time()),
    paste("Analysis:", sprintf("MAGIC (K=%d, Methods=%s)", mixture_model$k_comp, paste(methods, collapse = ","))),
    paste("Statistical methods:", paste(methods, collapse = ", ")),
    paste("Method: Mixture-based Analysis with Bayesian Inference"),
    paste("Mixture model dir:", args$mixture_model),
    paste("Mixture components:", mixture_model$k_comp),
    paste("Input data:", args$input),
    paste("Data type:", data_type),
    paste("Features analyzed:", format(nrow(results), big.mark = ",")),
    paste("Runtime:", sprintf("%.1f min", run_time)),
    paste("Cov threshold:", args$min_reads, "reads"),
    paste("Min samps per group type:", args$min_samp)
  )

  optional_filters <- c()
  if (is_tiled) optional_filters <- c(optional_filters, paste("CpG composition filter:", sprintf("%.0f%% min fraction", args$min_cpg_fraction * 100)))
  if (!is.null(args$subset_group) && !is.null(args$subset_group_type)) optional_filters <- c(optional_filters, paste("Samp subset:", sprintf("%s == %s", args$subset_group, args$subset_group_type)))
  if (!is.null(args$chr)) optional_filters <- c(optional_filters, paste("Chromosome filter:", args$chr))
  if (!is.null(args$data_subset) && args$data_subset != "full") optional_filters <- c(optional_filters, paste("Data subset:", args$data_subset))
  
  summary_lines <- c(summary_lines, optional_filters)
  
  if (!grp_info$is_correlation) {
    summary_lines <- c(summary_lines, paste("Group types compared:", sprintf("%s vs %s", grp_info$ref_type, grp_info$test_type)))
  }

  if (nrow(results) > 0) {
    included_methods <- intersect(methods, names(method_descriptions))
    if ("all" %in% methods) included_methods <- names(method_descriptions)
    
    if (length(included_methods) > 0) {
      summary_lines <- c(summary_lines, "", "Statistical Method Details:")
      for (method in included_methods) {
        summary_lines <- c(summary_lines, method_descriptions[[method]])
      }
    }
  }

  writeLines(summary_lines, file.path(output_dir_name, "magic_summary.txt"))
}

option_list <- list(
  make_option(c("-i", "--input"), type = "character", 
              help = "BSseq input file (.qs format) [required]"),
  make_option(c("-g", "--group"), type = "character", 
              help = "Group column name in patient data [required]"),
  make_option(c("-c", "--compare"), type = "character", 
              help = "Group types to compare (comma-separated, e.g. 'case,control')"),
  make_option(c("-s", "--subset"), type = "character", 
              help = "Subset group column name for patient filtering"),
  make_option(c("-k", "--keep"), type = "character", 
              help = "Subset group type to retain (filters to only this type)"),
  make_option(c("-r", "--min_reads"), type = "integer", default = 5L,
              help = "Minimum read coverage per site [default: %default]"),
  make_option(c("-m", "--min_samp"), type = "integer", default = 3L,
              help = "Minimum samples per group type [default: %default]"),
  make_option(c("-n", "--cores"), type = "integer", default = 1L,
              help = "Number of CPU cores [default: %default]"),
  make_option(c("--min_var"), type = "numeric", default = 0.01,
              help = "Minimum variance threshold [default: %default]"),
  make_option(c("--data_subset"), type = "character", default = "full",
              help = "Data subset size. Options: tiny, small, medium, large, full [default: %default]"),
  make_option(c("-f", "--min_cpg_fraction"), type = "numeric", default = 0.8,
              help = "Minimum CpG fraction for tiled data [default: %default]"),
  make_option(c("--chr"), type = "character",
              help = "Chromosome to analyze (e.g. 'chr1' or '1')"),
  make_option(c("-z", "--randomize"), action = "store_true", default = FALSE,
              help = "Randomize group assignments (for null testing)"),
  make_option(c("-q", "--quiet"), action = "store_true", default = FALSE,
              help = "Suppress progress messages"),
  make_option(c("--mixture_model"), type = "character", 
              help = "Path to mixture model directory [required]"),
  make_option(c("--prior_strength"), type = "numeric", default = 1.0,
              help = "Dirichlet prior concentration parameter [default: %default]"),
  make_option(c("--methods"), type = "character", default = "all",
              help = "Statistical methods: mixture, magic, all [default: %default]")
)

main <- function() {
  args <- parse_args(OptionParser(option_list = option_list))

  if (is.null(args$input)) stop("Error: --input parameter is required")
  if (is.null(args$group)) stop("Error: --group parameter is required") 
  if (is.null(args$mixture_model)) stop("Error: --mixture_model parameter is required")
  if (!file.exists(args$input)) stop("Error: input file does not exist: ", args$input)
  if (!dir.exists(args$mixture_model)) stop("Error: mixture model directory does not exist: ", args$mixture_model)
  
  # Map old parameter names for backward compatibility
  if (is.null(args$subset) && !is.null(args$subset_group)) {
    args$subset <- args$subset_group
  }
  if (is.null(args$keep) && !is.null(args$subset_group_type)) {
    args$keep <- args$subset_group_type
  }

  args$subset_group <- args$subset
  args$subset_group_type <- args$keep
  args$group_types <- args$compare
  
  if (!is.null(args$randomize) && is.na(args$randomize)) args$randomize <- FALSE

  start_time <- Sys.time()
  verbose <- !args$quiet

  method_list <- trimws(strsplit(args$methods, ",")[[1]])
  
  cat("MAGIC - Mixture-based Analysis with Bayesian Inference v3.0 SIMPLIFIED\n")
  cat(paste0(rep("=", 60L), collapse = ""), "\n")

  print_progress("Loading mixture model", verbose)
  mixture_model <- load_mixture_model(args$mixture_model)

  print_progress(sprintf("Config: %s | %d min cov | %d cores | MAGIC K=%d | Methods: %s", 
                        basename(args$input), args$min_reads, args$cores, 
                        mixture_model$k_comp, paste(method_list, collapse = ",")), verbose)

  # 1. LOAD DATA
  print_progress("Loading BSseq data", verbose)
  bsseq_data <- qs::qread(args$input)

  # 2. AGGRESSIVE CORE FILTERING
  print_progress("Applying aggressive core filtering", verbose)
  
  bs <- aggressive_sample_filter(bsseq_data, args, verbose)
  rm(bsseq_data); gc()
  
  bs <- aggressive_site_filter(bs, args, verbose)
  
  cpg_result <- apply_cpg_filter(bs, args, verbose)
  bs <- cpg_result$bs
  is_tiled <- cpg_result$is_tiled

  # 3. OPTIONAL USER-REQUESTED SUBSETTING
  bs <- apply_subset_filter(bs, args, verbose)
  bs <- apply_chromosome_filter(bs, args, verbose)
  bs <- apply_data_subset(bs, args, verbose)

  # 4. FINAL PROCESSING
  group_result <- establish_groups(bs, args, verbose)
  bs <- group_result$bs
  grp_info <- group_result$grp_info

  if (!grp_info$is_correlation) {
    if (length(grp_info$type1_samp) < args$min_samp || length(grp_info$type2_samp) < args$min_samp) {
      if (verbose) {
        print_progress("Analysis terminated: insufficient samps after filtering", verbose)
        print_progress(sprintf("Group type 1 (%s): %d samps (need %d)", 
                              grp_info$ref_type, length(grp_info$type1_samp), args$min_samp), verbose)
        print_progress(sprintf("Group type 2 (%s): %d samps (need %d)", 
                              grp_info$test_type, length(grp_info$type2_samp), args$min_samp), verbose)
      }
      return(invisible(NULL))
    }
  }

  bs <- randomize_groups(bs, args, verbose)

  data_type <- if (is_tiled) "tiles" else "CpG sites"
  print_progress(sprintf("Final dataset: %d %s × %d samples", nrow(bs), data_type, ncol(bs)), verbose)

  # 5. GLOBAL PRIOR CALCULATION (if needed for mixture test)
  global_prior <- list(m0 = -3.0, tau = 1.0)
  if (!grp_info$is_correlation && any(c("mixture", "all") %in% method_list)) {
    print_progress("Calculating global dispersion prior", verbose)
    
    meth_matrix_full <- getCoverage(bs, type = "M")
    total_matrix_full <- getCoverage(bs, type = "Cov") 
    
    global_prior <- estimate_global_prior_from_matrix(meth_matrix_full, total_matrix_full)
    
    print_progress(sprintf("Global prior: m0=%.3f, tau=%.3f", global_prior$m0, global_prior$tau), verbose)
    
    rm(meth_matrix_full, total_matrix_full); gc()
  }

  # 6. CHUNKED ANALYSIS
  print_progress("Starting chunked matrix-based analysis", verbose)
  
  cpp_results <- process_features_chunked(bs, grp_info, args, mixture_model, method_list, 
                                         global_prior, verbose)

  results <- format_results(cpp_results, bs, grp_info, method_list, 
                           is_correlation = grp_info$is_correlation)

  if (grp_info$is_correlation) {
    valid_cols <- intersect(c("pearP", "spearP", "magicP"), names(results))
  } else {
    valid_cols <- intersect(c("bfProb", "mixtureP"), names(results))
  }
  
  if (length(valid_cols) > 0) {
    valid_results <- !is.na(results[[valid_cols[1]]])
    results <- results[valid_results]
  }

  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

  if (grp_info$is_correlation) {
    print_progress(sprintf("Analysis complete: %d %s | MAGIC correlation (K=%d) | %.1f min", 
                   nrow(results), data_type, mixture_model$k_comp, run_time), verbose)
  } else {
    print_progress(sprintf("Analysis complete: %d %s | MAGIC (K=%d, %s) | %.1f min", 
                   nrow(results), data_type, mixture_model$k_comp, paste(method_list, collapse = ","), run_time), verbose)
  }

  summarize_results(results, method_list, mixture_model, verbose, is_correlation = grp_info$is_correlation)

  save_results_and_summary(results, mixture_model, "cat", args, grp_info, 
                           is_tiled, data_type, run_time, method_list, verbose)

  print_progress("MAGIC analysis completed successfully", verbose)
}

if (!interactive()) {
  main()
}
