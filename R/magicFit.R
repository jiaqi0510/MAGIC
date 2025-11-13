#!/usr/bin/env Rscript
#
# MAGIC: Methylation Analysis with Genomic Inferred Contexts
#
# Module: Mixture model optimization interface
# Description: R interface for beta-binomial mixture model optimization. Handles
#              data loading, filtering, validation, parameter initialization, and
#              convergence diagnostics. Provides command-line interface for model
#              fitting with multiple random initializations.
#
# Authors: Jiaqi Han, Michael Thompson, Matteo Pellegrini
# Contact: mjthompson69@gmail.com
# License: MIT
# Date: 2025

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(matrixStats)
  library(data.table)
  library(optparse)
  library(Rcpp)
})

EPSILON <- 1e-10
DEFAULT_HOLDOUT_SEED <- 12345
DEFAULT_CONV_SEED <- 54321

if (!exists("compute_ll")) {
  args_list <- commandArgs(trailingOnly = FALSE)
  script_arg <- args_list[grep("^--file=", args_list)]
  
  if (length(script_arg) > 0) {
    script_path <- sub("^--file=", "", script_arg[1])
    script_dir <- dirname(normalizePath(script_path))
  } else {
    script_dir <- getwd()
  }
  
  cpp_file <- file.path(dirname(script_dir), "src", "magicFit.cpp")
  
  if (!file.exists(cpp_file)) {
    stop(sprintf("Cannot find C++ module: %s", cpp_file))
  }
  
  Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
  Sys.setenv("PKG_LIBS" = "-fopenmp")
  sourceCpp(cpp_file)
}

required_cpp_functions <- c("compute_ll", "compute_gradients")
missing_functions <- required_cpp_functions[!sapply(required_cpp_functions, exists)]
if (length(missing_functions) > 0) {
  stop(sprintf("Missing required C++ functions: %s", paste(missing_functions, collapse = ", ")))
}

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "BSseq QS file"),
  make_option(c("-k", "--k_components"), type = "integer", help = "Number of mixture components"),
  make_option(c("--max_iters"), type = "integer", default = 500L),
  make_option(c("-t", "--tol"), type = "numeric", default = 1e-6),
  make_option(c("-r", "--min_reads"), type = "integer", default = 5L),
  make_option(c("-m", "--min_samples"), type = "integer", default = 10L),
  make_option(c("-v", "--min_var"), type = "numeric", default = 0.001),
  make_option(c("--chr"), type = "character", default = NULL),
  make_option(c("-q", "--quiet"), action = "store_true", default = FALSE),
  make_option(c("--filter_sample_outliers"), action = "store_true", default = FALSE),
  make_option(c("--sample_thresh"), type = "numeric", default = 3.0),
  make_option(c("--threads"), type = "integer", default = NULL),
  make_option(c("--no_sex_chr"), action = "store_true", default = FALSE),
  make_option(c("--synthetic"), action = "store_true", default = FALSE, help = "Skip filtering for synthetic data"),
  make_option(c("--subset_trait"), type = "character", default = NULL),
  make_option(c("--subset_trait_val"), type = "character", default = NULL),
  make_option(c("--n_conv_runs"), type = "integer", default = 1),
  make_option(c("--holdout"), action = "store_true", default = FALSE),
  make_option(c("--holdout_frac"), type = "numeric", default = 0.2),
  make_option(c("--conv_seed"), type = "integer", default = DEFAULT_CONV_SEED),
  make_option(c("--seed_methylation_file"), type = "character", default = NULL)
)

args <- parse_args(OptionParser(option_list = option_list))
if (is.null(args)) quit(status = 0)

if (!is.null(args$threads)) {
  Sys.setenv("OMP_NUM_THREADS" = as.character(args$threads))
  cat(sprintf("Set OpenMP threads to %d\n", args$threads))
} else {
  n_cores <- NULL
  if (!is.na(Sys.getenv("OMP_NUM_THREADS", unset = NA))) {
    n_cores <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
  } else {
    n_cores <- min(parallel::detectCores(), 8)
  }
  Sys.setenv("OMP_NUM_THREADS" = as.character(n_cores))
  cat(sprintf("Using %d OpenMP threads\n", n_cores))
}

get_omp_threads <- function() {
  omp_threads <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  if (!is.na(omp_threads)) return(as.integer(omp_threads))
  return(NA)
}

if (is.null(args$input) || is.null(args$k_components)) {
  stop("Required arguments missing: input, k_components")
}
if (!file.exists(args$input)) stop(sprintf("Input file does not exist: %s", args$input))
if (args$k_components <= 0) stop("k_components must be positive")

natural_to_log <- function(alpha, beta, pi) {
  k_comp <- length(alpha)
  if (any(alpha <= 0) || any(beta <= 0)) {
    stop("Alpha and beta params must be positive")
  }
  if (any(pi <= 0) || abs(sum(pi) - 1.0) > 1e-6) {
    stop("Pi params must be positive and sum to 1")
  }
  
  log_alpha <- log(alpha)
  log_beta <- log(beta) 
  
  params <- c(log_alpha, log_beta)
  
  if (k_comp > 1) {
    log_pi <- log(pi[1:(k_comp-1)])
    params <- c(params, log_pi)
  }
  
  return(params)
}

log_to_natural <- function(log_params, k_comp) {
  exp_size <- if (k_comp == 1) 2 else (3 * k_comp - 1)
  
  if (length(log_params) != exp_size) {
    stop(sprintf("Parameter vector size mismatch: expected %d, got %d", exp_size, length(log_params)))
  }
  
  alpha <- exp(log_params[1:k_comp])
  beta <- exp(log_params[(k_comp+1):(2*k_comp)])
  
  if (k_comp > 1) {
    offset <- 2 * k_comp
    log_pi_free <- log_params[(offset+1):(offset+k_comp-1)]
    sum_exp <- 1.0 + sum(exp(log_pi_free))
    pi <- numeric(k_comp)
    pi[1:(k_comp-1)] <- exp(log_pi_free) / sum_exp
    pi[k_comp] <- 1.0 / sum_exp
  } else {
    pi <- 1.0
  }
  
  return(list(alpha = alpha, beta = beta, pi = pi))
}

read_seed <- function(seed_file, k_comp) {
  if (is.null(seed_file) || !file.exists(seed_file)) return(NULL)
  
  seed_data <- read.csv(seed_file, stringsAsFactors = FALSE)
  
  required_cols <- c("mixWeight", "meanMeth", "methDisp")
  if ("alpha" %in% names(seed_data) && "beta" %in% names(seed_data) && "pi" %in% names(seed_data)) {
    if (nrow(seed_data) != k_comp) {
      stop(sprintf("Seed file has %d components, but K=%d requested", nrow(seed_data), k_comp))
    }
    alpha <- seed_data$alpha
    beta <- seed_data$beta
    pi <- seed_data$pi
  } else if (all(required_cols %in% names(seed_data))) {
    if (nrow(seed_data) != k_comp) {
      stop(sprintf("Seed file has %d components, but K=%d requested", nrow(seed_data), k_comp))
    }
    mu <- seed_data$meanMeth
    phi <- seed_data$methDisp
    if (any(mu <= 0 | mu >= 1)) stop("meanMeth values must be between 0 and 1 (exclusive)")
    if (any(phi <= 0 | phi >= 1)) stop("methDisp values must be between 0 and 1 (exclusive)")
    
    alpha_plus_beta <- (1 - phi) / phi
    alpha <- mu * alpha_plus_beta
    beta <- (1 - mu) * alpha_plus_beta
    pi <- seed_data$mixWeight / sum(seed_data$mixWeight)
  } else {
    stop(sprintf("Seed file missing required columns. Need either (alpha, beta, pi) or (%s)", 
                 paste(required_cols, collapse = ", ")))
  }
  
  if (any(alpha <= 0 | beta <= 0 | !is.finite(alpha) | !is.finite(beta))) {
    stop("Invalid alpha/beta params in seed file")
  }
  
  cat(sprintf("Loaded seed params: K=%d\n", k_comp))
  return(list(alpha = alpha, beta = beta, pi = pi))
}

sort_params <- function(params) {
  k_comp <- length(params$alpha)
  if (k_comp == 1) return(params)
  
  mean_meth <- params$alpha / (params$alpha + params$beta)
  sort_order <- order(mean_meth)
  
  sorted_params <- list(
    alpha = params$alpha[sort_order],
    beta = params$beta[sort_order],
    pi = params$pi[sort_order]
  )
  
  return(sorted_params)
}

ab_to_mu_phi <- function(alpha, beta) {
  if (any(alpha <= 0 | beta <= 0 | !is.finite(alpha) | !is.finite(beta))) {
    stop("Invalid alpha or beta params")
  }
  sum_ab <- alpha + beta
  mu <- alpha / sum_ab
  phi <- 1.0 / (sum_ab + 1.0)
  return(c(mu, phi))
}

objective <- function(log_params, meth_mat, total_mat, valid_mask, k_comp, n_feat, n_samples) {
  feature_ll <- compute_ll(meth_mat, total_mat, valid_mask, log_params)
  total_ll <- sum(feature_ll[is.finite(feature_ll)])
  return(-total_ll)
}

gradient <- function(log_params, meth_mat, total_mat, valid_mask, k_comp, n_feat, n_samples) {
  grad_log_params <- compute_gradients(meth_mat, total_mat, valid_mask, log_params)
  return(-grad_log_params)
}

set_bounds <- function(k_comp) {
  lower_bounds <- rep(-10.0, 2 * k_comp)   
  upper_bounds <- rep(10.0, 2 * k_comp)   
  
  if (k_comp > 1) {
    lower_bounds <- c(lower_bounds, rep(-10.0, k_comp - 1))  
    upper_bounds <- c(upper_bounds, rep(5.0, k_comp - 1))     
  }
  
  exp_size <- if (k_comp == 1) 2 else (3 * k_comp - 1)
  
  if (length(lower_bounds) != exp_size || length(upper_bounds) != exp_size) {
    stop(sprintf("Bounds size mismatch: expected %d, got lower=%d, upper=%d", 
         exp_size, length(lower_bounds), length(upper_bounds)))
  }
  
  return(list(lower = lower_bounds, upper = upper_bounds))
}

init_params <- function(meth_mat, total_mat, valid_mask, k_comp, seed = NULL) {
  if (!is.null(seed)) {
    alpha_init <- seed$alpha
    beta_init <- seed$beta
    pi_init <- seed$pi
  } else {
    feat_meth_counts <- rowSums(meth_mat * valid_mask)
    feat_total_counts <- rowSums(total_mat * valid_mask)
    valid_features <- feat_total_counts > 0
    feat_props <- ifelse(valid_features, feat_meth_counts / feat_total_counts, 0.5)
    
    valid_props <- feat_props[feat_props > 0 & feat_props < 1 & is.finite(feat_props)]
    if (length(valid_props) == 0) stop("No valid methylation proportions found")
    
    alpha_init <- numeric(k_comp)
    beta_init <- numeric(k_comp)
    
    if (k_comp == 1) {
      init_mu <- median(valid_props)
      concentration <- 9
      alpha_init[1] <- init_mu * concentration
      beta_init[1] <- (1 - init_mu) * concentration
    } else {
      quantile_points <- seq(0.1, 0.9, length.out = k_comp)
      init_mus <- quantile(valid_props, quantile_points, na.rm = TRUE)
      init_mus <- sort(init_mus)
      
      for (k in 2:k_comp) {
        if (init_mus[k] - init_mus[k-1] < 0.05) {
          init_mus[k] <- init_mus[k-1] + 0.05
        }
      }
      init_mus <- pmax(0.05, pmin(0.95, init_mus))
      
      base_concentration <- 8
      for (k in 1:k_comp) {
        concentration <- base_concentration + (k-1) * 0.5
        alpha_init[k] <- init_mus[k] * concentration
        beta_init[k] <- (1 - init_mus[k]) * concentration
      }
    }
    
    pi_init <- rep(1/k_comp, k_comp)
  }
  
  log_params <- natural_to_log(alpha_init, beta_init, pi_init)
  return(log_params)
}

check_parameter_boundaries <- function(params, k_comp, context = "optimization") {
  MIN_PARAM <- 0.01
  MAX_PARAM <- 1000.0
  MIN_PI <- 0.001
  MAX_PI <- 0.999
  
  warnings_issued <- FALSE
  
  for (k in 1:k_comp) {
    if (params$alpha[k] < MIN_PARAM + 0.01) {
      warning(sprintf("[%s] Alpha parameter for component %d near lower boundary (%.3f)\n  Suggestion: Component may be poorly constrained - consider reducing K or adding data",
                     context, k, params$alpha[k]))
      warnings_issued <- TRUE
    } else if (params$alpha[k] > MAX_PARAM - 1) {
      warning(sprintf("[%s] Alpha parameter for component %d near upper boundary (%.1f)\n  Suggestion: Component may be poorly constrained - consider reducing K or adding data",
                     context, k, params$alpha[k]))
      warnings_issued <- TRUE
    }
    
    if (params$beta[k] < MIN_PARAM + 0.01) {
      warning(sprintf("[%s] Beta parameter for component %d near lower boundary (%.3f)\n  Suggestion: Component may be poorly constrained - consider reducing K or adding data",
                     context, k, params$beta[k]))
      warnings_issued <- TRUE
    } else if (params$beta[k] > MAX_PARAM - 1) {
      warning(sprintf("[%s] Beta parameter for component %d near upper boundary (%.1f)\n  Suggestion: Component may be poorly constrained - consider reducing K or adding data",
                     context, k, params$beta[k]))
      warnings_issued <- TRUE
    }
    
    if (params$pi[k] < MIN_PI + 0.001) {
      warning(sprintf("[%s] Mixing proportion for component %d near lower boundary (%.4f)\n  Suggestion: Component may be unidentifiable - consider reducing K",
                     context, k, params$pi[k]))
      warnings_issued <- TRUE
    } else if (params$pi[k] > MAX_PI - 0.001) {
      warning(sprintf("[%s] Mixing proportion for component %d near upper boundary (%.4f)\n  Suggestion: May indicate single dominant component - consider reducing K",
                     context, k, params$pi[k]))
      warnings_issued <- TRUE
    }
  }
  
  pi_sum <- sum(params$pi)
  if (abs(pi_sum - 1.0) > 1e-4) {
    warning(sprintf("[%s] Mixing proportions sum to %.6f (should be 1.0)\n  This indicates a numerical issue in parameter transformation",
                   context, pi_sum))
  }
  
  return(!warnings_issued)
}

optimize_model <- function(meth_mat, total_mat, k_comp, max_iters, tol, min_reads = 5, seed = NULL) {
  n_feat <- nrow(meth_mat)
  n_samples <- ncol(meth_mat)
  
  if (n_feat <= 0 || n_samples <= 0 || k_comp <= 0) stop("Invalid dimensions or k_comp")
  
  valid_mask <- total_mat >= min_reads
  
  omp_threads <- get_omp_threads()
  thread_info <- if (!is.na(omp_threads)) paste0(omp_threads, " threads") else "threads unknown"
  
  cat(sprintf("Optimizing MAGIC K=%d (%s)\n", k_comp, thread_info))
  cat(sprintf("Data: %d features × %d samples, %d valid observations\n", 
              n_feat, n_samples, sum(valid_mask)))
  
  bounds <- set_bounds(k_comp)
  log_params_init <- init_params(meth_mat, total_mat, valid_mask, k_comp, seed)
  
  start_time <- Sys.time()
  
  control <- list(
    maxit = max_iters,
    factr = tol / .Machine$double.eps,
    pgtol = 1e-6,
    lmm = 15,
    trace = 1,
    REPORT = 25
  )
  
  result <- optim(par = log_params_init,
        fn = objective,
        gr = gradient,
        method = "L-BFGS-B",
        lower = bounds$lower,
        upper = bounds$upper,
        control = control,
        meth_mat = meth_mat,
        total_mat = total_mat,
        valid_mask = valid_mask,
        k_comp = k_comp,
        n_feat = n_feat,
        n_samples = n_samples)
  
  runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  final_conv <- (result$convergence == 0 || result$convergence == 1 || result$convergence %in% c(51, 52))
  
  final_log_params <- result$par
  final_params <- log_to_natural(final_log_params, k_comp)
  final_params <- sort_params(final_params)
  
  check_parameter_boundaries(final_params, k_comp, context = "OPTIMIZATION")
  
  if (result$convergence != 0 && result$convergence != 1) {
    warning(sprintf("[OPTIMIZATION] Non-standard convergence code: %d\n  Suggestion: Results may be unreliable - try increasing iterations or checking data quality",
                   result$convergence))
  }
  
  data_ll <- -result$value
  
  if (!is.finite(data_ll)) {
    stop(sprintf("[ERROR] Log-likelihood is not finite (%.6f)\n  This indicates numerical instability - check data quality and parameter initialization",
                data_ll))
  }
  
  n_params <- if (k_comp == 1) 2 else (3 * k_comp - 1)
  
  n_obs <- n_feat
  bic <- -2 * data_ll + n_params * log(n_obs)
  aic <- -2 * data_ll + 2 * n_params
  
  cat(sprintf("Conv in %d iters, %.1fs (code: %d)\n", result$counts[1], runtime, result$convergence))
  
  return(list(
    k = k_comp,
    params = final_params,
    ll = data_ll,
    bic = bic,
    aic = aic,
    n_params = n_params,
    n_feat = n_feat,
    n_samples = n_samples,
    conv = final_conv,
    n_iters = result$counts[1],
    conv_code = result$convergence,
    runtime = runtime,
    model_type = "MAGIC",
    iters = result$counts[1],
    function_evals = result$counts[2]
  ))
}

check_input <- function(bsseq_data) {
  if (!inherits(bsseq_data, "BSseq")) stop("Input must be a BSseq object")
}

validate_methylation_data <- function(meth_data, cov_data, context = "data") {
  n_features <- nrow(meth_data)
  n_samples <- ncol(meth_data)
  n_obs_total <- n_features * n_samples
  
  cat(sprintf("\nValidating %s: %d features x %d samples (%s observations)\n", 
              context, n_features, n_samples, format(n_obs_total, big.mark = ",")))
  
  invalid_count <- 0
  
  if (any(!is.finite(meth_data))) {
    n_invalid_meth <- sum(!is.finite(meth_data))
    warning(sprintf("[VALIDATION] Found %s NaN/Inf values in methylation counts (%.2f%%)",
                   format(n_invalid_meth, big.mark = ","),
                   100 * n_invalid_meth / n_obs_total))
    invalid_count <- invalid_count + n_invalid_meth
  }
  
  if (any(!is.finite(cov_data))) {
    n_invalid_cov <- sum(!is.finite(cov_data))
    warning(sprintf("[VALIDATION] Found %s NaN/Inf values in coverage counts (%.2f%%)",
                   format(n_invalid_cov, big.mark = ","),
                   100 * n_invalid_cov / n_obs_total))
    invalid_count <- invalid_count + n_invalid_cov
  }
  
  if (any(meth_data < 0, na.rm = TRUE)) {
    n_negative_meth <- sum(meth_data < 0, na.rm = TRUE)
    warning(sprintf("[VALIDATION] Found %s negative methylation counts (%.2f%%)",
                   format(n_negative_meth, big.mark = ","),
                   100 * n_negative_meth / n_obs_total))
    invalid_count <- invalid_count + n_negative_meth
  }
  
  if (any(cov_data <= 0, na.rm = TRUE)) {
    n_nonpositive_cov <- sum(cov_data <= 0, na.rm = TRUE)
    warning(sprintf("[VALIDATION] Found %s non-positive coverage counts (%.2f%%)",
                   format(n_nonpositive_cov, big.mark = ","),
                   100 * n_nonpositive_cov / n_obs_total))
    invalid_count <- invalid_count + n_nonpositive_cov
  }
  
  if (any(meth_data > cov_data, na.rm = TRUE)) {
    violations <- meth_data > cov_data
    violations[is.na(violations)] <- FALSE
    n_violations <- sum(violations)
    warning(sprintf("[VALIDATION] Found %s observations where methylation > coverage (%.2f%%)",
                   format(n_violations, big.mark = ","),
                   100 * n_violations / n_obs_total))
    cat("  Suggestion: Verify that methylation and coverage matrices are correctly aligned\n")
    invalid_count <- invalid_count + n_violations
  }
  
  valid_obs <- is.finite(meth_data) & is.finite(cov_data) & 
               (meth_data >= 0) & (cov_data > 0) & (meth_data <= cov_data)
  
  n_valid <- sum(valid_obs)
  pct_valid <- 100 * n_valid / n_obs_total
  
  if (pct_valid < 50) {
    stop(sprintf("[ERROR] Only %.1f%% of observations are valid (need >= 50%%)\n  Suggestion: Check input data quality and format",
                pct_valid))
  }
  
  if (invalid_count > 0) {
    cat(sprintf("  Setting %s invalid observations to NA\n", 
                format(invalid_count, big.mark = ",")))
    meth_data[!valid_obs] <- NA
    cov_data[!valid_obs] <- NA
  }
  
  valid_per_feature <- rowSums(valid_obs)
  min_valid_feature <- min(valid_per_feature)
  if (min_valid_feature < 3) {
    n_low_features <- sum(valid_per_feature < 3)
    cat(sprintf("  Note: %d features have < 3 valid observations (will be filtered)\n",
                n_low_features))
  }
  
  valid_per_sample <- colSums(valid_obs)
  min_valid_sample <- min(valid_per_sample)
  if (min_valid_sample < 1000) {
    n_low_samples <- sum(valid_per_sample < 1000)
    cat(sprintf("  Note: %d samples have < 1000 valid observations (may affect mixture fitting)\n",
                n_low_samples))
  }
  
  mean_cov <- mean(cov_data[valid_obs])
  mean_meth_prop <- mean(meth_data[valid_obs] / cov_data[valid_obs])
  var_meth_prop <- var(meth_data[valid_obs] / cov_data[valid_obs])
  
  cat(sprintf("  Valid observations: %s (%.1f%%)\n", 
              format(n_valid, big.mark = ","), pct_valid))
  cat(sprintf("  Mean coverage: %.2fx\n", mean_cov))
  cat(sprintf("  Mean methylation: %.4f\n", mean_meth_prop))
  cat(sprintf("  Methylation variance: %.6f\n", var_meth_prop))
  
  if (mean_cov < 3) {
    warning(sprintf("[VALIDATION] Very low mean coverage (%.2fx) - results may have high uncertainty\n  Suggestion: Consider using data with higher sequencing depth",
                   mean_cov))
  } else if (mean_cov < 5) {
    warning(sprintf("[VALIDATION] Low mean coverage (%.2fx) - results may have high uncertainty",
                   mean_cov))
  }
  
  if (mean_meth_prop < 0.01) {
    warning(sprintf("[VALIDATION] Extremely low mean methylation (%.4f) - data may be problematic\n  Suggestion: Verify that input represents methylation, not unmethylated counts",
                   mean_meth_prop))
  } else if (mean_meth_prop > 0.99) {
    warning(sprintf("[VALIDATION] Extremely high mean methylation (%.4f) - data may be problematic\n  Suggestion: Verify that input represents methylation, not total counts",
                   mean_meth_prop))
  }
  
  if (var_meth_prop < 0.001) {
    warning(sprintf("[VALIDATION] Very low methylation variance (%.6f) - mixture models may be unnecessary\n  Suggestion: Data shows limited heterogeneity; consider simpler analysis methods",
                   var_meth_prop))
  }
  
  cat("\n")
  
  return(list(
    meth_data = meth_data,
    cov_data = cov_data,
    valid_mask = valid_obs,
    n_valid = n_valid,
    pct_valid = pct_valid
  ))
}

filter_data <- function(bsseq_data, args) {
  check_input(bsseq_data)
  cat(sprintf("Initial: %s sites, %d samples\n", format(length(bsseq_data), big.mark = ","), ncol(bsseq_data)))

  if (args$synthetic) {
    cat("Synthetic mode: skipping all filtering\n")
    cov_data <- getCoverage(bsseq_data, type = "Cov")
    meth_data <- getCoverage(bsseq_data, type = "M")
    
    if (nrow(meth_data) == 0 || ncol(meth_data) == 0) stop("No data available")
    
    cat(sprintf("Unfiltered: %s sites\n", format(nrow(meth_data), big.mark = ",")))
    return(list(meth_data = meth_data, cov_data = cov_data))
    }

  if (args$no_sex_chr) {
    chr_names <- as.character(seqnames(bsseq_data))
    sex_chr_patterns <- c("chrX", "chrY", "X", "Y")
    sex_chr_mask <- chr_names %in% sex_chr_patterns
    n_sex_chr_sites <- sum(sex_chr_mask)
    
    if (n_sex_chr_sites > 0) {
      cat(sprintf("Removing %s sites on sex chromosomes\n", format(n_sex_chr_sites, big.mark = ",")))
      bsseq_data <- bsseq_data[!sex_chr_mask, ]
    }
  }
  
  if (!is.null(args$subset_trait)) {
    sample_data <- pData(bsseq_data)
    if (!args$subset_trait %in% names(sample_data)) {
      stop(sprintf("Subset trait '%s' not found in sample metadata", args$subset_trait))
    }
    
    sample_traits <- sample_data[[args$subset_trait]]
    if (!args$subset_trait_val %in% sample_traits) {
      available_vals <- unique(sample_traits[!is.na(sample_traits)])
      stop(sprintf("Trait value '%s' not found. Available: %s", args$subset_trait_val, paste(available_vals, collapse = ", ")))
    }
    
    trait_mask <- sample_traits == args$subset_trait_val & !is.na(sample_traits)
    if (sum(trait_mask) == 0) stop(sprintf("No samples found for trait value '%s'", args$subset_trait_val))
    
    cat(sprintf("Filtering to trait '%s' = '%s': %d samples\n", args$subset_trait, args$subset_trait_val, sum(trait_mask)))
    bsseq_data <- bsseq_data[, trait_mask]
  }
  
  if (!is.null(args$chr)) {
    target_chr <- if (grepl("^chr", args$chr)) args$chr else paste0("chr", args$chr)
    bsseq_data <- bsseq_data[as.character(seqnames(bsseq_data)) == target_chr, ]
  }
  
  cov_data <- getCoverage(bsseq_data, type = "Cov")
  meth_data <- getCoverage(bsseq_data, type = "M")
  
  if (nrow(meth_data) == 0 || ncol(meth_data) == 0) stop("No data remaining after filtering")
  
  validated <- validate_methylation_data(meth_data, cov_data, context = "input data")
  meth_data <- validated$meth_data
  cov_data <- validated$cov_data
  
  if (args$filter_sample_outliers) {
    subset_sites <- sample(nrow(cov_data), min(10000, nrow(cov_data)))
    mean_cov_per_sample <- colMeans(cov_data[subset_sites, , drop = FALSE] + 1, na.rm = TRUE)
    
    low_cov_samples <- mean_cov_per_sample < args$sample_thresh
    if (sum(low_cov_samples) > 0) {
      cat(sprintf("Removing %d low-coverage samples\n", sum(low_cov_samples)))
      cov_data <- cov_data[, !low_cov_samples, drop = FALSE]
      meth_data <- meth_data[, !low_cov_samples, drop = FALSE]
    }
  }
  
  valid_mask <- cov_data >= args$min_reads
  valid_mask[is.na(valid_mask)] <- FALSE
  samples_per_feat <- rowSums(valid_mask)
  valid_feats <- samples_per_feat >= args$min_samples
  
  if (sum(valid_feats, na.rm = TRUE) == 0) stop("No features pass filtering criteria")
  
  meth_data <- meth_data[valid_feats, , drop = FALSE]
  cov_data <- cov_data[valid_feats, , drop = FALSE]
  valid_mask <- valid_mask[valid_feats, , drop = FALSE]
  
  props_matrix <- meth_data / cov_data
  props_matrix[!valid_mask] <- NA
  vars <- rowVars(props_matrix, na.rm = TRUE)
  vars[is.na(vars)] <- 0
  
  valid_vars <- vars[vars > 0 & is.finite(vars)]
  var_thresh <- max(args$min_var, quantile(valid_vars, 0.25, na.rm = TRUE))
  var_pass <- vars >= var_thresh & is.finite(vars)
  var_pass[is.na(var_pass)] <- FALSE
  
  cat(sprintf("Filtered: %s sites\n", format(sum(var_pass, na.rm = TRUE), big.mark = ",")))
  
  filtered_result <- list(meth_data = meth_data[var_pass, , drop = FALSE], 
                         cov_data = cov_data[var_pass, , drop = FALSE])
  
  n_final_features <- nrow(filtered_result$meth_data)
  n_final_samples <- ncol(filtered_result$meth_data)
  
  if (n_final_features == 0) {
    stop("[ERROR] No sites remain after filtering")
  }
  
  if (n_final_features < 100) {
    stop(sprintf("[ERROR] Only %d features remain after filtering (need >= 100 for mixture modeling)\n  Suggestion: Reduce filtering stringency (--min-var, --min-reads, --min-samples)",
                n_final_features))
  }
  
  if (n_final_samples < 2) {
    stop(sprintf("[ERROR] Only %d samples remain (need >= 2)\n  Suggestion: Check sample filtering criteria (--sample-thresh)",
                n_final_samples))
  }
  
  if (n_final_features < 500) {
    warning(sprintf("[VALIDATION] Low feature count (%d) - model may be less stable\n  Suggestion: Consider reducing filtering stringency for better estimates",
                   n_final_features))
  }
  
  return(filtered_result)
}

make_holdout_split <- function(meth_data, total_data, use_holdout = FALSE, frac = 0.2) {
  if (!use_holdout) {
    return(list(
      train_meth = meth_data, train_total = total_data,
      test_meth = meth_data, test_total = total_data,
      train_indices = 1:nrow(meth_data), test_indices = 1:nrow(meth_data),
      strategy = "none"
    ))
  }
  
  set.seed(DEFAULT_HOLDOUT_SEED)
  n_feat <- nrow(meth_data)
  n_test <- floor(n_feat * frac)
  test_indices <- sample(n_feat, n_test)
  train_indices <- setdiff(1:n_feat, test_indices)
  
  if (length(train_indices) == 0 || length(test_indices) == 0) {
    stop("Holdout split resulted in empty train or test set")
  }
  
  cat(sprintf("Holdout split: %d train, %d test sites\n", length(train_indices), length(test_indices)))
  
  return(list(
    train_meth = meth_data[train_indices, , drop = FALSE],
    train_total = total_data[train_indices, , drop = FALSE],
    test_meth = meth_data[test_indices, , drop = FALSE],
    test_total = total_data[test_indices, , drop = FALSE],
    train_indices = train_indices, test_indices = test_indices, strategy = "random"
  ))
}

test_ll <- function(test_meth, test_total, params, k_comp, min_reads = 5) {
  valid_mask <- test_total >= min_reads
  log_params <- natural_to_log(params$alpha, params$beta, params$pi)
  feature_ll <- compute_ll(test_meth, test_total, valid_mask, log_params)
  return(sum(feature_ll[is.finite(feature_ll)]))
}

run_attempt <- function(split_data, k_comp, max_iters, tol, min_reads, 
                       run_idx, conv_seed, verbose, seed) {
  if (verbose) cat(sprintf("  Run %d... ", run_idx))
  
  set.seed(conv_seed + run_idx * 1000)
  
  result <- optimize_model(
    meth_mat = split_data$train_meth, total_mat = split_data$train_total,
    k_comp = k_comp, max_iters = max_iters, tol = tol, min_reads = min_reads, seed = seed
  )
  
  test_ll_val <- NA
  overfitting <- NA
  
  if (split_data$strategy != "none") {
    test_ll_val <- test_ll(
      test_meth = split_data$test_meth, test_total = split_data$test_total,
      params = result$params, k_comp = k_comp, min_reads = min_reads
    )
    
    train_per_sample <- result$ll / nrow(split_data$train_meth)
    test_per_sample <- test_ll_val / nrow(split_data$test_meth)
    overfitting <- train_per_sample - test_per_sample
  }
  
  result$run_id <- run_idx
  result$train_ll <- result$ll
  result$test_ll <- test_ll_val
  result$overfitting <- overfitting
  result$holdout_strat <- split_data$strategy
  result$n_train_feat <- nrow(split_data$train_meth)
  result$n_test_feat <- nrow(split_data$test_meth)
  
  if (verbose) {
    cat(sprintf("%.1fs, train_ll=%.1f", result$runtime, result$train_ll))
    if (!is.na(test_ll_val)) cat(sprintf(", test_ll=%.1f", test_ll_val))
    if (!is.na(overfitting)) cat(sprintf(", overfit=%.3f", overfitting))
    cat(sprintf(", %s\n", if(result$conv) "OK" else "FAIL"))
  }
  
  return(result)
}

select_best_result <- function(all_results) {
  conv_results <- all_results[sapply(all_results, function(x) x$conv)]
  if (length(conv_results) == 0) stop("All conv attempts failed")
  train_ll <- sapply(conv_results, function(x) x$train_ll)
  return(conv_results[[which.max(train_ll)]])
}

prep_result <- function(best_result, all_results, holdout_strat) {
  result <- best_result
  result$n_runs <- length(all_results)
  result$n_conv <- sum(sapply(all_results, function(x) x$conv))
  result$all_results <- all_results
  result$holdout_strat <- holdout_strat
  return(result)
}

run_optimization <- function(split_data, k_comp, max_iters, tol, min_reads = 5, 
                            n_runs = 1, conv_seed = DEFAULT_CONV_SEED, seed = NULL) {
  if (n_runs > 1) {
    cat(sprintf("Running %d conv attempts for MAGIC K=%d\n", n_runs, k_comp))
  }
  
  all_results <- list()
  for (run_idx in 1:n_runs) {
    all_results[[run_idx]] <- run_attempt(
      split_data, k_comp, max_iters, tol, min_reads, 
      run_idx, conv_seed, n_runs > 1, seed
    )
  }
  
  best_result <- select_best_result(all_results)
  result <- prep_result(best_result, all_results, split_data$strategy)
  
  if (n_runs > 1) {
    cat(sprintf("Best run: %d/%d conv, train_ll=%.1f", 
                result$n_conv, n_runs, best_result$train_ll))
    if (!is.na(best_result$test_ll)) {
      cat(sprintf(", test_ll=%.1f, overfit=%.3f", best_result$test_ll, best_result$overfitting))
    }
    cat("\n")
  }
  
  return(result)
}

make_output_dir <- function(input_file, args) {
  base_name <- basename(tools::file_path_sans_ext(input_file))
  base_name <- sub("^BSseq_", "", base_name)
  base_name <- sub("_mergedStrands", "", base_name)
  base_name <- sub("_minCov\\d+", "", base_name)
  base_name <- sub("_maxCov\\d+", "", base_name)
  base_name <- sub("_minSamples\\d+", "", base_name)
  base_name <- sub("pediatric", "ped", base_name)
  
  dir_parts <- c(base_name, sprintf("K%d", args$k_components))
  if (args$synthetic) dir_parts <- c(dir_parts, "synth")
  if (args$no_sex_chr) dir_parts <- c(dir_parts, "noXY")
  if (!is.null(args$subset_trait)) dir_parts <- c(dir_parts, sprintf("trait_%s", args$subset_trait_val))
  if (args$filter_sample_outliers) dir_parts <- c(dir_parts, "clean")
  
  dir_parts <- c(dir_parts, sprintf("S%d", args$min_samples))
  
  if (args$n_conv_runs > 1) dir_parts <- c(dir_parts, sprintf("n%d", args$n_conv_runs))
  if (args$holdout) dir_parts <- c(dir_parts, "ho")
  
  if (!is.null(args$seed_methylation_file)) {
    seed_suffix <- "methSeed"
    dir_parts <- c(dir_parts, seed_suffix)
  }
  
  tol_str <- sprintf("e%02d", -log10(args$tol))
  iters_str <- sprintf("it%d", args$max_iters)
  
  chr_part <- NULL
  if (!is.null(args$chr)) {
    chr_part <- if (grepl("^chr", args$chr)) args$chr else paste0("chr", args$chr)
  }
  
  all_parts <- c(dir_parts, tol_str, iters_str, chr_part)
  all_parts <- all_parts[!is.null(all_parts) & all_parts != ""]
  
  dir_name <- paste(all_parts, collapse = "_")
  
  if (dir.exists(dir_name)) {
    timestamp <- format(Sys.time(), "%m%d_%H%M")
    dir_name <- paste0(dir_name, "_", timestamp)
  }
  
  dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)
  return(dir_name)
}

write_results <- function(output_dir, results, args) {
  if (!is.null(results$params)) {
    mu_phi_list <- lapply(1:length(results$params$alpha), function(k) {
      ab_to_mu_phi(results$params$alpha[k], results$params$beta[k])
    })
    
    components_df <- data.frame(
      component = 1:length(results$params$pi),
      mixWeight = round(results$params$pi, 4),
      meanMeth = round(sapply(mu_phi_list, function(x) x[1]), 4),
      methDisp = round(sapply(mu_phi_list, function(x) x[2]), 5)
    )
    
    write.csv(components_df, file.path(output_dir, "methComps.csv"), row.names = FALSE)
    
    model_df <- data.frame(
      component = 1:length(results$params$pi),
      alpha = round(results$params$alpha, 4),
      beta = round(results$params$beta, 4),
      pi = round(results$params$pi, 4)
    )
    
    write.csv(model_df, file.path(output_dir, "optModel.csv"), row.names = FALSE)
  }
  
  n_feat <- results$n_feat
  n_samples <- results$n_samples
  n_train_feat <- results$n_train_feat
  n_test_feat <- results$n_test_feat
  
  summary_df <- data.frame(
    k = args$k_components, model_type = results$model_type,
    bic = round(results$bic, 1), aic = round(results$aic, 1),
    train_ll = round(results$train_ll, 1),
    test_ll = if(!is.na(results$test_ll)) round(results$test_ll, 1) else NA,
    overfitting = if(!is.na(results$overfitting)) round(results$overfitting, 4) else NA,
    nParams = results$n_params, nFeat = n_feat, nSamples = n_samples,
    nTrainFeat = n_train_feat, nTestFeat = n_test_feat,
    conv = results$conv, n_runs = results$n_runs, n_conv = results$n_conv,
    holdout_strat = results$holdout_strat, runtime = round(results$runtime, 1),
    nIters = results$n_iters, convCode = results$conv_code,
    sample_outliers_filtered = args$filter_sample_outliers,
    sample_thresh = if(args$filter_sample_outliers) args$sample_thresh else NA,
    no_sex_chr = args$no_sex_chr,
    subset_filter = if(!is.null(args$subset_trait)) args$subset_trait_val else "none",
    subset_trait = if(!is.null(args$subset_trait)) args$subset_trait else NA,
    seed_file = if(!is.null(args$seed_methylation_file)) basename(args$seed_methylation_file) else NA
  )
  write.csv(summary_df, file.path(output_dir, "runSum.csv"), row.names = FALSE)
  
  all_conv_data <- list()
  for (run_result in results$all_results) {
    conv_data <- data.frame(
      k = run_result$k, run_id = run_result$run_id, model_type = run_result$model_type,
      conv = run_result$conv, convCode = run_result$conv_code,
      train_ll = round(run_result$train_ll, 1),
      test_ll = if(!is.na(run_result$test_ll)) round(run_result$test_ll, 1) else NA,
      overfitting = if(!is.na(run_result$overfitting)) round(run_result$overfitting, 4) else NA,
      bic = round(run_result$bic, 1), aic = round(run_result$aic, 1),
      nParams = run_result$n_params, nIters = run_result$n_iters,
      runtime = round(run_result$runtime, 1), holdout_strat = run_result$holdout_strat
    )
    all_conv_data[[paste0("run", run_result$run_id)]] <- conv_data
  }
  
  all_conv_df <- do.call(rbind, all_conv_data)
  write.csv(all_conv_df, file.path(output_dir, "allConvRuns.csv"), row.names = FALSE)
  
  param_est <- list()
  for (i in 1:length(results$all_results)) {
    run_result <- results$all_results[[i]]
    if (run_result$conv) {
      params <- run_result$params
      param_data <- data.frame(
        run_id = i, component = 1:length(params$pi),
        pi = round(params$pi, 3), 
        alpha = round(params$alpha, 3), 
        beta = round(params$beta, 3),
        mean_meth = round(params$alpha / (params$alpha + params$beta), 3)
      )
      
      param_est[[i]] <- param_data
    }
  }
  
  if (length(param_est) > 0) {
    all_param_df <- do.call(rbind, param_est)
    write.csv(all_param_df, file.path(output_dir, "allParamEst.csv"), row.names = FALSE)
  }
  
  result_df <- data.frame(
    k = args$k_components,
    conv_code = results$conv_code,
    iters = results$iters,
    function_evals = results$function_evals
  )
  write.csv(result_df, file.path(output_dir, "optResult.csv"), row.names = FALSE)
}

if (sys.nframe() == 0) {
  start_time <- Sys.time()
  
  cat("MAGIC Methylation Analysis with Genomic Inferred Contexts\n")
  
  seed <- read_seed(args$seed_methylation_file, args$k_components)
  
  if (!is.null(seed)) {
    cat("Mode: Using seed params (INITIAL), optimizing all params\n")
  } else {
    cat("Mode: Default init for all params\n")
  }
  
  cat(sprintf("Loading %s\n", basename(args$input)))
  bsseq_data <- qs::qread(args$input)
  
  if (!inherits(bsseq_data, "BSseq")) stop("Input file must contain a BSseq object")
  
  filters <- c()
  if (args$no_sex_chr) filters <- c(filters, "excluding sex chromosomes")
  if (!is.null(args$subset_trait)) filters <- c(filters, sprintf("trait=%s", args$subset_trait_val))
  if (!is.null(args$chr)) filters <- c(filters, sprintf("chr=%s", args$chr))
  
  if (length(filters) > 0) {
    cat(sprintf("Filtering options: %s\n", paste(filters, collapse = ", ")))
  }
  
  filtered <- filter_data(bsseq_data, args)
  rm(bsseq_data); gc(verbose = FALSE)
  
  cat(sprintf("Data: %s sites × %d samples\n", 
              format(nrow(filtered$meth_data), big.mark = ","), ncol(filtered$meth_data)))
  
  split_data <- make_holdout_split(
    meth_data = filtered$meth_data, total_data = filtered$cov_data,
    use_holdout = args$holdout, frac = args$holdout_frac
  )
  
  output_dir <- make_output_dir(args$input, args)
  cat(sprintf("Output: %s\n", output_dir))
  
  info <- sprintf("%d runs", args$n_conv_runs)
  if (args$holdout) info <- sprintf("%s, holdout", info)
  if (!is.null(seed)) info <- sprintf("%s, seeded meth", info)
  
  cat(sprintf("Optimization: K=%d (%s)\n", args$k_components, info))
  
  result <- run_optimization(
    split_data = split_data, k_comp = args$k_components, max_iters = args$max_iters,
    tol = args$tol, min_reads = args$min_reads, n_runs = args$n_conv_runs,
    conv_seed = args$conv_seed, seed = seed
  )
  
  write_results(output_dir, result, args)
  
  cat(sprintf("Final: BIC=%.1fK, %d/%d conv, %.1fs\n", 
              result$bic/1000, result$n_conv, result$n_runs, result$runtime))
  
  if (!is.null(result$params)) {
    cat("\nFinal model summary:\n")
    for (k in 1:length(result$params$pi)) {
      mu_phi <- ab_to_mu_phi(result$params$alpha[k], result$params$beta[k])
      cat(sprintf("  Component %d: π=%.3f, μ=%.3f, φ=%.4f\n", 
                  k, result$params$pi[k], mu_phi[1], mu_phi[2]))
    }
  }
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("Complete: %.1fs\n", total_time))
}
