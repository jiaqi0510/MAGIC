# MAGIC

**Methylation Analysis with Genomic Inferred Contexts**

Beta-binomial mixture models for DNA methylation analysis, featuring robust parameter estimation and Bayesian differential testing.

## Overview

MAGIC implements a statistical framework for analyzing DNA methylation data using beta-binomial mixture models. The method:

- Estimates mixture model parameters via maximum likelihood with numerical optimization
- Performs Bayesian differential methylation testing between conditions
- Provides component-wise dispersion estimation for improved statistical power
- Utilizes OpenMP parallelization for efficient large-scale analysis

## Features

- **Mixture Model Optimization**: Maximum likelihood estimation of beta-binomial mixture components
- **Differential Testing**: Bayesian framework for identifying differentially methylated regions
- **Model Comparison**: Utilities for evaluating parameter recovery and model performance
- **Scalable Computing**: OpenMP parallelization for multi-core systems
- **Convergence Diagnostics**: Multiple initialization runs with convergence tracking
- **Flexible Testing**: Support for case-control comparisons and correlation analyses

## Installation

### System Requirements

- **R** ≥ 4.0.0
- **C++ compiler** with C++11 support and OpenMP
  - Linux: GCC ≥ 4.9 or Clang with OpenMP
  - macOS: GCC via Homebrew recommended (see notes below)
  - Windows: Rtools compatible with your R version

### Required R Packages

```r
# CRAN packages
install.packages(c("Rcpp", "optparse", "data.table", "ggplot2", "gridExtra"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("bsseq", "GenomicRanges", "matrixStats"))
```

### Installation from GitHub

```bash
git clone https://github.com/Pickledzebra/MAGIC.git
cd MAGIC
```

The C++ code will be compiled automatically by Rcpp when first run. Alternatively, you can pre-compile:

```bash
R CMD SHLIB src/magicFit.cpp -fopenmp
R CMD SHLIB src/magicTest.cpp -fopenmp
```

### macOS Installation Notes

Apple Clang does not include OpenMP by default. Install via Homebrew:

```bash
brew install libomp

# Set environment variables before running R:
export CPPFLAGS="-Xpreprocessor -fopenmp"
export LDFLAGS="-lomp"
```

Alternatively, install GCC via Homebrew and specify it as your compiler:

```bash
brew install gcc
# Then use gcc-13 (or your version) when compiling
```

## Quick Start

**Important:** All commands should be run from the MAGIC repository root directory.

```bash
cd MAGIC  # Make sure you're in the repository root
```

### Quick Test (No Data Required)

Test your installation with synthetic data:

```bash
Rscript examples/example_synthetic_workflow.R
```

This generates test data, fits a model, and validates parameter recovery in ~3-5 minutes.

### Parameter Estimation

Estimate a K=3 component mixture model:

```bash
Rscript R/magicFit.R \
  --input data/BSseq_data.qs \
  --k_components 3 \
  --max_iters 500 \
  --tol 1e-6 \
  --threads 4
```

### Differential Methylation Testing

Test for differential methylation between conditions:

```bash
Rscript R/magicTest.R \
  --input data/BSseq_data.qs \
  --group condition \
  --compare case,control \
  --mixture_model results/data_K3_S10_e06_it500/ \
  --cores 4
```

### Generate Synthetic Data (for testing)

Generate synthetic methylation data with known parameters:

```bash
Rscript scripts/generate_synthetic_data.R \
  -K 3 \
  -n 100000 \
  -s 50 \
  -o synthetic_K3.qs
```

This creates test data with known ground truth for validating the method.

### Model Comparison (for simulations)

Compare fitted parameters to true values:

```bash
Rscript scripts/model_comparison.R \
  -t simulations/true_models/ \
  -f results/fitted_models/ \
  -o comparison_summary.csv \
  -p performance_plots.pdf
```

## Usage

### Optimization Module

```bash
Rscript R/magicFit.R [options]
```

**Required Arguments:**
- `-i, --input <file>` - BSseq data file (`.qs` format)
- `-k, --k_components <int>` - Number of mixture components

**Optional Arguments:**
- `--max_iters <int>` - Maximum optimization iterations (default: 500)
- `-t, --tol <float>` - Convergence tolerance (default: 1e-6)
- `-r, --min_reads <int>` - Minimum read coverage per site (default: 5)
- `-m, --min_samples <int>` - Minimum samples per component (default: 10)
- `--threads <int>` - Number of OpenMP threads (default: auto-detect)
- `--n_conv_runs <int>` - Number of convergence runs (default: 1)
- `--holdout` - Use holdout validation
- `--holdout_frac <float>` - Fraction for holdout set (default: 0.2)
- `--chr <string>` - Restrict to specific chromosome
- `--seed_methylation_file <file>` - Initialize from seed parameters
- `-q, --quiet` - Suppress progress messages

### Testing Module

```bash
Rscript R/magicTest.R [options]
```

**Required Arguments:**
- `-i, --input <file>` - BSseq data file (`.qs` format)
- `-g, --group <string>` - Group column name in sample metadata
- `--mixture_model <dir>` - Path to fitted mixture model directory

**Optional Arguments:**
- `-c, --compare <string>` - Group types to compare (comma-separated, e.g., 'case,control')
- `-s, --subset <string>` - Subset group column for filtering
- `-k, --keep <string>` - Subset group type to retain
- `-r, --min_reads <int>` - Minimum read coverage (default: 5)
- `-m, --min_samp <int>` - Minimum samples per group (default: 3)
- `-n, --cores <int>` - Number of CPU cores (default: 1)
- `--chr <string>` - Restrict to specific chromosome
- `--methods <string>` - Statistical methods: mixture, magic, all (default: all)
- `--prior_strength <float>` - Dirichlet prior concentration (default: 1.0)
- `-z, --randomize` - Randomize group assignments (null testing)
- `-q, --quiet` - Suppress progress messages

### Model Comparison Utility

```bash
Rscript scripts/model_comparison.R [options]
```

**Required Arguments:**
- `-t, --true-dir <dir>` - Directory with true model files (*_trueModel.csv)
- `-f, --fitted-dir <dir>` - Directory with fitted model subdirectories (*_K*_*)

**Optional Arguments:**
- `-o, --output <file>` - Save comparison summary to CSV
- `-p, --plots <file>` - Generate performance plots (PDF)
- `-h, --help` - Show help message

## Output Files

### Optimization Outputs

Output directory: `[prefix]_K[n]_[options]/`

- **optModel.csv** - Fitted parameters in natural space
  - Columns: `component`, `alpha`, `beta`, `pi`
  
- **methComps.csv** - Interpretable parameter representation
  - Columns: `component`, `mixWeight`, `meanMeth`, `methDisp`
  
- **runSum.csv** - Optimization summary
  - BIC, AIC, log-likelihood, convergence status, runtime
  
- **allParamEst.csv** - Parameters from all convergence runs
  - Useful for assessing convergence stability
  
- **allConvRuns.csv** - Convergence diagnostics for each run
  - Convergence codes, iteration counts, objective values

### Testing Outputs

Output directory: `[analysis]_magic_results_[timestamp]/`

- **results.csv** - Main results table
  - Genomic coordinates, mean methylation by group
  - P-values, Bayes factors, effect sizes
  - Dominant mixture components per group
  
- **magic_summary.txt** - Analysis summary
  - Dataset characteristics, filtering steps
  - Method descriptions, significance thresholds

### Comparison Outputs

- **comparison_summary.csv** - RMSE metrics for each parameter
- **performance_plots.pdf** - Visualization of parameter recovery

## Examples

See the `examples/` directory for complete working examples:

- **example_synthetic_workflow.R** - Complete workflow with synthetic data (start here!)
- **example_optimization.R** - Basic optimization workflow
- **example_testing.R** - Differential methylation testing
- **example_comparison.R** - Model comparison for simulations

## Method

MAGIC models methylation read counts using a beta-binomial mixture:

For site *i* with *m* methylated reads out of *n* total reads:

```
m_i ~ Σ_k π_k * BetaBinomial(n_i, α_k, β_k)
```

Where:
- *π_k* are mixture weights (sum to 1)
- *α_k*, *β_k* are beta-binomial shape parameters
- Mean methylation: *μ_k = α_k / (α_k + β_k)*
- Dispersion: *φ_k = 1 / (α_k + β_k + 1)*

**Parameter Estimation:**
- Maximum likelihood via L-BFGS-B optimization
- Log-space parameterization for numerical stability
- Multiple random initializations for global optimization

**Differential Testing:**
- Bayesian framework with mixture-based priors
- Component-wise dispersion estimation
- Wald tests for mean differences

See `docs/method.md` for mathematical details.

## Performance

MAGIC is designed for genome-scale analysis:

- **Parallelization**: OpenMP multi-threading for both optimization and testing
- **Memory Efficiency**: Chunked processing of large datasets
- **Numerical Stability**: Log-space calculations, specialized functions
- **Convergence**: Multiple runs with best-fit selection

**Typical Runtime** (4 cores, ~100K CpG sites):
- Optimization (K=3): 2-5 minutes
- Differential testing: 5-10 minutes

## Citation

If you use MAGIC in your research, please cite:

```bibtex
@article{MAGIC2025,
  author = {[Authors]},
  title = {MAGIC: Methylation Analysis with Genomic Inferred Contexts},
  journal = {[Journal]},
  year = {2025},
  note = {Manuscript in preparation},
  url = {https://github.com/Pickledzebra/MAGIC}
}
```

See `CITATION.bib` for the complete citation.

## Documentation

- **Installation Guide**: `docs/installation.md` - Detailed compilation instructions
- **Usage Guide**: `docs/usage.md` - Comprehensive parameter descriptions
- **Method Description**: `docs/method.md` - Mathematical formulation

## License

[SPECIFY LICENSE HERE - e.g., MIT, GPL-3, BSD-3-Clause]

See `LICENSE` file for details.

## Contributing

This software accompanies a manuscript submission. For questions or issues:

- Open an issue on GitHub: https://github.com/Pickledzebra/MAGIC/issues
- Contact: [maintainer email]

## Authors

[AUTHOR LIST WITH AFFILIATIONS]

## Acknowledgments

[FUNDING SOURCES, COLLABORATORS, ETC.]
