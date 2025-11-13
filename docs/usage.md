# Usage Guide

Comprehensive guide to using MAGIC for methylation analysis.

## Overview

MAGIC consists of three main modules:

1. **Optimization** (`R/magicFit.R`) - Estimate mixture model parameters
2. **Testing** (`R/magicTest.R`) - Differential methylation analysis
3. **Comparison** (`scripts/model_comparison.R`) - Evaluate parameter recovery

## Typical Workflow

```
1. Prepare BSseq data → 2. Optimize model → 3. Test differences → 4. Interpret results
```

---

---

## Synthetic Data Generation

For testing and validation, generate synthetic methylation data with known parameters.

### Basic Usage

```bash
Rscript scripts/generate_synthetic_data.R \
  -K 3 \
  -n 100000 \
  -s 50 \
  -o synthetic_data.qs
```

### Complete Parameter Reference

#### Required Parameters

**`-K, --components <int>`**
- Number of mixture components
- Typical values: 2-5

#### Optional Parameters

**`-n, --cpgs <int>`** (default: 1000000)
- Number of CpG sites to generate
- Reduce for faster testing

**`-s, --samples <int>`** (default: 100)
- Number of samples

**`-o, --output <file>`** (default: synthetic_bb_K{K}.qs)
- Output filename
- Also creates {filename}_trueModel.csv with true parameters

**`--seed <int>`** (default: 12345)
- Random seed for reproducibility

**`--no-filter`**
- Skip magicFit filtering
- By default, applies same filtering as optimization

**`--min-reads <int>`** (default: 5)
- Minimum read coverage for filtering

**`--min-samples <int>`** (default: 6)
- Minimum samples per site for filtering

**`--min-var <float>`** (default: 0.001)
- Minimum variance threshold for filtering

**`--sparse`**
- Generate sparse coverage (WGBS-like)
- Creates realistic missing data patterns

**`--sparsity <float>`** (default: 0.8)
- Fraction of zero coverage entries
- Only used with `--sparse`

**`--sex-diff`**
- Add sex-specific differential methylation
- Creates chrX with dosage compensation patterns

**`--chrx-fraction <float>`** (default: 0.05)
- Fraction of CpGs designated as chrX sites
- Only used with `--sex-diff`

### Output Files

1. **{output}.qs** - BSseq object with synthetic data
   - Contains methylation and coverage matrices
   - Includes sample metadata
   - True parameters stored in metadata

2. **{output}_trueModel.csv** - True parameters
   ```
   component,alpha,beta,pi
   1,1.0,15.0,0.333
   2,3.0,3.0,0.333
   3,15.0,1.0,0.333
   ```

### Component Parameters

The generator creates realistic beta-binomial parameters:

**K=2:**
- Component 1: Low methylation (α=1, β=10)
- Component 2: High methylation (α=10, β=1)

**K=3:**
- Component 1: Unmethylated (α=1, β=15)
- Component 2: Partially methylated (α=3, β=3)
- Component 3: Highly methylated (α=15, β=1)

**K>3:**
- Evenly spaced from low to high methylation

### Common Usage Patterns

#### Quick Test Data
```bash
Rscript scripts/generate_synthetic_data.R \
  -K 3 \
  -n 10000 \
  -s 20
```

#### Realistic WGBS
```bash
Rscript scripts/generate_synthetic_data.R \
  -K 3 \
  -n 500000 \
  -s 100 \
  --sparse \
  --sparsity 0.7
```

#### Sex-Specific Patterns
```bash
Rscript scripts/generate_synthetic_data.R \
  -K 3 \
  -n 200000 \
  -s 80 \
  --sex-diff \
  --chrx-fraction 0.05
```

#### Multiple K Values
```bash
for k in 2 3 4 5; do
  Rscript scripts/generate_synthetic_data.R \
    -K $k \
    -n 100000 \
    -o synthetic_K${k}.qs
done
```

### Use Cases

**Method Validation:**
- Compare fitted vs true parameters
- Assess parameter recovery accuracy
- Validate convergence

**Testing Installation:**
- Quick functionality check
- No real data required
- Known ground truth

**Benchmarking:**
- Test computational performance
- Evaluate scaling with data size
- Compare optimization settings

**Development:**
- Test new features
- Debug issues
- Rapid iteration

---

## Module 1: Parameter Optimization

### Basic Usage

```bash
Rscript R/magicFit.R \
  --input data.qs \
  --k_components 3 \
  --threads 4
```

### Complete Parameter Reference

#### Required Parameters

**`-i, --input <file>`**
- BSseq data file in `.qs` format
- Must contain methylation counts and coverage
- Created using `bsseq` package in R

**`-k, --k_components <int>`**
- Number of mixture components
- Typical values: 2-5
- Higher K increases model flexibility but risks overfitting

#### Optimization Parameters

**`--max_iters <int>`** (default: 500)
- Maximum optimization iterations
- Increase if convergence not reached
- Typical range: 200-1000

**`-t, --tol <float>`** (default: 1e-6)
- Convergence tolerance
- Smaller values = stricter convergence
- Typical range: 1e-4 to 1e-8

**`--n_conv_runs <int>`** (default: 1)
- Number of independent optimization runs
- Uses different random initializations
- Recommended: 3-10 for important analyses
- Best fit (by BIC) is returned

**`--conv_seed <int>`** (default: 54321)
- Seed for random initialization
- Use for reproducibility
- Different seeds explore different local optima

#### Data Filtering Parameters

**`-r, --min_reads <int>`** (default: 5)
- Minimum read coverage per CpG site
- Sites with fewer reads are excluded
- Higher values = more stringent filtering

**`-m, --min_samples <int>`** (default: 10)
- Minimum number of samples with valid data per site
- Higher values = more reliable estimates
- Balance with dataset size

**`-v, --min_var <float>`** (default: 0.001)
- Minimum methylation variance threshold
- Filters invariant sites
- Reduces computational burden

**`--chr <string>`**
- Restrict analysis to specific chromosome
- Format: "chr1" or "1"
- Useful for testing or memory constraints
- Example: `--chr chr21`

**`--no_sex_chr`**
- Exclude sex chromosomes (chrX, chrY)
- Recommended for mixed-sex cohorts
- No additional argument needed (flag only)

**`--filter_sample_outliers`**
- Remove sample outliers based on global methylation
- Uses `--sample_thresh` parameter
- Can improve model fit

**`--sample_thresh <float>`** (default: 3.0)
- Z-score threshold for outlier detection
- Only used with `--filter_sample_outliers`
- Typical values: 2.5-3.5

**`--subset_trait <string>`**
- Filter samples by metadata column
- Must exist in BSseq colData
- Example: `--subset_trait tissue`

**`--subset_trait_val <string>`**
- Value to keep for `--subset_trait`
- Example: `--subset_trait_val liver`

#### Validation Parameters

**`--holdout`**
- Use holdout validation
- Splits data into training and testing sets
- Provides overfitting assessment

**`--holdout_frac <float>`** (default: 0.2)
- Fraction of sites for holdout set
- Only used with `--holdout`
- Typical values: 0.1-0.3

#### Initialization Parameters

**`--seed_methylation_file <file>`**
- Initialize from pre-specified parameters
- CSV file with columns: `component`, `alpha`, `beta`, `pi`
- Or: `mixWeight`, `meanMeth`, `methDisp`
- Useful for warm starts or constrained optimization

#### Computational Parameters

**`--threads <int>`**
- Number of OpenMP threads
- Auto-detects if not specified
- Recommended: number of physical cores
- Example: `--threads 8`

#### Other Parameters

**`-q, --quiet`**
- Suppress progress messages
- Useful for batch processing

**`--synthetic`**
- Skip filtering for synthetic data
- Assumes data is already clean

### Output Structure

Creates directory: `[basename]_K[k]_[options]/`

Example: `data_K3_S10_e06_it500/`

**Output Files:**

1. **optModel.csv** - Final parameters
   ```
   component,alpha,beta,pi
   1,0.523,12.456,0.234
   2,5.123,8.234,0.456
   3,15.234,2.345,0.310
   ```

2. **methComps.csv** - Interpretable parameters
   ```
   component,mixWeight,meanMeth,methDisp
   1,0.234,0.040,0.072
   2,0.456,0.384,0.071
   3,0.310,0.867,0.054
   ```

3. **runSum.csv** - Summary statistics
   - BIC, AIC, log-likelihood
   - Convergence status, iterations
   - Runtime, parameters

4. **allParamEst.csv** - All convergence runs
   - Parameters from each run
   - Useful for stability assessment

5. **allConvRuns.csv** - Convergence diagnostics
   - Per-run convergence codes
   - Iteration counts, objective values

### Parameter Selection Guidelines

#### Choosing K (Number of Components)

Run multiple values and compare BIC:

```bash
for k in 2 3 4 5; do
  Rscript R/magicFit.R \
    --input data.qs \
    --k_components $k \
    --n_conv_runs 5
done

# Compare BIC values in runSum.csv files
```

**Guidelines:**
- Start with K=2-3 for initial exploration
- Lower BIC indicates better model
- Balance fit quality with interpretability
- Diminishing returns typically beyond K=5

#### Setting Convergence Parameters

**For exploratory analysis:**
```bash
--max_iters 200 --tol 1e-4 --n_conv_runs 1
```

**For publication:**
```bash
--max_iters 500 --tol 1e-6 --n_conv_runs 10
```

**For difficult datasets:**
```bash
--max_iters 1000 --tol 1e-8 --n_conv_runs 20
```

#### Memory and Performance

**Large datasets** (>500K sites):
```bash
--chr chr1  # Process one chromosome at a time
--min_samples 15  # More stringent filtering
--threads 4  # Fewer threads may reduce memory pressure
```

**Small datasets** (<50K sites):
```bash
--min_samples 5  # Less stringent filtering
--threads 8  # More parallelization
```

### Common Usage Patterns

#### Standard Analysis
```bash
Rscript R/magicFit.R \
  --input BSseq_data.qs \
  --k_components 3 \
  --n_conv_runs 5 \
  --threads 8
```

#### Chromosome-Specific
```bash
Rscript R/magicFit.R \
  --input BSseq_data.qs \
  --k_components 3 \
  --chr chr21 \
  --threads 4
```

#### With Validation
```bash
Rscript R/magicFit.R \
  --input BSseq_data.qs \
  --k_components 3 \
  --holdout \
  --holdout_frac 0.2 \
  --n_conv_runs 10
```

#### Subset Analysis
```bash
Rscript R/magicFit.R \
  --input BSseq_data.qs \
  --k_components 3 \
  --subset_trait tissue \
  --subset_trait_val liver \
  --no_sex_chr
```

---

## Module 2: Differential Testing

### Basic Usage

```bash
Rscript R/magicTest.R \
  --input data.qs \
  --group condition \
  --compare case,control \
  --mixture_model fitted_model_dir/
```

### Complete Parameter Reference

#### Required Parameters

**`-i, --input <file>`**
- BSseq data file (`.qs` format)
- Same format as optimization module

**`-g, --group <string>`**
- Column name in sample metadata
- Defines groups for comparison
- Must exist in `colData(BSseq)`

**`--mixture_model <dir>`**
- Path to fitted mixture model directory
- Output from optimization module
- Must contain `optModel.csv` or `methComps.csv`

#### Comparison Parameters

**`-c, --compare <string>`**
- Groups to compare (comma-separated)
- Format: "group1,group2"
- Example: `--compare case,control`
- Omit for correlation analysis

#### Filtering Parameters

**`-r, --min_reads <int>`** (default: 5)
- Minimum read coverage
- Same as optimization module

**`-m, --min_samp <int>`** (default: 3)
- Minimum samples per group
- Sites with fewer samples excluded

**`--min_var <float>`** (default: 0.01)
- Minimum methylation variance

**`-f, --min_cpg_fraction <float>`** (default: 0.8)
- For tiled data: minimum CpG fraction
- Filters tiles with low CpG density

**`--chr <string>`**
- Restrict to specific chromosome
- Example: `--chr chr1`

**`-s, --subset <string>`** (formerly `--subset_group`)
- Subset samples by metadata column

**`-k, --keep <string>`** (formerly `--subset_group_type`)
- Value to retain for subset filter

#### Statistical Parameters

**`--methods <string>`** (default: "all")
- Statistical methods to use
- Options: "mixture", "magic", "all"
- Comma-separated for multiple
- Example: `--methods mixture,magic`

**`--prior_strength <float>`** (default: 1.0)
- Dirichlet prior concentration
- Controls prior strength in Bayesian test
- Higher values = stronger prior
- Typical range: 0.5-2.0

#### Data Subsetting

**`--data_subset <string>`** (default: "full")
- Subset data for testing
- Options: "tiny" (1K), "small" (10K), "medium" (50K), "large" (100K), "full"
- Useful for development/testing

#### Computational Parameters

**`-n, --cores <int>`** (default: 1)
- Number of CPU cores for parallel processing
- Recommended: 4-8 for typical datasets

#### Other Parameters

**`-z, --randomize`**
- Randomize group assignments
- For null distribution testing
- Useful for FDR calibration

**`-q, --quiet`**
- Suppress progress messages

### Output Structure

Creates directory: `[analysis]_magic_results_[timestamp]/`

Example: `case_vs_control_magic_results_20250110_143022/`

**Output Files:**

1. **results.csv** - Main results table
   ```
   chr,start,end,mean1,mean2,difference,
   bfProb,mixtureP,domComp1,domComp2,entropy1,entropy2
   ```
   
   Columns:
   - Genomic coordinates (chr, start, end)
   - Mean methylation by group (mean1, mean2)
   - Effect size (difference)
   - Bayes factor probability (bfProb)
   - Mixture p-value (mixtureP)
   - Dominant components (domComp1, domComp2)
   - Shannon entropy (entropy1, entropy2)

2. **magic_summary.txt** - Analysis summary
   - Dataset characteristics
   - Filtering parameters
   - Method descriptions
   - Significance thresholds

### Interpreting Results

#### Bayes Factor Probability (bfProb)

- Range: 0-1
- **>0.9**: Strong evidence for difference
- **0.5-0.9**: Moderate evidence
- **<0.5**: Weak/no evidence
- Null value: 0.5 (no difference)

#### Mixture P-value (mixtureP)

- Standard p-value from Wald test
- Adjust for multiple testing (e.g., Benjamini-Hochberg)
- Typical threshold: 0.05 after FDR correction

#### Effect Size (difference)

- Mean methylation difference between groups
- Range: -1 to 1
- **>0.2**: Moderate effect
- **>0.5**: Large effect

#### Component Assignment (domComp1, domComp2)

- Dominant mixture component per group
- Indicates methylation state:
  - Component 1: Typically unmethylated
  - Component 2: Partially methylated
  - Component 3: Highly methylated
- Different components → different states

#### Entropy

- Shannon entropy of component assignment
- Range: 0 (certain) to log(K) (uniform)
- Lower entropy → clearer assignment
- Higher entropy → ambiguous state

### Statistical Methods

#### "mixture" Method
- Mixture-weighted dispersion estimation
- Component-specific parameters
- Accounts for heterogeneity

#### "magic" Method
- Bayesian framework with mixture priors
- Incorporates component structure
- Provides Bayes factors

### Common Usage Patterns

#### Standard Case-Control
```bash
Rscript R/magicTest.R \
  --input BSseq_data.qs \
  --group condition \
  --compare case,control \
  --mixture_model results/data_K3/ \
  --cores 8
```

#### Subset Analysis
```bash
Rscript R/magicTest.R \
  --input BSseq_data.qs \
  --group condition \
  --compare case,control \
  --mixture_model results/data_K3/ \
  --subset tissue \
  --keep liver \
  --cores 4
```

#### Chromosome-Specific
```bash
Rscript R/magicTest.R \
  --input BSseq_data.qs \
  --group condition \
  --compare case,control \
  --mixture_model results/data_K3/ \
  --chr chr21 \
  --cores 4
```

#### Null Testing
```bash
Rscript R/magicTest.R \
  --input BSseq_data.qs \
  --group condition \
  --compare case,control \
  --mixture_model results/data_K3/ \
  --randomize \
  --cores 8
```

---

## Module 3: Model Comparison

For simulation studies, compare fitted vs true parameters.

### Basic Usage

```bash
Rscript scripts/model_comparison.R \
  -t simulations/true_models/ \
  -f results/fitted_models/ \
  -o comparison.csv \
  -p plots.pdf
```

### Parameters

**`-t, --true-dir <dir>`** (required)
- Directory containing true model files
- Files named: `*_trueModel.csv`

**`-f, --fitted-dir <dir>`** (required)
- Directory containing fitted model subdirectories
- Subdirectories named: `*_K[n]_*`

**`-o, --output <file>`**
- Save comparison summary to CSV

**`-p, --plots <file>`**
- Generate performance plots (PDF)

**`-h, --help`**
- Show help message

### Output

**CSV file** contains:
- RMSE for α, β, π parameters
- Mean absolute relative errors
- Per-component metrics

**PDF file** contains:
- RMSE vs K plots
- Parameter recovery plots
- Convergence stability plots

---

## Tips and Best Practices

### Data Preparation

1. **Quality Control**
   - Remove low-coverage samples
   - Filter sex chromosomes if needed
   - Check for batch effects

2. **Coverage Thresholds**
   - Balance stringency with sample size
   - Typical: 5-10 reads minimum
   - Higher for publication-quality analyses

3. **Sample Size**
   - Minimum 10-15 samples per group
   - More samples → better parameter estimates
   - Consider power analysis

### Model Selection

1. **Start Simple**
   - Begin with K=2-3
   - Increase K if poor fit

2. **Compare Models**
   - Use BIC for model selection
   - Consider interpretability
   - Check convergence stability

3. **Validate**
   - Use holdout validation
   - Check multiple convergence runs
   - Examine parameter estimates

### Performance Optimization

1. **Threading**
   - Use physical cores, not hyperthreads
   - Test optimal thread count
   - More threads ≠ always faster

2. **Memory**
   - Process chromosomes separately if needed
   - Use stringent filtering for large datasets
   - Monitor memory usage

3. **Batch Processing**
   - Process multiple chromosomes in parallel
   - Use job scheduler for HPC
   - Save intermediate results

### Troubleshooting

**Convergence Issues**
- Increase `--max_iters`
- Decrease `--tol`
- Try more `--n_conv_runs`
- Check data quality

**Memory Issues**
- Use `--chr` to process one chromosome
- Increase `--min_samples` threshold
- Reduce `--threads`

**Slow Performance**
- Check thread scaling
- Use SSD storage
- Filter more aggressively
- Process in chunks

---

## Example Workflows

### Complete Analysis Pipeline

```bash
#!/bin/bash

# 1. Optimize mixture model
Rscript R/magicFit.R \
  --input data/BSseq_cohort.qs \
  --k_components 3 \
  --n_conv_runs 10 \
  --threads 8 \
  --holdout

# 2. Test for differential methylation
Rscript R/magicTest.R \
  --input data/BSseq_cohort.qs \
  --group disease_status \
  --compare case,control \
  --mixture_model results/BSseq_cohort_K3_S10_e06_it500/ \
  --cores 8

# 3. Extract significant sites
awk -F',' '$7 < 0.05' results/case_vs_control_magic_results_*/results.csv \
  > significant_sites.csv
```

### Chromosome-by-Chromosome

```bash
#!/bin/bash

for chr in {1..22}; do
  echo "Processing chr$chr"
  
  Rscript R/magicFit.R \
    --input data.qs \
    --k_components 3 \
    --chr chr$chr \
    --threads 4 &
done

wait  # Wait for all to complete
```

### Systematic Model Selection

```bash
#!/bin/bash

for k in 2 3 4 5; do
  Rscript R/magicFit.R \
    --input data.qs \
    --k_components $k \
    --n_conv_runs 10 \
    --threads 8
done

# Compare BIC values
find results/ -name "runSum.csv" -exec \
  awk -F',' 'NR>1 {print $1,$3}' {} \;
```

---

## Getting Help

For additional help:
- Check examples in `examples/` directory
- See method documentation in `docs/method.md`
- Open an issue: https://github.com/Pickledzebra/MAGIC/issues
