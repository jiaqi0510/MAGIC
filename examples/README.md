# MAGIC Examples

This directory contains working examples demonstrating the main MAGIC workflows.

## ⚠️ Important: Working Directory

**All example scripts must be run from the MAGIC root directory, not from the examples/ directory.**

```bash
# Correct:
cd MAGIC
Rscript examples/example_synthetic_workflow.R

# Incorrect:
cd MAGIC/examples
Rscript example_synthetic_workflow.R  # Will not work!
```

The scripts use relative paths like `R/magicFit.R` and `scripts/...` which require running from the repository root.

---

## Example Scripts

### 0. example_synthetic_workflow.R
**Complete Workflow with Synthetic Data (Recommended Starting Point)**

The easiest way to test MAGIC! Demonstrates the complete pipeline:
- Generate synthetic data with known parameters
- Fit mixture model
- Compare fitted vs true parameters
- Evaluate parameter recovery

**Usage:**
```bash
Rscript examples/example_synthetic_workflow.R
```

**Prerequisites:**
- None! Creates its own test data

**Runtime:** ~3-5 minutes

This example is perfect for:
- Testing your installation
- Understanding the complete workflow
- Validating the method works correctly
- Learning before applying to real data

### 1. example_optimization.R
**Mixture Model Parameter Estimation**

Demonstrates:
- Basic model optimization
- Model selection across different K values
- Convergence diagnostics and stability assessment
- Chromosome-specific analysis
- Holdout validation
- Reading and interpreting fitted parameters

**Usage:**
```bash
Rscript examples/example_optimization.R
```

**Prerequisites:**
- BSseq data file: `data/BSseq_example.qs`

### 2. example_testing.R
**Differential Methylation Testing**

Demonstrates:
- Case-control comparisons
- Interpreting Bayes factors and p-values
- Multiple testing correction (FDR)
- Chromosome-specific testing
- Subset analysis by metadata
- Null distribution testing
- Exporting results for visualization

**Usage:**
```bash
Rscript examples/example_testing.R
```

**Prerequisites:**
- BSseq data file: `data/BSseq_example.qs`
- Fitted mixture model: `results/data_K3_S10_e06_it500/`
  (Run `example_optimization.R` first)

### 3. example_comparison.R
**Model Comparison for Simulations**

Demonstrates:
- Comparing fitted vs true parameters
- RMSE and CV metrics
- Component-wise parameter evaluation
- Analyzing multiple simulation replicates
- Creating publication-quality figures

**Usage:**
```bash
Rscript examples/example_comparison.R
```

**Prerequisites:**
- True model files: `simulations/*_trueModel.csv`
- Fitted model directories: `results/*_K*_*/`

### 4. example_synthetic_with_DML.R
**Model Evaluation with True DML**

Demonstrates:
- Generate synthetic data with controlled DML injection
- Fit MAGIC mixture model
- Test for differential methylation
- Visualize methylation distributions by component
- Create ROC curves comparing methods
- Validate Type I error control with QQ plots

**Usage:**
```bash
Rscript examples/example_synthetic_with_DML.R
```

**Outputs:**
- Synthetic BSseq data with ground truth DML labels
- MAGIC differential methylation results
- Methylation distribution plots by mixture component
- ROC curves (MAGIC vs baseline methods)
- QQ plot for null simulation validation
- Performance metrics (sensitivity, precision, F1)

This example is perfect for:
- Benchmarking MAGIC performance

## Data Requirements

### Creating Synthetic Data

If you don't have real data, generate synthetic test data:

```bash
# Generate K=3 mixture with 100K CpGs and 50 samples
Rscript scripts/generate_synthetic_data.R \
  -K 3 \
  -n 100000 \
  -s 50 \
  -o data/synthetic_K3.qs

# This creates:
# - data/synthetic_K3.qs (BSseq object)
# - data/synthetic_K3_trueModel.csv (true parameters)
```

The generator supports:
- Realistic beta-binomial mixture parameters
- Sparse coverage patterns (like WGBS)
- Sex-specific differential methylation
- Automatic filtering to match magicFit behavior

### BSseq Format

All examples expect BSseq objects saved in `.qs` format with:

```r
library(bsseq)
library(qs)

# Your BSseq object should have:
# - Methylation counts: getCoverage(bs, type = "M")
# - Total coverage: getCoverage(bs, type = "Cov")
# - Sample metadata: colData(bs)

# Save for use with MAGIC:
qsave(bsseq_object, "data/BSseq_example.qs")
```

### True Model Format (for comparisons)

For simulation studies, true model files should be CSV with columns:
- `component`: Component number (1 to K)
- `alpha`: Alpha parameter
- `beta`: Beta parameter
- `pi`: Mixing proportion

Example:
```
component,alpha,beta,pi
1,0.500,12.000,0.300
2,5.000,8.000,0.400
3,15.000,2.000,0.300
```

### Ground Truth DML Table (`_truth_DML.csv`)
Example: `synthetic_bb_dml_100k_s10_d01_cov20_truth_DML.csv`

CSV file containing ground truth information for all injected DML sites:

- `CpG_ID`: Row index of the CpG site in the BSseq object
- `chr`: Chromosome name 
- `pos`: Genomic position 
- `Effect`: Injected effect size 
- `Direction`: Direction of methylation change (increase/decrease) 
- `mean_meth_control`: Mean methylation level in Control group 
- `mean_meth_treatment`: Mean methylation level in Treatment group 
- `min_cov_control`: Minimum coverage across Control samples 
- `min_cov_treatment`: Minimum coverage across Treatment samples 
- `delta`: Observed methylation difference (|treatment - control|) 


## Running Examples

### Quick Start

If you don't have your own data yet, the examples will print helpful error messages explaining what data is needed and how to prepare it.

### With Your Own Data

1. Prepare BSseq object:
   ```r
   library(bsseq)
   library(qs)
   
   # Load your methylation data into BSseq format
   # ... (see bsseq documentation)
   
   # Save for MAGIC
   qsave(bs, "data/my_data.qs")
   ```

2. Edit example scripts to use your data:
   ```r
   # Change this line:
   input_file <- "data/BSseq_example.qs"
   
   # To:
   input_file <- "data/my_data.qs"
   ```

3. Run the examples:
   ```bash
   Rscript examples/example_optimization.R
   Rscript examples/example_testing.R
   ```

## Example Outputs

### From example_optimization.R

Creates directories like:
```
data_K3_S10_e06_it500/
├── optModel.csv           # Fitted parameters
├── methComps.csv          # Interpretable format
├── runSum.csv             # Model fit statistics
├── allParamEst.csv        # All convergence runs
└── allConvRuns.csv        # Convergence diagnostics
```

### From example_testing.R

Creates directories like:
```
case_vs_control_magic_results_20250110_143022/
├── results.csv            # Main results table
├── results_with_fdr.csv   # With FDR correction
├── significant_sites_fdr05.csv
├── magic_summary.txt      # Analysis summary
└── [visualization files]
```

### From example_comparison.R

Creates:
```
comparison_summary.csv     # RMSE metrics
comparison_plots.pdf       # Visualization
parameter_recovery.png     # Custom plots
convergence_stability.png
publication_figure.png     # High-res for papers
publication_figure.pdf     # Vector format
```

### From example_synthetic_with_DML.R

Creates:
```
synthetic_example/
├── synthetic_bb_100k_s10_d01_cov10.qs              # synthetic data with true DMLs
├── synthetic_bb_100k_s10_d01_cov10_truth_DML.csv   # true DMLs
├── synthetic_bb_100k_s10_d01_cov10_trueModel.csv   # model parameters
├── magic_results/
│   └── magic_results.csv
├── methylation_histogram_by_component.png          # visualization of methylation levels
├── roc_curve_comparison.png
├── null_simulation.qs                              # synthetic data with no DMLs
└── qq_plot_null_simulation.png                     # evaluation of false positive control
```

## Customization

All examples are heavily commented and designed to be modified. Key customization points:

- **Input files**: Change file paths to your data
- **Parameters**: Adjust K, convergence settings, filtering
- **Metadata columns**: Update `--group`, `--subset` arguments
- **Output formats**: Modify plotting and export code

## Getting Help

If examples don't work:

1. Check that prerequisites are met
2. Verify data format matches requirements
3. Read error messages carefully
4. See full documentation in `docs/`
5. Open an issue on GitHub

## Additional Resources

- **Installation**: `docs/installation.md`
- **Usage Guide**: `docs/usage.md`
- **Method Description**: `docs/method.md`
- **Main README**: `../README.md`
