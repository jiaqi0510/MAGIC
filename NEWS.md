# MAGIC Release Notes

## Version 1.0.0 (2025-XX-XX)

### Initial Release

First public release of MAGIC (Methylation Analysis with Genomic Inferred Contexts).

#### Features

**Core Functionality:**
- Beta-binomial mixture model parameter estimation
- Maximum likelihood optimization via L-BFGS-B
- Multiple random initializations for global optimization
- Bayesian differential methylation testing
- Component-wise dispersion estimation
- Mixture-based statistical framework

**Computational:**
- OpenMP parallelization for multi-core systems
- Efficient log-space calculations for numerical stability
- Chunked processing for memory efficiency
- Support for genome-scale datasets

**Model Selection:**
- BIC and AIC for model comparison
- Holdout validation support
- Convergence diagnostics
- Multi-run stability assessment

**Testing:**
- Case-control differential testing
- Correlation analysis support
- Bayes factor computation
- Mixture-based p-values
- FDR correction utilities

**Utilities:**
- Model comparison for simulation studies
- Parameter recovery assessment
- RMSE and CV metrics
- Visualization utilities

#### Documentation

- Comprehensive README with installation and usage
- Detailed installation guide for Linux, macOS, Windows
- Complete usage guide with all parameters
- Mathematical method description
- Three working examples:
  - Mixture model optimization
  - Differential methylation testing
  - Model comparison for simulations

#### Platform Support

- Linux (primary platform)
- macOS (with libomp)
- Windows (via Rtools)

#### Dependencies

**R Packages:**
- Rcpp (≥ 1.0.0)
- bsseq, GenomicRanges, matrixStats (Bioconductor)
- optparse, data.table, ggplot2, gridExtra (CRAN)

**System:**
- R ≥ 4.0.0
- C++ compiler with C++11 and OpenMP support

#### Known Limitations

- Requires OpenMP for performance (serial execution not supported)
- Memory usage scales with number of sites and samples
- Large datasets may require chromosome-by-chromosome processing

#### Performance

Typical runtime on 4 cores for ~100K CpG sites:
- Optimization (K=3): 2-5 minutes
- Differential testing: 5-10 minutes

---

## Future Releases

Planned features for future versions:

- Additional statistical methods
- GPU acceleration
- Streaming processing for very large datasets
- Integration with standard methylation pipelines
- R package distribution via Bioconductor

---

## Version Numbering

MAGIC follows semantic versioning (MAJOR.MINOR.PATCH):
- MAJOR: Incompatible API changes
- MINOR: New functionality, backward compatible
- PATCH: Bug fixes, backward compatible

---

## Citation

When using MAGIC, please cite:

[Authors] (2025). MAGIC: Methylation Analysis with Genomic Inferred Contexts. 
[Journal]. Manuscript in preparation.
