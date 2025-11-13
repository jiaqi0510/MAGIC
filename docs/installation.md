# Installation Guide

This document provides detailed installation instructions for MAGIC.

## System Requirements

### Operating Systems

MAGIC has been tested on:
- Linux (Ubuntu 20.04+, CentOS 7+, Debian 10+)
- macOS (10.15+)
- Windows 10/11 (via Rtools)

**Primary Platform**: Linux with GCC and OpenMP

### Software Dependencies

**Required:**
- R ≥ 4.0.0
- C++ compiler with C++11 standard support
- OpenMP library for parallelization

**Recommended:**
- 8+ GB RAM for typical datasets
- 4+ CPU cores for parallel processing
- SSD storage for large datasets

## R Package Dependencies

### CRAN Packages

```r
install.packages(c(
  "Rcpp",        # C++ integration
  "optparse",    # Command-line parsing
  "data.table",  # Fast data manipulation
  "ggplot2",     # Plotting
  "gridExtra"    # Plot arrangements
))
```

### Bioconductor Packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "bsseq",          # BSseq data structures
  "GenomicRanges",  # Genomic intervals
  "matrixStats"     # Matrix operations
))
```

### Optional Packages

```r
install.packages("qs")  # For fast serialization (recommended for large datasets)
```

## Platform-Specific Instructions

### Linux Installation

#### Ubuntu/Debian

```bash
# Install R
sudo apt-get update
sudo apt-get install r-base r-base-dev

# Install compiler and OpenMP
sudo apt-get install build-essential
sudo apt-get install libgomp1

# Clone MAGIC
git clone https://github.com/Pickledzebra/MAGIC.git
cd MAGIC

# Install R packages
Rscript -e "install.packages(c('Rcpp', 'optparse', 'data.table', 'ggplot2', 'gridExtra'))"
Rscript -e "BiocManager::install(c('bsseq', 'GenomicRanges', 'matrixStats'))"

# Test compilation
R CMD SHLIB src/magicFit.cpp -fopenmp
```

#### CentOS/RHEL

```bash
# Install R
sudo yum install epel-release
sudo yum install R R-devel

# Install compiler
sudo yum groupinstall "Development Tools"

# Continue with package installation as above
```

### macOS Installation

macOS requires special handling for OpenMP.

#### Option 1: Install libomp (Recommended)

```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install OpenMP library
brew install libomp

# Set environment variables for R session
export CPPFLAGS="-Xpreprocessor -fopenmp"
export LDFLAGS="-lomp"
export CXXFLAGS="-Xpreprocessor -fopenmp"

# Clone and install
git clone https://github.com/Pickledzebra/MAGIC.git
cd MAGIC

# Install R packages
Rscript -e "install.packages(c('Rcpp', 'optparse', 'data.table', 'ggplot2', 'gridExtra'))"
Rscript -e "BiocManager::install(c('bsseq', 'GenomicRanges', 'matrixStats'))"
```

To make these permanent, add to your `~/.bash_profile` or `~/.zshrc`:
```bash
export CPPFLAGS="-Xpreprocessor -fopenmp"
export LDFLAGS="-lomp"
export CXXFLAGS="-Xpreprocessor -fopenmp"
```

#### Option 2: Install GCC via Homebrew

```bash
# Install GCC (includes OpenMP)
brew install gcc

# Check GCC version
gcc-13 --version  # Version number may vary

# Create ~/.R/Makevars file
mkdir -p ~/.R
cat > ~/.R/Makevars << EOF
CC=gcc-13
CXX=g++-13
CXX11=g++-13
SHLIB_OPENMP_CFLAGS=-fopenmp
SHLIB_OPENMP_CXXFLAGS=-fopenmp
EOF
```

### Windows Installation

```bash
# Install R from CRAN
# Install Rtools from: https://cran.r-project.org/bin/windows/Rtools/

# In R console:
install.packages(c("Rcpp", "optparse", "data.table", "ggplot2", "gridExtra"))
BiocManager::install(c("bsseq", "GenomicRanges", "matrixStats"))

# Clone repository
git clone https://github.com/Pickledzebra/MAGIC.git
cd MAGIC
```

**Note**: Rtools includes GCC with OpenMP support. Ensure Rtools is in your PATH.

## Compilation

### Automatic Compilation

C++ code is compiled automatically when first running the R scripts:

```bash
cd MAGIC
Rscript R/magicFit.R --help
```

Rcpp will compile `src/magicFit.cpp` on first use.

### Manual Compilation

To pre-compile the C++ modules:

```bash
cd MAGIC

# Compile optimization module
R CMD SHLIB src/magicFit.cpp -fopenmp

# Compile testing module
R CMD SHLIB src/magicTest.cpp -fopenmp
```

This creates `.so` (Linux), `.dylib` (macOS), or `.dll` (Windows) files.

### Compilation Flags

The R scripts automatically set:
```r
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")
```

For manual compilation with specific flags:
```bash
R CMD SHLIB src/magicFit.cpp \
  -fopenmp \
  -O3 \
  -march=native
```

## Verification

### Test Installation

```bash
cd MAGIC

# Check that C++ modules can be found
Rscript R/magicFit.R --help

# Should display usage information
```

### Test Compilation

Create a simple test:

```r
library(Rcpp)
sourceCpp("src/magicFit.cpp")

# Should load without errors
# Check that functions are available
exists("compute_ll")  # Should return TRUE
exists("compute_gradients")  # Should return TRUE
```

### Test OpenMP

Check thread count:

```bash
# Set thread count
export OMP_NUM_THREADS=4

# Run a script
Rscript R/magicFit.R --threads 4 [other options]

# Should report using 4 threads
```

## Troubleshooting

### OpenMP Not Found

**Error**: `fatal error: omp.h: No such file or directory`

**Solution (Linux)**:
```bash
sudo apt-get install libgomp1
```

**Solution (macOS)**:
```bash
brew install libomp
export CPPFLAGS="-Xpreprocessor -fopenmp"
export LDFLAGS="-lomp"
```

**Solution (Windows)**: Ensure Rtools is installed and in PATH.

### C++ Compilation Errors

**Error**: `R CMD SHLIB` fails with compiler errors

**Diagnosis**:
```bash
# Check compiler version
g++ --version  # Should be ≥ 4.9

# Check R configuration
R CMD config CXX11
R CMD config CXX11FLAGS
```

**Solution**: Update compiler or specify compatible flags in `~/.R/Makevars`.

### Rcpp Compilation Issues

**Error**: `Rcpp.h: No such file or directory`

**Solution**:
```r
install.packages("Rcpp")
```

**Error**: Version mismatch

**Solution**:
```r
update.packages("Rcpp")
```

### Permission Denied

**Error**: Cannot create `.so` file

**Solution**: Check write permissions in repository directory.

### Missing Bioconductor Packages

**Error**: `bsseq` not found

**Solution**:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("bsseq")
```

### OpenMP Thread Count Issues

**Problem**: Using fewer cores than specified

**Diagnosis**:
```bash
echo $OMP_NUM_THREADS
```

**Solution**:
```bash
export OMP_NUM_THREADS=8
# Or specify in script:
Rscript R/magicFit.R --threads 8
```

### Memory Issues

**Error**: Memory allocation errors with large datasets

**Solution**:
- Use chunked processing (automatic in testing module)
- Increase system memory
- Filter to smaller genomic regions with `--chr` flag
- Use `--data_subset` option in testing module

## Performance Optimization

### Compiler Optimizations

For maximum performance, compile with optimization flags:

```bash
R CMD SHLIB src/magicFit.cpp \
  -fopenmp \
  -O3 \
  -march=native \
  -ffast-math
```

**Warning**: `-ffast-math` may reduce numerical precision. Use with caution.

### Thread Scaling

Test optimal thread count:

```bash
for threads in 1 2 4 8; do
  echo "Testing with $threads threads"
  time Rscript R/magicFit.R \
    --threads $threads \
    [other options]
done
```

Optimal thread count is typically:
- Number of physical cores (not hyperthreads)
- For I/O-bound tasks: fewer threads may be faster

### Memory Management

For large datasets:
- Use `.qs` format (faster than `.rds`)
- Filter early with `--chr` flag
- Monitor memory with `top` or `htop`

## Docker Installation (Advanced)

For reproducible environments:

```dockerfile
FROM rocker/r-ver:4.2.0

RUN apt-get update && apt-get install -y \
    libgomp1 \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('Rcpp', 'optparse', 'data.table', 'ggplot2', 'gridExtra'))"
RUN R -e "BiocManager::install(c('bsseq', 'GenomicRanges', 'matrixStats'))"

COPY . /MAGIC
WORKDIR /MAGIC

RUN R CMD SHLIB src/magicFit.cpp -fopenmp
RUN R CMD SHLIB src/magicTest.cpp -fopenmp
```

## Version Information

To check installed versions:

```r
# R version
R.version.string

# Package versions
packageVersion("Rcpp")
packageVersion("bsseq")

# Compiler version
system("g++ --version")
```

## Getting Help

If installation fails:

1. Check the error message against this troubleshooting guide
2. Verify all dependencies are installed
3. Check compiler and OpenMP availability
4. Open an issue: https://github.com/Pickledzebra/MAGIC/issues

Include in your issue:
- Operating system and version
- R version (`R.version.string`)
- Compiler version (`g++ --version`)
- Complete error message
- Output of `R CMD config CXX11`
