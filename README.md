# scSLIDE

**Single-cell Sample-Level Integration using Density Estimation**

scSLIDE is an R package to perform sample-level analysis for multi-sample single-cell RNA sequencing data. It leverages a semi-supervised dimensional reduction framework to embed cells into a latent space that robustly retains both their underlying type- and state-identity as well as phenotype-driven differences. Each sample is then represented as a probability distribution of cellular states, yielding a sample-level representation that can be directly used for clustering, trajectory inference, and integrative analyses.

## Key Features

### Dimensionality Reduction
- **RunPLS**: Partial Least Squares (PLS) dimensionality reduction with support for plsr, spls, and cppls methods
- **RunDiffusionMap**: Diffusion map analysis for trajectory inference 

### Sample-Level Analysis
- **PrepareSampleObject**: Comprehensive workflow for preparing single-cell data for sample-level analysis
- **GenerateSampleObject**: Generate sample-level count matrices from single-cell data

### Novel Differential Expression Test
- **TrajDETest**: Trajectory-based differential expression analysis using negative binomial regression

## Installation

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor dependencies first
BiocManager::install(c("edgeR", "destiny"))

# Install CRAN dependencies
install.packages(c("pls", "spls", "glmGamPoi", "ggplot2", "dplyr", "tidyr", "RColorBrewer"))

# Install Seurat and SeuratObject
install.packages("Seurat")

# Install scSLIDE from GitHub
devtools::install_github("longmanz/scSLIDE")
```

## Dependencies

scSLIDE depends on:
- **Seurat** (>= 5.0.0): Core single-cell analysis framework
- **SeuratObject** (>= 5.0.2): Seurat object structure
- **pls**: Partial least squares regression
- **spls**: Sparse partial least squares
- **glmGamPoi**: Gamma-Poisson regression for DE analysis
- **edgeR**: Gene filtering and preprocessing
- **destiny**: Diffusion map analysis
- **ggplot2, dplyr, tidyr, RColorBrewer**: Visualization and data manipulation

## Citation

If you use scSLIDE in your research, please cite:

```
[Citation information to be added]
```

## License

MIT License - see LICENSE file for details.

## Support

For questions and support, please open an issue on the GitHub repository.
