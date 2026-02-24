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

# Install Bioconductor dependencies manually
BiocManager::install(c("glmGamPoi", "destiny"))

# Install SeuratObject 
install.packages("SeuratObject")

# Install Seurat from a developmental branch that is compatible with scSLIDE (built upon v5.3.1) 
remotes::install_github("satijalab/seurat", "v5.3.1_scSLIDE_compatible")

# Install scSLIDE from GitHub
devtools::install_github("satijalab/scSLIDE")
```

## Dependencies

scSLIDE depends on:
- **Seurat** (>= 5.3.1): Core single-cell analysis framework
- **SeuratObject** (>= 5.2.0): Seurat object structure
- **pls**: Partial least squares regression
- **spls**: Sparse partial least squares
- **glmGamPoi**: Gamma-Poisson regression for DE analysis
- **destiny**: Diffusion map analysis
- **ggplot2, dplyr, tidyr, RColorBrewer**: Visualization and data manipulation

## Tutorial
A tutorial for running the scSLIDE package is provided [here](https://satijalab.github.io/scSLIDE/articles/scSLIDE_COVID19.html). 

## Citation

If you use scSLIDE in your research, please cite:

```
https://doi.org/10.64898/2025.12.10.693462
```

## License

MIT License - see LICENSE file for details.

## Support

For questions and support, please open an issue on the GitHub repository.
