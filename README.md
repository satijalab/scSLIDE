# scSLIDE

**Single Cell Sketching and Landmark Integrated Dimensional Embedding**

scSLIDE is an R package that extends Seurat with advanced functionality for single-cell RNA sequencing analysis. It provides enhanced methods for dimensionality reduction, sketching, trajectory analysis, and sample-level aggregation for large-scale single-cell datasets.

## Features

### Dimensionality Reduction
- **RunPLS**: Partial Least Squares (PLS) dimensionality reduction with support for plsr, spls, and cppls methods
- **RunDiffusionMap**: Diffusion map analysis for trajectory inference

### Enhanced Sketching
- **SketchDataByGroup**: Group-aware sketching that maintains representation across cell types and conditions
- **FindmmNN**: Multi-modal nearest neighbor finding for integrated analysis

### Sample-Level Analysis
- **PrepareSampleObject**: Comprehensive workflow for preparing single-cell data for sample-level analysis
- **GenerateSampleObject**: Generate sample-level count matrices from single-cell data
- **NormalizeChiSquared**: Chi-squared normalization for sample-level data

### Trajectory Analysis
- **TrajDETest**: Trajectory-based differential expression analysis using negative binomial regression
- **QuickCorTest**: Fast correlation testing between genes and response variables

### Visualization
- **BuildLandmarkObject**: Build landmark objects with correlation analysis and UMAP embedding
- **PlotLandmarkObject**: Visualize landmark-trajectory correlations
- **SampleLevelDimPlot**: Sample-level visualization combining correlation and density plots
- **RunAndProjectUMAP**: Generate UMAP for sketched data and project to full dataset

## Installation

```r
# Install from GitHub (when available)
# devtools::install_github("yourusername/scSLIDE")

# For now, install dependencies
install.packages(c("pls", "spls", "glmGamPoi", "edgeR", "destiny", 
                   "ggplot2", "dplyr", "tidyr", "RColorBrewer"))

# Install Seurat and SeuratObject
install.packages("Seurat")
```

## Quick Start

```r
library(scSLIDE)
library(Seurat)

# Load your Seurat object
seurat_obj <- your_seurat_object

# Run PLS dimensionality reduction
seurat_obj <- RunPLS(seurat_obj, Y = "condition", ncomp = 10)

# Enhanced sketching by group
sketched_obj <- SketchDataByGroup(seurat_obj, 
                                  group.by = "cell_type", 
                                  ncells = 500)

# Trajectory analysis
de_results <- TrajDETest(seurat_obj, traj.var = "pseudotime")

# Sample-level analysis workflow
prepared_obj <- PrepareSampleObject(seurat_obj, 
                                    Y = "condition",
                                    group.by.Sketch = "cell_type")

sample_obj <- GenerateSampleObject(prepared_obj, 
                                   group.by = "donor_id")
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