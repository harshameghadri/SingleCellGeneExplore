
# SingleCellGeneExplore

<!-- badges: start -->
[![R-CMD-check](https://github.com/harshameghadri/SingleCellGeneExplore/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/harshameghadri/SingleCellGeneExplore/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R >= 4.1.0](https://img.shields.io/badge/R-≥%204.1.0-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

## Overview

SingleCellGeneExplore is a comprehensive R package designed for single-cell RNA sequencing (scRNA-seq) data analysis. This package provides a collection of robust, publication-ready functions for data processing, differential expression analysis, pathway enrichment analysis, and advanced visualization.

## Key Features

### Core Analysis Functions
- **Differential Expression Analysis**: Advanced Wilcoxon rank-sum tests with multiple correction methods
- **Pseudobulk Analysis**: Comprehensive workflow for pseudobulk differential expression
- **Cell Cycle Analysis**: Tools for cell cycle phase assignment and visualization
- **Marker Gene Discovery**: Identification and visualization of cell type-specific markers

### Pathway and Enrichment Analysis
- **Multi-platform Enrichment**: Integration with Enrichr, Gene Ontology, KEGG, and Reactome
- **Pathway Visualization**: Advanced pathview integration for KEGG pathway mapping
- **Parallel Processing**: Optimized enrichment analysis with parallel computing support

### Advanced Visualization
- **Volcano Plots**: Publication-ready differential expression volcano plots with caching
- **Heatmaps**: Unity-normalized and Z-score heatmaps with ComplexHeatmap integration
- **Violin Plots**: Split violin plots for comparative analysis
- **Dot Plots**: Faceted dot plots for multi-condition comparisons
- **Box Plots**: Cell frequency and gene expression box plots

### Data Processing Utilities
- **Data Splitting**: Efficient cell type-based data partitioning
- **Tibble Processing**: Batch processing of analysis results with progress tracking
- **File Management**: Automated result saving and organization

## Installation

You can install the development version of SingleCellGeneExplore from GitHub:

```r
# Install devtools if you haven't already
if (!require(devtools)) install.packages("devtools")

# Install SingleCellGeneExplore
devtools::install_github("harshameghadri/SingleCellGeneExplore")
```

## Dependencies

### Required R Packages
The package depends on several key R packages:

- **Core**: Seurat (≥ 4.0.0), SeuratObject, Matrix
- **Data Manipulation**: dplyr (≥ 1.0.0), tidyr, rlang
- **Visualization**: ggplot2 (≥ 3.3.0), ComplexHeatmap, patchwork
- **Statistical Analysis**: broom, rstatix, epitools
- **Utilities**: progress, parallel, scales, viridis

### Optional Dependencies (Suggested)
For full functionality, consider installing:

```r
# Enrichment analysis
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "pathview"))
install.packages("enrichR")

# Enhanced visualization
install.packages(c("EnhancedVolcano", "cowplot"))

# Parallel processing
install.packages("future.apply")
```

## Quick Start

```r
library(SingleCellGeneExplore)

# Load your Seurat object
# seurat_obj <- readRDS("path/to/your/seurat_object.rds")

# 1. Perform differential expression analysis
de_results <- test_wilcox_on_average(your_data)

# 2. Generate volcano plots
volcano_plots <- generate_differential_volcano_plot_with_caching(
  seurat_obj = seurat_obj,
  cell_type_to_analyze = c("T_cells", "B_cells")
)

# 3. Run comprehensive enrichment analysis
enrichment_results <- perform_enrichment_dual_pathview(
  de_results_list = de_results,
  outdir = "enrichment_output"
)

# 4. Create unity heatmaps
heatmap_result <- generate_unity_heatmap(
  seurat_obj = seurat_obj,
  genes_to_plot = top_genes
)

# 5. Process multiple datasets
results_list <- calculate_wilcox_on_tibble_list(
  tibble_list = your_tibble_list,
  test_wilcox_fn = test_wilcox_on_average
)
```

## Main Functions

### Differential Expression
- `test_wilcox_on_average()`: Wilcoxon test with fold change calculation
- `process_tibbles_calculate_average_wilcox()`: Batch processing of multiple datasets

### Visualization
- `generate_differential_volcano_plot_with_caching()`: Advanced volcano plots
- `generate_unity_heatmap()`: Unity-normalized heatmaps
- `plot_cell_frequency_boxplots()`: Cell frequency analysis
- `plot_faceted_dotplots_cell_types()`: Multi-condition dot plots
- `split_violin_plot()`: Comparative violin plots

### Enrichment Analysis
- `perform_enrichment_dual_pathview()`: Comprehensive pathway analysis
- `find_cluster_markers_with_dor()`: Marker gene identification

### Data Processing
- `split_data_by_celltype()`: Cell type-based data partitioning
- `run_unified_pseudobulk_workflow()`: End-to-end pseudobulk analysis
- `calculate_wilcox_on_tibble_list()`: Batch tibble processing

## Backward Compatibility

This package maintains backward compatibility with previous function names:

```r
# New (recommended)
results <- test_wilcox_on_average(data)

# Old (still works)
results <- Test.wilcox.onAverage(data)
```

## Documentation

Comprehensive documentation is available for all functions:

```r
# View function help
?test_wilcox_on_average
?perform_enrichment_dual_pathview

# Browse package vignettes
browseVignettes("SingleCellGeneExplore")
```

## Citation

If you use SingleCellGeneExplore in your research, please cite:

```
Meghadri, S.H. (2025). SingleCellGeneExplore: Comprehensive Tools for Single-Cell 
RNA-seq Analysis. R package version 1.0.0. 
https://github.com/harshameghadri/SingleCellGeneExplore
```

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details on how to contribute to this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Support

- **Issues**: Report bugs or request features via [GitHub Issues](https://github.com/harshameghadri/SingleCellGeneExplore/issues)
- **Documentation**: Comprehensive help available via `?function_name`
- **Examples**: See function documentation and package vignettes

## Development Status

This package is actively maintained and developed. Version 1.0.0 represents a major refactoring for publication readiness with:

- ✅ Standardized function naming (snake_case)
- ✅ Comprehensive documentation
- ✅ Backward compatibility
- ✅ Extensive test coverage
- ✅ Publication-ready visualizations
- ✅ Optimized performance with parallel processing

---

*SingleCellGeneExplore: Making single-cell RNA-seq analysis accessible and reproducible.*
