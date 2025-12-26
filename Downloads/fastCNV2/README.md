# fastCNV2

## Introduction

`fastCNV2` is an enhanced version of fastCNV that helps you detect, plot and analyse the putative Copy Number Variations (CNVs) in single cell (scRNA-seq) data or Spatial Transcriptomics (ST) data, including Visium HD. Built on `SeuratObject`, it is easily integrated into scRNA-seq or ST pipelines.

**Key improvements in fastCNV2:**
- Enhanced compatibility with both Seurat 4.x and 5.x
- Replaced SCTransform with standard normalization workflow for better stability
- Improved error handling and user experience
- Fully automated dependency management

## Installation

To install fastCNV2, run in R:

```r
# Install from GitHub
remotes::install_github("Zaoqu-Liu/fastCNV2")
```

## Quick Start

```r
library(fastCNV2)

# Run CNV analysis
result <- fastCNV(
  seuratObj = your_seurat_object,
  sampleName = "Sample1",
  referenceVar = "celltype",
  referenceLabel = c("Normal_cells"),
  prepareCounts = TRUE,
  getCNVPerChromosomeArm = TRUE,
  getCNVClusters = TRUE,
  doPlot = TRUE
)

# Check results
table(result$cnv_clusters)
```

## Usage with CellScope

fastCNV2 is designed to work seamlessly with the CellScope package:

```r
library(CellScope)

# Load data
data("pancreas_sub")

# Run fast CNV detection
result <- RunFastCNV(
  srt = pancreas_sub,
  reference_celltype = "CellType",
  reference_labels = "Ductal",
  window_size = 150,
  generate_plot = FALSE,
  verbose = TRUE
)

# View results
head(result@meta.data)
table(result$cnv_clusters)
```

## Features

### 1. Data Types Supported
- **scRNA-seq**: Single-cell RNA sequencing data
- **Visium**: Spatial transcriptomics data
- **Visium HD**: High-definition spatial transcriptomics data

### 2. Key Functions
- `fastCNV()`: Main function for CNV analysis
- `fastCNV_10XHD()`: Specialized function for Visium HD data
- `CNVClassification()`: Cell classification based on CNV patterns
- `plotCNVResults()`: Visualization of CNV heatmaps

### 3. Performance
- ~1 minute for small datasets (~4,000 cells)
- ~40 minutes for Visium HD slide at 16 µm resolution

## Improvements Over Original fastCNV

1. **Compatibility**: Works with both Seurat 4.x and 5.x
2. **Stability**: Replaced SCTransform with standard normalization
3. **Error Handling**: Better error messages and fallback mechanisms
4. **User Experience**: Automatic dependency management
5. **Maintenance**: Active development and bug fixes

## FAQ

### Does it work on mouse data?
Currently works best on **human data**. Mouse data support is in development.

### Can I run without reference cells?
Yes, but we **highly recommend** using reference cells (normal/healthy cells) for better accuracy.

### What if I have multiple samples?
You can process multiple samples together by providing a list:

```r
fastCNV(
  seuratObj = list(sample1, sample2, sample3),
  sampleName = c("sample1", "sample2", "sample3"),
  referenceVar = "annotations",
  referenceLabel = c("Normal1", "Normal2")
)
```

## Citation

If you use fastCNV2 in your research, please cite the original fastCNV:

```
Cabrejas G, Groeneveld C, et al. fastCNV: Fast and accurate inference of copy number 
variations from single-cell RNA sequencing data. bioRxiv 2025.

Original fastCNV: https://github.com/must-bioinfo/fastCNV
```

This is a modified version with enhanced Seurat compatibility maintained at:
https://github.com/Zaoqu-Liu/fastCNV2

## License

GPL-3

## Authors

- **Gadea Cabrejas** (Original Author) - gadea.cabrejas-saiz@u-paris.fr
- **Clarice Groeneveld** (Original Author) - clarice.groeneveld@inserm.fr
- **Zaoqu Liu** (Maintainer, Modified Version) - liuzaoqu@163.com

## Links

- GitHub: https://github.com/Zaoqu-Liu/fastCNV2
- Original fastCNV: https://github.com/must-bioinfo/fastCNV
- CellScope: https://github.com/Zaoqu-Liu/CellScope
