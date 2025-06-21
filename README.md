# Shiny scRNA-seq Explorer

This repository contains a Shiny application for interactive visualization and analysis of single-cell RNA sequencing (scRNA-seq) data using Seurat objects.

## Features

- **Cluster Visualization**
  - Interactive `DimPlot` (t-SNE, UMAP, PCA)
  - Group or split by multiple metadata fields
  - Downloadable plots

- **Gene Expression**
  - `FeaturePlot` for single or multiple genes
  - Split plots by metadata
  - Raster or vector rendering options
  - Downloadable plots

- **Pathway Analysis**
  - Overlay pathway activity on t-SNE plots
  - Generate heatmaps grouped by metadata
  - Downloadable pathway plots

- **Custom Differential Expression**
  - Select custom metadata groupings for DE analysis
  - Volcano plots with automatic gene labeling
  - Downloadable DE tables and plots

## Requirements

- **R version:** >= 4.0
- **Packages:**
  - `shiny`
  - `shinydashboard`
  - `shinycssloaders`
  - `dplyr`
  - `Seurat`
  - `ggplot2`
  - `ggrepel`
  - `viridis`
  - `stringr`

## Files

- `app.R` — Main Shiny app code (UI + server logic)
- `fgsea_sets.RDS` — Precomputed pathway gene sets
- `sample_seuratobj.rds` — Sample seurat object containing single-cell data
