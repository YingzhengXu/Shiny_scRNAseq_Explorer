Shiny scRNA-seq Explorer
This repository contains a Shiny application for interactive visualization and analysis of single-cell RNA sequencing (scRNA-seq) data. The app is designed to support exploration of clustering, gene expression, and pathway activity from Seurat objects.

Features
Cluster visualization:
Interactive DimPlot for t-SNE, UMAP, or PCA, with options to group or split by metadata (e.g., clustering resolutions, sample IDs).

Gene expression visualization:
FeaturePlot for single or multiple genes, supporting split views and adjustable rendering (fast raster or high-quality vector).

Pathway analysis:
View pathway activity as t-SNE overlays and heatmaps using custom gene sets (e.g., fgsea results).

Custom differential expression:
Compare selected groups and generate volcano plots and downloadable DE tables.

Requirements
R (>= 4.0)
R packages:

shiny

shinydashboard

shinycssloaders

dplyr

Seurat

ggplot2

ggrepel

viridis

stringr

Files
app.R: Main Shiny app code (UI + server)

fgsea_sets.RDS: RDS file containing gene sets used for pathway visualization

singlet.df.rds: Seurat object for scRNA-seq data
