rm(list = ls())
library(Seurat)
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'downsample-sce.Rdata')
DimPlot(sce,reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(sce,c('SELL','SELPLG'))
VlnPlot(sce,c('SELL','SELPLG'))
VlnPlot(sce,c('CD74','MIF'))
VlnPlot(sce,c('TNFRSF1B','GRN'))


