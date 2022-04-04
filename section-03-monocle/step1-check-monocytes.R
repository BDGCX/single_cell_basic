## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2021-06-26 16:13:22
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log:  2021-06-26  First version
###
### ---------------

rm(list = ls())
library(Seurat)
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = '../section-01-cluster/basic.sce.pbmc.Rdata')
pbmc

DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
ggsave('DimPlot-umap.pdf',
       units = 'cm', height = 8, width = 16)

sce=pbmc 
table( Idents(sce ))
table(sce@meta.data$seurat_clusters) 
table(sce@meta.data$orig.ident) 

# 取子集
levels(Idents(sce))
sce = sce[, Idents(sce) %in% 
            c( "FCGR3A+ Mono", "CD14+ Mono"  )] # CD16
sce
levels(Idents(sce))
markers_df <- FindMarkers(object = sce, 
                          ident.1 = 'FCGR3A+ Mono',
                          ident.2 = 'CD14+ Mono',
                          #logfc.threshold = 0,
                          min.pct = 0.25)
head(markers_df)
# 在seurat V4里面， 是 avg_log2FC ， 但是如果是V3版本包，仍然是 avg_logFC
#cg_markers_df=markers_df[abs(markers_df$avg_log2FC) >1,]
cg_markers_df=markers_df[abs(markers_df$avg_logFC) >1,]
dim(cg_markers_df) 
cg_markers_df=cg_markers_df[order(cg_markers_df$avg_logFC),]
DotPlot(sce,
        features = rownames(cg_markers_df)) + theme(axis.text.x = element_text(angle = 45, 
                                                                               vjust = 0.5, hjust=0.5))
ggsave('DotPlot-cg_markers_df-monocyte.pdf')
DoHeatmap(sce,
        features = rownames(cg_markers_df)) 
ggsave('DoHeatmap-cg_markers_df-monocyte.pdf')

save(sce,file = 'sce-monocyte.Rdata')





