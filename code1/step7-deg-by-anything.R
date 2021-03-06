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
load(file = 'basic.sce.pbmc.Rdata')

DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
sce=pbmc 

if (file.exists('sce.markers.all_10_celltype.Rdata')) {
  load('sce.markers.all_10_celltype.Rdata')
}else {
  sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE,
                                min.pct = 0.25, 
                                thresh.use = 0.25)
  save(sce.markers,file = 'sce.markers.all_10_celltype.Rdata')
  
}

head(sce.markers)
table(sce.markers$cluster)
# 首先挑选基因
kp=grepl('Mono',sce.markers$cluster)
table(kp)
cg_sce.markers = sce.markers [ kp ,]

# 然后挑选细胞
kp=grepl('Mono', Idents(sce ) )
table(kp)
sce=sce[,kp]
sce
table( Idents(sce ))
# 在seurat V4里面， 是 avg_log2FC ， 但是如果是V3版本包，仍然是 avg_logFC
cg_sce.markers=cg_sce.markers[cg_sce.markers$avg_logFC>2,]
dim(cg_sce.markers)
DoHeatmap(subset(sce, downsample = 15),
          unique(cg_sce.markers$gene),
          slot = 'counts',
          size=3 ) 

levels(Idents(sce))
markers_df <- FindMarkers(object = sce, 
                          ident.1 = 'FCGR3A+ Mono',
                          ident.2 = 'CD14+ Mono',
                          #logfc.threshold = 0,
                          min.pct = 0.25)
head(markers_df)
# 在seurat V4里面， 是 avg_log2FC ， 但是如果是V3版本包，仍然是 avg_logFC
cg_markers_df=markers_df[abs(markers_df$avg_log2FC) >1,]
dim(cg_markers_df)
DoHeatmap(subset(sce, downsample = 15),
          slot = 'counts',
          unique(rownames(cg_markers_df)),size=3) 

intersect( rownames(cg_markers_df) ,
          cg_sce.markers$gene)


# drop-out
#想自己划分一个亚群，例如某个基因高表达是一个亚群
highCells= colnames(subset(x = sce, subset = FCGR3A > 1,
                           slot = 'counts')) 
highORlow=ifelse(colnames(sce) %in% highCells,'high','low')
table(highORlow)
table(sce.markers$cluster)
table(Idents(sce))
table(Idents(sce),highORlow)
sce@meta.data$highORlow=highORlow


markers <- FindMarkers(sce, ident.1 = "high", 
                       group.by = 'highORlow' )
head(x = markers)
# 在seurat V4里面， 是 avg_log2FC ， 但是如果是V3版本包，仍然是 avg_logFC
cg_markers=markers[abs(markers$avg_log2FC) >1,]
dim(cg_markers)
DoHeatmap(subset(sce, downsample = 15),
          rownames(cg_markers) ,
          size=3 ) 

DoHeatmap(subset(sce, downsample = 15),
          unique(cg_sce.markers$gene),
          size=3 ) 
DoHeatmap(subset(sce, downsample = 15),
          unique(rownames(cg_markers_df)),size=3) 

 

