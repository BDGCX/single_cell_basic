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
levels(Idents(pbmc))
# 首先提取T细胞子集
DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(pbmc, features = c("CD3D","CD3E"))
sce=pbmc
table(Idents(sce))
t_sce = sce[, Idents(sce) %in% 
              c(  'CD8 T' )]
# 然后进行标准的降维聚类分群
# 代码不要变动
sce=t_sce

sce <- NormalizeData(sce, normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst',
                            nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 1 )
# Look at cluster IDs of the first 5 cells
head(Idents(sce), 5)
table(sce$seurat_clusters) 
sce <- RunUMAP(sce, dims = 1:10)
DimPlot(sce, reduction = 'umap')

sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers)
# write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
save(sce.markers,file = 'sce.markers-for-cd8-subsets.Rdata')

library(dplyr) 
# 在seurat V4里面， 是 avg_log2FC ， 但是如果是V3版本包，仍然是 avg_logFC
top10 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DoHeatmap(sce,top10$gene,size=3)  
p <- DotPlot(sce, features = unique(top10$gene),
             assay='RNA' )  + coord_flip()

p
ggsave('check-top5-for-cd8-subsets.pdf' )


# cytotoxicity (GZMB, PRF1),  0 
# naive (LEF1, SELL, TCF7), 1 
# LTB T cells , 2
# RPl/s 

# 参考： https://mp.weixin.qq.com/s/1YjoX8lTIB0e2vV0Eqey2Q 

marker_genes= c("MKI67","TOP2A",'TNFRSF9','MX1',
                "SELL","IL7R","CD40LG","ANXA1","FOS"
                )


p <- DotPlot(sce, features = marker_genes,
             assay='RNA'  )  + coord_flip() +
  theme(axis.text.x = element_text(angle = 90))

p
ggsave('check_g1_markers_by_Tcell-SubType.pdf')

marker_genes= c("MKI67","TOP2A",'TNFRSF9','MX1')
FeaturePlot(pbmc, features = marker_genes )





