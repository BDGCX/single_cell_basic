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

# 参考： https://mp.weixin.qq.com/s/enGx9_Sv5wKLdtygL7b4Jw 

features= c('IL7R', 'CCR7','CD14', 'LYZ',  'IL7R', 'S100A4',"MS4A1", "CD8A",'FOXP3',
            'FCGR3A', 'MS4A7', 'GNLY', 'NKG7',
            'FCER1A', 'CST3','PPBP')
DotPlot(sce, features = unique(features)) + RotatedAxis()

if (file.exists('sce.markers.all_10_celltype.Rdata')) {
  load('sce.markers.all_10_celltype.Rdata')
}else {
  sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE,
                                min.pct = 0.25, 
                                thresh.use = 0.25)
  save(sce.markers,file = 'sce.markers.all_10_celltype.Rdata')
  
}


DT::datatable(sce.markers) 
library(dplyr) 
# 不同seurat版本的 avg_logFC 不一样 
top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DoHeatmap(sce,top5$gene,size=3)  
p <- DotPlot(sce, features = unique(top5$gene) ,
             assay='RNA' )  + coord_flip()

p+ theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))

head(top5)
top5=top5[!duplicated(top5$gene),]
select_genes_all=split(top5$gene,top5$cluster)
select_genes_all
DotPlot(object = sce, 
        features=select_genes_all, 
        assay = "RNA") + theme(axis.text.x = element_text(angle = 45, 
                                                          vjust = 0.5, hjust=0.5))

# 作业： https://mp.weixin.qq.com/s/-17oUL0-GProZb9apiZJkg
# 文章是《High-Throughput Single-Cell Transcriptome Profiling of Plant Cell Types》，里面的图 
