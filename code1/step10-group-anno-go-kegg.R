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


rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(gplots)
library(ggplot2)
library(clusterProfiler)

pro='all'
load('sce.markers.all_10_celltype.Rdata')
table(sce.markers$cluster)
library(org.Hs.eg.db)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster)
gcSample # entrez id , compareCluster 
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
p=dotplot(xx) 
p+ theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
p

# msigdb， 单细胞相关基因集
source('group_kegg_go.R')
group_kegg_go(gcSample, pro )




