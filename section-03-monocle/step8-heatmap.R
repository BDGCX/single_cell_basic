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
library(monocle)
library(dplyr)
load(file = 'output_of_monocle.Rdata')
cds=my_cds_subset
phe=pData(cds)
colnames(phe)
plot_cell_trajectory(cds)
counts = Biobase::exprs(cds)
dim(counts)

library(dplyr) 
my_pseudotime_de %>% arrange(qval) %>% head(100) %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene=my_pseudotime_gene[,1]
my_pseudotime_gene


library(pheatmap)
n=t(scale(t( counts[my_pseudotime_gene,] ))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=phe[,c(10,16,17)]
head(ac)
rownames(ac)=colnames(n)
dim(n)
n[1:4,1:4]
pheatmap(n,show_colnames =F,
         show_rownames = F,
         annotation_col=ac)
od=order(ac$Pseudotime)
pheatmap(n[,od],show_colnames =F,
         show_rownames = F,cluster_cols = F,
         annotation_col=ac[od,])





