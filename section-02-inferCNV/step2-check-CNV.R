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
library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
library(Seurat)
library(infercnv)
library(miscTools)

#  Import inferCNV dendrogram
infercnv.dend <- read.dendrogram(file = "inferCNV_output/infercnv.observations_dendrogram.txt")
# Cut tree 
infercnv.labels <- cutree(infercnv.dend, k = 6, 
                          order_clusters_as_data = FALSE)
table(infercnv.labels)
# Color labels
the_bars <- as.data.frame(tableau_color_pal("Tableau 20")(20)[infercnv.labels])
colnames(the_bars) <- "inferCNV_tree"
the_bars$inferCNV_tree <- as.character(the_bars$inferCNV_tree)

infercnv.dend %>% set("labels",rep("", nobs(infercnv.dend)) )  %>% plot(main="inferCNV dendrogram") %>%
  colored_bars(colors = as.data.frame(the_bars), dend = infercnv.dend, sort_by_labels_order = FALSE, add = T, y_scale=10, y_shift = 0)

infercnv.labels=as.data.frame(infercnv.labels)
groupFiles='groupFiles.txt'   
meta=read.table(groupFiles,sep = '\t')
infercnv.labels$V1=rownames(infercnv.labels)
meta=merge(meta,infercnv.labels,by='V1')
table(meta[,2:3]) 

if( ! file.exists(  "cnv_scores.csv")){
  tmp=read.table("inferCNV_output/infercnv.references.txt", header=T)
  down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
  up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
  oneCopy=up-down
  oneCopy
  a1= down- 2*oneCopy
  a2= down- 1*oneCopy
  down;up
  a3= up +  1*oneCopy
  a4= up + 2*oneCopy 
  
  cnv_table <- read.table("inferCNV_output/infercnv.observations.txt", header=T)
  # Score cells based on their CNV scores 
  # Replicate the table 
  cnv_score_table <- as.matrix(cnv_table)
 
  cnv_score_mat <- as.matrix(cnv_table)
  # Scoring
  cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
  cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
  # Check
  table(cnv_score_table[,1])
  # Replace with score 
  cnv_score_table_pts <- cnv_table
  rm(cnv_score_mat)
  # 
  cnv_score_table_pts[cnv_score_table == "A"] <- 2
  cnv_score_table_pts[cnv_score_table == "B"] <- 1
  cnv_score_table_pts[cnv_score_table == "C"] <- 0
  cnv_score_table_pts[cnv_score_table == "D"] <- 1
  cnv_score_table_pts[cnv_score_table == "E"] <- 2
  cnv_score_table_pts[cnv_score_table == "F"] <- 2
  
  # Scores are stored in “cnv_score_table_pts”. Use colSums to add up scores for each cell and store as vector 
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  colnames(cell_scores_CNV) <- "cnv_score"
  head(cell_scores_CNV)
  write.csv(x = cell_scores_CNV, file = "cnv_scores.csv")
  
} 


# 除去了reference后的走inferCNV的细胞
cell_scores_CNV=read.csv('cnv_scores.csv',row.names = 1)
head(cell_scores_CNV) 


load(file = '../section-01-cluster/basic.sce.pbmc.Rdata') 
sce=pbmc 
phe=sce@meta.data
phe$celltype=Idents(sce)
head(rownames(phe))
head(rownames(cell_scores_CNV)) 

#rownames(phe)=paste0('X',rownames(phe))
rownames(phe)=gsub('-','.',rownames(phe))

head(rownames(phe)) 
head(rownames(cell_scores_CNV))

head(rownames(phe))
phe=phe[rownames(phe) %in% rownames(cell_scores_CNV),]
identical(rownames(phe),rownames(cell_scores_CNV))

infercnv.labels <- cutree(infercnv.dend, k = 6, order_clusters_as_data = FALSE)
phe$inferCNV= infercnv.labels[match(rownames(phe), names(infercnv.labels) )]

phe$cnv_scores  =  cell_scores_CNV[rownames(phe),]

table(phe$celltype,phe$inferCNV) 
head(rownames(phe))
dim(phe)
library(ggpubr)
p1=ggboxplot(phe,'celltype','cnv_scores', fill = "celltype") 
p1= p1+ theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
p2=ggboxplot(phe,'inferCNV','cnv_scores', fill = "inferCNV")
library(patchwork)
p1+p2
ggsave(filename = 'anno_CNVscore.pdf')
table(phe$celltype,phe$inferCNV)


#rownames(phe)=gsub('X','',rownames(phe))
rownames(phe)=gsub('[.]','-',rownames(phe))
head(rownames(phe))

sce
tail(rownames(sce@meta.data))
head(rownames(phe))
sce$celltype=Idents(sce)
table(sce$celltype)
sce=subset(sce,celltype %in% c('epithelial' )) 
sce
kp=rownames(sce@meta.data) %in% rownames(phe)
table(kp)
sce=sce[,kp]

phe=phe[rownames(sce@meta.data),]

sce@meta.data=phe
head(phe)
save(sce, file = 'epi_sce_annoCNV.Rdata')

