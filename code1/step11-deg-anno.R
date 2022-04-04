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

levels(Idents(sce))
sce = sce[, Idents(sce) %in% 
              c( "FCGR3A+ Mono", "CD14+ Mono"  )] # CD16

levels(Idents(sce))
markers_df <- FindMarkers(object = sce, 
                          ident.1 = 'FCGR3A+ Mono',
                          ident.2 = 'CD14+ Mono',
                          #logfc.threshold = 0,
                          min.pct = 0.25)
head(markers_df)
# 在seurat V4里面， 是 avg_log2FC ， 但是如果是V3版本包，仍然是 avg_logFC
cg_markers_df=markers_df[abs(markers_df$avg_logFC) >1,]
dim(cg_markers_df)
DoHeatmap(subset(sce, downsample = 150),
          slot = 'counts',
          unique(rownames(cg_markers_df)),size=3) 

?FindMarkers 
deg=FindMarkers(object = sce, 
            ident.1 = 'FCGR3A+ Mono',
            ident.2 = 'CD14+ Mono', 
            test.use='MAST' )
head(deg)
save(deg,file = 'deg-by-MAST-for-mono-2-cluster.Rdata')

# 作业：绘制火山图，热图

gene_up=rownames(deg[deg$avg_logFC>0,])
gene_down=rownames(deg[deg$avg_logFC < 0,])


library(org.Hs.eg.db)
#把SYMBOL改为ENTREZID，这里会损失一部分无法匹配到的
gene_up=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                   keys = gene_up,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))
gene_down=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                     keys = gene_down,
                                                     columns = 'ENTREZID',
                                                     keytype = 'SYMBOL')[,2]))

# 人家的包就是 entrez ID  设计
library(clusterProfiler)
#function里是kegg和GO以及GSEA分析
source('functions.R')
#直接利用之前保存的脚本中的代码对正、负相关基因进行kegg富集分析
run_kegg(gene_up,gene_down,pro='test')
# 下面的 run_go 会比较慢，根据需要
# run_go(gene_up,gene_down,pro='shRNA')

#对正相关基因进行富集分析
go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
library(ggplot2)
library(stringr)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")+ 
  ggsave('gene_up_GO_all_barplot.png') 
go=DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
write.csv( go@result,file = 'gene_up_GO_all_barplot.csv')

#对负相关基因进行富集分析
go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")+ 
  ggsave('gene_down_GO_all_barplot.png') 
go=DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
write.csv( go@result,file = 'gene_down_GO_all_barplot.csv')



