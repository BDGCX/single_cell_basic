rm(list = ls())  
library(SCENIC)
packageVersion("SCENIC")
library(Seurat)
library(gplots)
library(ggplot2) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA) # BiocManager::install('GSVA')
library(GSEABase)

library(Seurat)
library(gplots)
library(ggplot2)


scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicOptions

scenicOptions@inputDatasetInfo

library(SCopeLoomR)
scenicLoomPath='output/scenic.loom'
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
exprMat_log[1:4,1:4] 
dim(exprMat_log)

regulons_incidMat <- get_regulons(loom, column.attr.name="MotifRegulons")
regulons <- regulonsToGeneLists(regulons_incidMat)

regulons

gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, regulons, names(regulons)))
gsc
geneset <- gsc

# method=c("gsva", "ssgsea", "zscore", "plage"),
# kcdf=c("Gaussian", "Poisson", "none"),
# mx.diff=FALSE: ES is calculated as the maximum distance of the random walk from 0. 
# mx.diff=TRUE (default): ES is calculated as the magnitude difference between the largest positive and negative random walk deviations.

es.max <- gsva( exprMat_log , geneset, 
                mx.diff=FALSE, verbose=FALSE, 
                parallel.sz=8)
 
es.max[1:4,1:4]
pheatmap::pheatmap(es.max)
load(file = 'sce-monocyte.Rdata')
sce 
ac=data.frame(group= as.character( Idents(sce))) 
rownames(ac)= colnames( es.max )
mat=es.max
mat[1:4,1:4]
pheatmap::pheatmap(mat,
                   show_rownames = T,
                   show_colnames =  F,
                   annotation_col = ac)
pos=order(ac$group)
pheatmap::pheatmap(mat[,pos],
                   show_rownames = T,
                   cluster_cols = F,
                   show_colnames =  F,
                   annotation_col = ac[pos,,drop=F])
 
pheatmap::pheatmap(mat,
                   show_rownames = T,show_colnames = F,
                   annotation_col = ac,
                   filename = 'gsva_tfs.pdf',
                   width = 10,height = 10
)
dev.off()



library(ggpubr)
es.max[1:4,1:4]
head(ac)
boxplot(es.max['CEBPB_extended',] ~ ac$group)
df=data.frame(value=es.max['CEBPB_extended',],
              group= ac$group)
ggboxplot(df, "group", "value",
          color = "group", #palette =c("#00AFBB", "#E7B800", "#FC4E07"),
          add = "jitter", shape = "group")
ggsave('gsea_value_for_CEBPB_extended.pdf',
       width = 4,height = 3)

# 作业：
# https://mp.weixin.qq.com/s/s25DLc-tj0lPAcsPurn89Q
# https://mp.weixin.qq.com/s/vlfcMw2UOeaFXkcwfxXYcQ
# 对这个pbmc数据集里面的 DC 和Platelet做同样的分析，monocle和SCENIC




