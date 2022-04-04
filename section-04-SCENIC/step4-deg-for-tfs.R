rm(list = ls())  
library(SCENIC)
packageVersion("SCENIC")

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
regulonAUC <- get_regulons_AUC(loom)
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom) 
#cellClusters <- get_clusterings(loom)
close_loom(loom)

rownames(regulonAUC)
names(regulons)

head(names(regulons))
regulons[[1]]
load(file = 'sce-monocyte.Rdata')
sce 

DimPlot(sce, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()

# 检查任意一个转录因子 
sg_list=lapply(regulons, function(x) x[x%in% rownames(sce)])
sg_list
th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
DotPlot(object = sce, 
        features=sg_list[[28]] ,  # [28] "HES1"
        assay = "RNA") + th
 
length(sg_list[[28]])
# 30 , "HES1 (16g)"  
getAUC(regulonAUC[30,])
hist(getAUC(regulonAUC[30,]))
boxplot( as.numeric(getAUC(regulonAUC[30,]) ) ~ 
           as.character( Idents(sce)))


# [7] "JUN_extended"
DotPlot(object = sce, 
        features=sg_list[[7]],  
        assay = "RNA") + th
getAUC(regulonAUC[1,])
hist(getAUC(regulonAUC[1,]))
boxplot( as.numeric(getAUC(regulonAUC[1,]) ) ~ 
           as.character( Idents(sce)))
library(ggpubr)
df=data.frame(value=as.numeric(getAUC(regulonAUC[1,]))  ,
              group= as.character( Idents(sce)))
ggboxplot(df, "group", "value",
          color = "group", 
          add = "jitter", shape = "group")+ stat_compare_means(method = "t.test")


library(pheatmap)
pheatmap(  getAUC(regulonAUC[,] ),show_colnames = F)
n=t(scale(t( getAUC(regulonAUC[,] )))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(group= as.character( Idents(sce)))
rownames(ac)=colnames(n)
pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col=ac)
pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col=ac,
         filename = 'heatmap_top_regulon.png')

dev.off()




