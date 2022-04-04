### 5、聚类，筛选marker基因，可视化----
#5.1 聚类
pc.num=1:20
#基于PCA数据
scRNA <- FindNeighbors(scRNA, dims = pc.num) 
# dims参数，需要指定哪些pc轴用于分析；这里利用上面的分析，选择20
scRNA <- FindClusters(scRNA, resolution = 0.5)#分辨率，0~1之间，值越大分的cluster越多
table(scRNA@meta.data$seurat_clusters)

scRNA = RunTSNE(scRNA, dims = pc.num)
DimPlot(scRNA, reduction = "tsne",label=T)
?RunTSNE
p3_1 <- DimPlot(scRNA, reduction = "tsne",label=T) +
  labs(tag = "E")
p3_1

#5.2 marker gene
#进行差异分析，一般使用标准化数据
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize")
#结果储存在"data"slot里
GetAssayData(scRNA,slot="data",assay="RNA")[1:8,1:4]

#if test.use is "negbinom", "poisson", or "DESeq2", slot will be set to "counts
diff.wilcox = FindAllMarkers(scRNA)##默认使用wilcox方法挑选差异基因，大概4-5min
load("../../tmp/diff.wilcox.Rdata")
head(diff.wilcox)
dim(diff.wilcox)
library(tidyverse)
all.markers = diff.wilcox %>% select(gene, everything()) %>%
  subset(p_val<0.05 & abs(diff.wilcox$avg_log2FC) > 0.5)
#An adjusted P value < 0.05and | log 2 [fold change (FC)] | > 0.5 
#were considered the 2 cutoff criteria for identifying marker genes.
dim(all.markers)
summary(all.markers)
save(all.markers,file = "../../tmp/markergene.Rdata")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(scRNA)) 
top10
length(top10)
length(unique(sort(top10)))

p3_2 <- DoHeatmap(scRNA, features = top10, group.by = "seurat_clusters")
p3_2
p3_1 | p3_2 #下图



