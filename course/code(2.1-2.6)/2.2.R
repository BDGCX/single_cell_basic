### 2、构建seurat对象，质控绘图----
# 2.1 构建seurat对象，质控
#In total, 2,343 cells from tumor cores were included in this analysis.
#quality controlstandards: 
#1) genes detected in < 3 cells were excluded; 筛选基因
#2) cells with < 50 total detected genes were excluded; 筛选细胞 
#3) cells with ≥ 5% of mitochondria-expressed genes were excluded. 筛选细胞
sessionInfo()
library("Seurat")
?CreateSeuratObject
sce.meta <- data.frame(Patient_ID=group$Patient_ID,
                       row.names = group$sample)
head(sce.meta)
table(sce.meta$Patient_ID)
# 这个函数 CreateSeuratObject 有多种多样的执行方式
scRNA = CreateSeuratObject(counts=a.filt,
                           meta.data = sce.meta,
                           min.cells = 3, 
                           min.features = 50)
#counts:a matrix-like object with unnormalized data with cells as columns and features as rows 
#meta.data:Additional cell-level metadata to add to the Seurat object
#min.cells: features detected in at least this many cells. 
#min.features:cells where at least this many features are detected.
head(scRNA@meta.data)
#nCount_RNA：the number of cell total counts
#nFeature_RNA：the number of cell's detected gene
summary(scRNA@meta.data)
scRNA@assays$RNA@counts[1:4,1:4]
# 可以看到，之前的counts矩阵存储格式发生了变化：4 x 4 sparse Matrix of class "dgCMatrix"

dim(scRNA)
#  20047  2342 仅过滤掉一个细胞
#接下来根据线粒体基因表达筛选低质量细胞
#Calculate the proportion of transcripts mapping to mitochondrial genes
table(grepl("^MT-",rownames(scRNA)))
#FALSE 


#20050 没有染色体基因
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
head(scRNA@meta.data)
summary(scRNA@meta.data)
#结果显示没有线粒体基因，因此这里过滤也就没有意义，但是代码留在这里
# 万一大家的数据里面有线粒体基因，就可以如此这般进行过滤啦。
pctMT=5 #≥ 5% of mitochondria-expressed genes
scRNA <- subset(scRNA, subset = percent.mt < pctMT)
dim(scRNA)


table(grepl("^ERCC-",rownames(scRNA)))
#FALSE  TRUE 
#19961    86  发现是有ERCC基因
#External RNA Control Consortium，是常见的已知浓度的外源RNA分子spike-in的一种
#指标含义类似线粒体含量，ERCC含量大，则说明total sum变小
scRNA[["percent.ERCC"]] <- PercentageFeatureSet(scRNA, pattern = "^ERCC-")
head(scRNA@meta.data)
summary(scRNA@meta.data)
rownames(scRNA)[grep("^ERCC-",rownames(scRNA))]
#可以看到有不少ERCC基因

sum(scRNA$percent.ERCC< 40)
#较接近原文过滤数量2149，但感觉条件有点宽松了，先做下去看看
#网上看了相关教程，一般ERCC占比不高于10%
sum(scRNA$percent.ERCC< 10)   #就只剩下460个cell，明显低于文献中的数量
pctERCC=40
scRNA <- subset(scRNA, subset = percent.ERCC < pctERCC)
dim(scRNA)
# 20047  2142   原文为19752  2149
dim(a.filt)
#23460  2343 未过滤前


# 2.2 可视化
#图A：观察不同组cell的counts、feature分布
col.num <- length(unique(scRNA@meta.data$Patient_ID))
library(ggplot2)

p1_1.1 <- VlnPlot(scRNA,
                features = c("nFeature_RNA"),
                group.by = "Patient_ID",
                cols =rainbow(col.num)) +
  theme(legend.position = "none") +
  labs(tag = "A")
p1_1.1
p1_1.2 <- VlnPlot(scRNA,
                features = c("nCount_RNA"),
                group.by = "Patient_ID",
                cols =rainbow(col.num)) +
  theme(legend.position = "none") 
p1_1.2
p1_1 <- p1_1.1 | p1_1.2
p1_1
VlnPlot(scRNA,
        features = c("nFeature_RNA","nCount_RNA","percent.ERCC"))
#图B：nCount_RNA与对应的nFeature_RNA关系
p1_2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                       group.by = "Patient_ID",pt.size = 1.3) +
  labs(tag = "B")
p1_2
FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ERCC")
sessionInfo()


