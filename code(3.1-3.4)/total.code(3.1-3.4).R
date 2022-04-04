### step3.1-3.4（对应下面的1-4） ###
load("../../tmp/scRNA.Rdata")
library(Seurat)
##1、detect special cells----
# empty droplet
# BiocManager::install("DropletUtils")
library(DropletUtils)
e.out <- emptyDrops(GetAssayData(scRNA,slot="counts",assay="RNA"))
#Error in testEmptyDrops(m, lower = lower, ...) : 
#no counts available to estimate the ambient profile
##https://support.bioconductor.org/p/123554/#123562
#如上回答所说，empty droplet往往在第一步就已经过滤掉了，而一般上传到GEO的也都是过滤掉空液滴的。

#double droplet
#https://osca.bioconductor.org/doublet-detection.html
# BiocManager::install("scran")
head(scRNA@meta.data)
library(scran)
#GetAssayData(scRNA,slot="counts",assay="RNA")[1:8,1:4]
?doubletCluster #检查有无double droplet聚在一起的类
db.test <- doubletCluster(GetAssayData(scRNA,slot="counts",assay="RNA"),
                          clusters=scRNA@meta.data$seurat_clusters)
head(db.test)
table(scRNA@meta.data$seurat_clusters)
library(scater)
chosen.doublet <- rownames(db.test)[isOutlier(db.test$N, 
                                              type="lower", log=TRUE)]
chosen.doublet #结果显示没有
#还有其它多种方法

##3、cell annotation-----
# 对肿瘤细胞来说，分群后的细胞亚群注释是不可行的
# 这里仅仅是演示 SingleR 做 cell annotation的流程
library(SingleR)
refdata <- get(load("../../rawdata/HumanPrimaryCellAtlasData.Rdata"))
#参考数据库，等待时间较长。建议下载成功后，储存为Rdata，以后方便使用。
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    # label.finea耗时比较长一点
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
rm(refdata, HumanPrimaryCellAtlasData, testdata) #珍惜内存
table(cellpred$labels)
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
table(celltype$ClusterID,celltype$celltype) #如下为singleR的细胞cluster鉴定结果。
#结合上述结果，给scRNA增添celltype注释信息
scRNA@meta.data$celltype = "NA"
#先新增列celltype，值均为NA，然后利用下一行代码循环填充
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 <- DimPlot(scRNA, group.by="celltype", label=F , reduction='tsne')
p1
ggsave("../../out/3.3celltype_anno.pdf", plot = p1, width = 18, height = 12)



##4、轨迹分析----
#  BiocManager::install('monocle')
library(monocle)
data <- as(as.matrix(scRNA@assays$RNA@counts), 'sparseMatrix')
# count矩阵
pd <- new('AnnotatedDataFrame', data = scRNA@meta.data)
# meta表转成特定格式
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# 基因名表转成特定格式
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
#expressionFamily参数用于指定表达矩阵的数据类型，有几个选项可以选择：
#稀疏矩阵用negbinomial.size()，FPKM值用tobit()，logFPKM值用gaussianff()
save(mycds, file = "../../tmp/mycds_raw.Rdata")
rm(list = ls())
load("../../tmp/mycds_raw.Rdata")
#library("monocle")
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
#完成数据导入和预处理后，就可以考虑选择特定基因代表细胞的发育特征
#这里可以选取我们之前挑选的marker gene
load("../../tmp/markergene.Rdata")
markers.gene <- all.markers$gene
mycds <- setOrderingFilter(mycds, markers.gene)

#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#耗时，耗内存
#排序
mycds <- orderCells(mycds)
save(mycds,file = "../../tmp/mycds_reduced.Rdata")

p1 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("../../out/3.4trajectory_1.pdf", plot = p1) 
p2 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("../../out/3.4trajectory_2.pdf", plot = p2) 
rm(list = ls())
