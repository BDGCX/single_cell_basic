##3、cell annotation-----
# 对肿瘤细胞来说，分群后的细胞亚群注释是不可行的
# 这里仅仅是演示 SingleR 做 cell annotation的流程
library(SingleR)
??singleR
refdata <- get(load("../../rawdata/HumanPrimaryCellAtlasData.Rdata"))
assay(refdata)[1:4,1:4]
head(refdata@colData)
head(refdata)
ref <- HumanPrimaryCellAtlasData()
#参考数据库，等待时间较长。建议下载成功后，储存为Rdata，以后方便使用。
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    # label.fine耗时比较长一点
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
save(cellpred,file = "../../tmp/cellpred.Rdata")
load("../../tmp/cellpred.Rdata")
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

ggsave("../../out/3.3celltype_anno.pdf", plot = p1, width = 18, height = 12) 

