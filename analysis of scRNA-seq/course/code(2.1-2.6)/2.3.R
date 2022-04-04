### 3、挑选hvg（高变）基因，可视化----
#highly Variable gene:简单理解sd大的
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 1500) 
#根据文献原图，挑选变化最大的1500个hvg
top10 <- head(VariableFeatures(scRNA), 10) 
top10
plot1 <- VariableFeaturePlot(scRNA) 
#标记top10 hvg
p1_3 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) +
  theme(legend.position = c(0.1,0.8)) +
  labs(tag = "C")
p1_3

#看看ERCC
ERCC <- rownames(scRNA)[grep("^ERCC-",rownames(scRNA))]
LabelPoints(plot = plot1, points = ERCC, repel = TRUE, 
            size=2.5,colour = "blue") +
  theme(legend.position = c(0.1,0.8)) +
  labs(tag = "C")
#可以直观看到ERCC均不是高变基因，而且部分的ERCC基因表达量确实很高
p1_1 | p1_2 | p1_3 #上图

# 这里展开介绍一下 scater 
# https://bioconductor.org/packages/release/bioc/html/scater.html
# Single-cell analysis toolkit工具箱 for expression 
#scater contains tools to help with the analysis of single-cell transcriptomic data, 
#focusing on low-level steps such as quality control, normalization and visualization.
#based on the SingleCellExperiment class (from the SingleCellExperiment package)
#关于sce对象，https://www.jianshu.com/p/9bba0214844b
library(scater)
ct=as.data.frame(scRNA@assays$RNA@counts)
pheno_data=scRNA@meta.data
sce <- SingleCellExperiment(
  assays = list(counts = ct), 
  colData = pheno_data
)
#SingleCellExperiment是SingleCellExperiment包的函数；在加载scater包时会一起加载
sce
?stand_exprs
stand_exprs(sce) <- log2(
  calculateCPM(sce) + 1)  #只对自己的文库的标准化
#assays(sce)
#sum(counts(sce)[,1])
#head(counts(sce)[,1])
#log2(1*10^6/422507+1)
#logcounts(sce)[1:4,1:4]
#exprs(sce)[1:4,1:4]
stand_exprs(sce)[1:4,1:4]

sce <- logNormCounts(sce) #可以考虑不同细胞的文库差异的标准化
assays(sce)
logcounts(sce)[1:4,1:4]

#https://osca.bioconductor.org/normalization.html#spike-norm
#基于ERCC的标准化方式也有许多优势（不同cell的量理论上是一样的），详见链接
#关于一些常见的FPKM等方式在番外篇会有简单的介绍与学习



#观察上面确定的高变top10基因在四个样本的分布比较
plotExpression(sce, top10 ,
               x = "Patient_ID",  colour_by = "Patient_ID", 
               exprs_values = "logcounts") 
# 下面的绘图非常耗时：(保存为本地文件查看比较高效，建议为pdf文件)
p1 <- plotHighestExprs(sce, exprs_values = "logcounts")###表达量最高的基因
ggsave("../../out/2.3HighestExprs.pdf", plot = p1, width = 15, height = 18) 
#如果按照ERCC 40%的过滤标准，ERCC表达量也十分大
?plotHighestExprs
# Sometimens few spike-in transcripts may also be present here, 
# though if all of the spike-ins are in the top 50, 
# it suggests that too much spike-in RNA was added

#后续步骤暂时还按照pctERCC=40的过滤标准的结果进行分析

#tips:后面主要还是基于Seurat对象，此处可以删除该变量，节约内存。
save(scRNA, file="../../tmp/2.3.Rdata")
rm(list=ls())
