rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)
## 首先载入文章的数据
load(file='../input.Rdata')
counts=a
counts[1:4,1:4];dim(counts)
library(stringr) 
meta=metadata
head(meta) 

options(warn=-1) # turn off warning message globally
suppressMessages(library(scater))
## 创建 scater 要求的对象
#由于scater对单细胞转录组学数据进行可视化是基于SingleCellExperiment类，从而能够与其它Bioconductor包如scan、scuttled等共同操作，所以使用SingleCellExperiment()创建scater对象，因为counts必须是矩阵，所以此处使用as.matrix()将其转化为矩阵，colData是临床信息，assays和colData这两个参数必须设置
#sce主要结构组成,主要储存了4组相关信息

#（1）Assays，即counts表达矩阵的标准化处理的矩阵(可以有任意多种，但常见的也就两三种)；
#（2）colData，即scRNA-seq的每个细胞的信息(例如批次信息、分组信息、表达概况信息)；
#（3）rowData，即scRNA-seq的每个基因的信息(例如表达概况、不同类基因名ID)；
#（4）reducedDims，即每个细胞的降维特征信息(主要有PCA、tSNE、uMAP三类)

#相关函数命令
#assays(sce) 查看当前sce对象的所有assays'name[暂且可理解一种表达矩阵称之为一个assa]
#assay(sce,"name") 查看sce指定一种assay的表达矩阵针对常见的assay，例如count assay、logcounts assay。

#counts(sce)提取SingleCellExperiment对象中的count表达矩阵
sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(counts)), 
  colData = meta
)
sce

##对数据进行过滤和质控
## 只有运行了下面的函数后才有各式各样的过滤指标
genes=rownames(rowData(sce))
genes[grepl('^MT-',genes)]
genes[grepl('^ERCC-',genes)]
sce <- addPerCellQC(sce, 
                        subsets = list(ERCC = grep('^ERCC',genes)))

##一些可视化
#sum：为每一个细胞得到的所有基因的表达量之和，即测序总量
#detected：高于阈值的观察数，即每一个细胞中基因表达量高于阈值（前面设置的threshold）的基因的数量。即测序总量中合格的序列
#图中每个点代表1个细胞，我们希望看到随着总的counts增加，探测到的基因数也不断增加。
plotColData(sce, x = "sum", y="detected", colour_by="g")

#如果一个基因在总基因表达量上的比例多，在ERCC中的比例少，就是正常细胞，反之是测序质量低的细胞
plotColData(sce, x = "sum", y="subsets_ERCC_percent", 
            other_fields="g") + facet_wrap(~g)

#plotScater是先在表达量最高的基因中选一部分(默认是500)，然后从高到低累加，看她们对每个细胞文库的贡献如何，它将不同细胞的不同表达分布绘制出来，为了看不通风细胞表达量的差异，可以利用colData中的数据进行分类。
plotScater(sce,nfeatures = 300,block1 = "plate",colour_by = "g",exprs_values = "counts")

#在基因层面上，我们可以看到(默认为50)高表达的前50个基因，图里的每一行代表一个基因，每一个bar(图例的黑色竖线)代表一个基因在一个细胞里的表达，
#圆圈代表这个基因表达的中位数，我们希望看到一些寻常的对象，如线粒体基因，actin,核糖体基因，MALAT1等。
#一些spike-in转录本(ERCC基因)也可能出现在top50里，这表明加入太多的spike-inRNA，大量的假基因或预测基因表明序列(alignment)比对有问题。
plotHighestExprs(sce, exprs_values = "counts")
#从上图看出top50里ERCC基因很多，可以去除ERCC基因
##基因层面的过滤
# 根据基因的表达过滤掉那些在所有细胞中表达量之和大于5的基因，该参数自主确定
keep_feature <- rowSums(counts(sce) > 0) > 5
table(keep_feature)
sce <- sce[keep_feature,]
#细胞层面的过滤
##视频教程
tf=sce$n_g
#n_g为section01-RNA-seq中的第一步处理数据计算得到的
#n_g = apply(a,2,function(x) sum(x>1)) #统计每个样本有表达的有多少行（基因）
boxplot(tf)
fivenum(tf)
table(tf>2000)
sce=sce[,tf > 2000 ]
sce
head(meta)
##细胞筛选官网函数有待进一步探索
#1.isOutlier函数提供了一种更为合适的方法进行细胞筛选。它定义距离中位数一定数量的中位数绝对偏差(MADs)的阈值，超过这个值的细胞被认为是离群的，从而被过滤掉
keep_total <- isOutlier(sce$sum,nmads = 3,type = "lower",log = TRUE)
filtered <- sce[,!keep_total]
#2.
per.cell <- perCellQCMetrics(sce, 
                         subsets = list(ERCC = grep('^ERCC',genes)))
qc.stats2 <- perCellQCFilters(per.cell, 
                              sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent")) 
colSums(as.matrix(qc.stats2))
filtered <- sce[,!qc.stats2$discard]
#3.quickpercellQC函数还可以更方便的检测几个常见指标异常值
#该方法使用总counts数、检测到的基因数和某个基因集的count百分比(如线粒体基因，spike-in转录本)来确定要舍弃哪些细胞
is.ercc <- grep("^ERCC-", rownames(sce))
filtered<- quickPerCellQC(sce, subsets=list(ERCC=is.ercc), sub.fields="subsets_ERCC_percent")

##对数据进行标准化
sce <- logNormCounts(sce)
assayNames(sce)
summary(librarySizeFactors(sce))
plotExplanatoryVariables(sce)
##变量水平的质控，这里没有啥可做的变量所以不做
vars <- getVarianceExplained(sce, 
                             variables=c("tissue", "age"))
head(vars)
plotExplanatoryVariables(vars)

##基因表达可视化理论上应该是跟384孔板 这个变量无关
plotExpression(sce, rownames(sce)[1:6])
plotExpression(sce, rownames(sce)[1:6], 
               x = "plate",
               exprs_values = "logcounts")
plotExpression(sce, rownames(sce)[1:6], 
               x = "g",exprs_values = "logcounts")
plotExpression(sce, rownames(sce)[1:6],
               x = rownames(sce)[10])
plotExpression(sce, rownames(example_sce)[1:6],
               x = "plate", colour_by="g")

##降维
#PCA
sce <- runPCA(sce)
str(reducedDim(sce, "PCA"))
plotPCA(sce, colour_by="g")
plotPCA(sce, colour_by="plate")
sce <- runPCA(sce, name="PCA2",
                      subset_row=rownames(sce)[1:1000],
                      ncomponents=25)
str(reducedDim(sce, "PCA2"))
## PCA分布图上面添加临床信息
plotReducedDim(sce, dimred = "PCA", colour_by = "g")
plotReducedDim(sce, dimred = "PCA", 
               shape_by= "plate", 
               colour_by= "g")
plotReducedDim(sce, dimred = "PCA", 
               shape_by= "plate", 
               colour_by= "g",
               size_by = "基因名字")#看特定基因在PCA中起的作用

sce <- runPCA(sce, ncomponents=20)
plotPCA(sce, ncomponents = 4, colour_by = "g")
#reducedDimNames(sce)

## 考虑 ERCC 影响后继续PCA
#sce2 <- runPCA(sce, 
               #feature_set = rowData(sce)$is_feature_control)

## PCA分布图上面添加临床信息--------------
plotReducedDim(sce2, dimred = "PCA", 
               shape_by= "plate", 
               colour_by= "g")

## 运行 tSNE 降维算法
set.seed(1000)
sce <- runTSNE(sce, perplexity=10)
head(reducedDim(sce, "TSNE"))
plotTSNE(sce, 
         shape_by= "plate", 
         colour_by= "g")
set.seed(1000)
sce <- runTSNE(sce, perplexity=50, 
              dimred="PCA", n_dimred=10)
head(reducedDim(example_sce, "TSNE"))
plotTSNE(sce, 
         shape_by= "plate", 
         colour_by= "g")
## 对tSNE降维后结果进行不同的聚类,一般写个循环，把centers = 2到10都画一下，选个最合适的
colData(sce)$tSNE_kmeans <- as.character(kmeans(sce@int_colData$reducedDims$TSNE,
                                                centers = 5)$clust)
head(sce@int_colData$reducedDims$TSNE)
hc=hclust(dist( sce@int_colData$reducedDims$TSNE ))
clus = cutree(hc, 5) 
colData(sce)$tSNE_hc <-  as.character(clus)
plotTSNE(sce,  colour_by = "tSNE_kmeans")
plotTSNE(sce,  colour_by = "tSNE_hc")
table(colData(sce)$tSNE_hc , colData(sce)$tSNE_kmeans)
for (i in 2:10){
      colData(sce)$tSNE_kmeans <- as.character(kmeans(sce@int_colData$reducedDims$TSNE,                                                               centers = i)$clust)
      head(sce@int_colData$reducedDims$TSNE)
      hc=hclust(dist( sce@int_colData$reducedDims$TSNE ))
      clus = cutree(hc, i) 
      colData(sce)$tSNE_hc <-  as.character(clus)
      p1=plotTSNE(sce,  colour_by = "tSNE_kmeans")
      p2=plotTSNE(sce,  colour_by = "tSNE_hc")
      print(p1,p2)
 }
## 同样是一直降维方式，不同的算法
library(destiny)
sce <- runDiffusionMap(sce)
plotDiffusionMap(sce,  
                 shape_by= "plate", 
                 colour_by= "g")

##UMAP
sce <- runUMAP(sce)
head(reducedDim(sce, "UMAP"))
plotUMAP(sce,colour_by="g")

sessionInfo()

library(SC3) # BiocManager::install('SC3')
sce <- sc3_estimate_k(sce)
metadata(sce)$sc3$k_estimation
rowData(sce)$feature_symbol=rownames(rowData(sce))
# 耗费时间
kn=4
sc3_cluster="sc3_4_clusters"
# 非常耗时
sce <- sc3(sce, ks = kn, biology = TRUE)

sc3_plot_consensus(sce, k = kn, show_pdata = c("g",sc3_cluster))
sc3_plot_expression(sce, k = kn, show_pdata =  c("g",sc3_cluster))
sc3_plot_markers(sce, k = kn, show_pdata =  c("g",sc3_cluster))
plotPCA(sce, shape_by= "g" , colour_by =  sc3_cluster )
sc3_interactive(sce)





