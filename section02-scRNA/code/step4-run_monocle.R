rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
## 首先载入文章的数据
load(file='../input.Rdata')
counts=a
# using raw counts is the easiest way to process data through Seurat.
counts[1:4,1:4];dim(counts)
library(stringr) 
meta=metadata
head(meta) 

fivenum(apply(counts,1,function(x) sum(x>0) ))
boxplot(apply(counts,1,function(x) sum(x>0) ))
fivenum(apply(counts,2,function(x) sum(x>0) ))
hist(apply(counts,2,function(x) sum(x>0) ))
# 上面检测了 counts 和 meta 两个变量，后面需要使用
library(monocle)
## 首先创建对象
#monocle需要的用来构建 CellDataSet 对象的三个数据集
#1.表达量矩阵exprs:数值矩阵 行名是基因, 列名是细胞编号.
#2.细胞的表型信息phenoData: 第一列是细胞编号，其他列是细胞的相关信息
#3.基因注释featureData: 第一列是基因编号, 其他列是基因对应的信息
#并且这三个数据集要满足如下要求:
#表达量矩阵必须：
#保证它的列数等于phenoData的行数
#保证它的行数等于featureData的行数
#而且phenoData的行名需要和表达矩阵的列名匹配
#featureData和表达矩阵的行名要匹配
#featureData至少要有一列“gene_short_name”, 就是基因的symbol
gene_ann <- data.frame(
  gene_short_name = row.names(counts), 
  row.names = row.names(counts)
)
sample_ann=meta
pd<- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
sc_cds <- newCellDataSet(
  as.matrix(counts), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds

library(dplyr)
colnames(phenoData(sc_cds)@data)
## 必要的归一化 
sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- estimateDispersions(sc_cds)
## 然后进行一定程度的过滤
#ERCC基因很多的话需要去除ERCC基因
## 接下来的分析，都是基于sc_cds对象

cds=sc_cds
cds
##基因层面的过滤
## 起初是 24582 features, 768 samples 
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 5))#在大于5个细胞中表达的基因
length(expressed_genes)
cds <- cds[expressed_genes,]
cds
# 过滤基因后是  14442 features, 768 samples 
print(head(pData(cds)))
tmp=pData(cds)
head(tmp)
# 这里简单过滤细胞，如果想做其它过滤，就自己摸索阈值，然后挑选细胞即可。
valid_cells <- row.names(tmp[tmp$num_genes_expressed>2000,] )#有大于2000个基因表达的细胞
cds <- cds[,valid_cells]
cds 

## 最后是 14442 features, 693 samples

## 挑选合适的基因进入下游分析，比如PCA或者tSNE

disp_table <- dispersionTable(cds)
#要理解disp_table矩阵里的4个指标然后根据自身要求挑选合适的基因进行分析
unsup_clustering_genes <- subset(disp_table, 
                                 mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
cds
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'

cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 5) 


plot_cell_clusters(cds, 1, 2, color = "g")
print(head(pData(cds)))
plot_cell_clusters(cds, 1, 2, color = "Cluster")

table(pData(cds)$Cluster,pData(cds)$g)
boxplot(pData(cds)$num_genes_expressed~pData(cds)$Cluster)
boxplot(pData(cds)$num_genes_expressed~pData(cds)$g)
## 去除检测到基因数量效应，可以不用去除
cds <- reduceDimension(cds, max_components = 2, num_dim = 2,
                       reduction_method = 'tSNE',
                       residualModelFormulaStr = "~num_genes_expressed",
                       verbose = T)
cds <- clusterCells(cds, num_clusters = 5)
plot_cell_clusters(cds, 1, 2, color = "Cluster")
table(pData(cds)$Cluster,pData(cds)$g)#比较层次聚类和monocle包细胞分组的不同
boxplot(pData(cds)$num_genes_expressed~pData(cds)$Cluster)

# pData(cds)$CellType=pData(cds)$Cluster

## 接着找差异基因

Sys.time()
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Cluster")
Sys.time()
# 可以看到运行耗时5分钟

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
dim(sig_genes)
head(sig_genes[,c("gene_short_name", "pval", "qval")] )




##  最后推断发育轨迹

## 首先挑选合适的基因
# 这里选取统计学显著的差异基因列表
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

# 然后降维
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
# 降维是为了更好的展示数据。
# 降维有很多种方法, 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法

# 接着对细胞进行排序
cds <- orderCells(cds)

## 最后两个可视化函数 
plot_cell_trajectory(cds, color_by = "Cluster")  
# 可以很明显看到细胞的发育轨迹

## 这里可以展现marker基因在发育轨迹推断的效果，本例子随便 选取了6个差异表达基因。
plot_genes_in_pseudotime(cds[head(sig_genes$gene_short_name),], 
                         color_by = "Cluster")











