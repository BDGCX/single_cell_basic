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
load(file = 'sce-monocyte.Rdata')
sce
# 因为 monocle包是 使用 CellDataSet 对象，所以重新构建
library(monocle)
# scater 
sample_ann <-  sce@meta.data  
sample_ann$celltype=Idents(sce)
head(sample_ann)
# rownames(sample_ann)=sample_ann[,1]
gene_ann <- data.frame(
  gene_short_name = rownames(sce@assays$RNA) , 
  row.names =  rownames(sce@assays$RNA) 
)
head(gene_ann)

pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
ct=as.data.frame(sce@assays$RNA@counts)
ct[1:4,1:4]

sc_cds <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds

# 这个 CellDataSet 对象一定要认识清楚，务必花两个小时去摸索它。


# 接下来仅仅是  monocle的标准流程而已
library(monocle)
sc_cds
sc_cds <- detectGenes(sc_cds, min_expr = 1) 
# 数值可以自行摸索
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]
sc_cds

cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 
# plyr 

# 并不是所有的基因都有作用，所以先进行挑选，合适的基因用来进行聚类。
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table,
                                 mean_expression >= 0.1)
unsup_clustering_genes
cds <- setOrderingFilter(cds, 
                         unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, 
                           return_all = F) # norm_method='log'
# 其中 num_dim 参数选择基于上面的PCA图
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 6) 
plot_cell_clusters(cds, 1, 2 )
table(pData(cds)$Cluster) 
colnames(pData(cds)) 
table(pData(cds)$seurat_clusters)
table(pData(cds)$Cluster,pData(cds)$seurat_clusters)
table(pData(cds)$Cluster,pData(cds)$celltype)
plot_cell_clusters(cds, 1, 2 )

# 可以看到 monocle 给细胞重新定义了亚群，亚群数量是自己选择的
# 整体来说，monocle和seurat 各自独立流程定义的亚群的一致性还不错。

# 只是跑流程而已
save(cds,file = 'input_cds.Rdata')

# 构建对象，seurat，monocle，scater
# monocle标准流程，降维聚类分群 



