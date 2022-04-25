#monocle构建CDS需要3个矩阵：expr.matrix、pd、fd
# 将Seurat中的对象转换为monocle识别的对象
#cds <- importCDS(GetAssayData(sce))
#cds <- importCDS(sce)  不知道为啥，使用imporCDS总报错
#从Seurat对象中选择做拟时序的细胞群cluster
Mono_tj<-subset(sce, idents = c(0,1,2,3))

#使用经过处理后的Seurat对象中的数据构建monocle对象
#monocle构建CDS需要3个矩阵：表达矩阵expr.matrix、细胞表型信息pd、基因注释fd
Mono_matrix<-as(as.matrix(GetAssayData(Mono_tj,slot = "counts")), 'sparseMatrix')
#构建featuredata，一般featuredata需要两个col，一个是gene_id,一个是gene_short_name,row对应counts的rownames
feature_ann<-data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann)<-rownames(Mono_matrix)
#
Mono_fd<-new("AnnotatedDataFrame", data = feature_ann)
#
#Seurat object中的@meta.data一般会存放表型相关的信息如cluster、sample的来源、group等，所以选择将metadata转换为phenodata
sample_ann<-Mono_tj@meta.data
#rownames(sample_ann)<-colnames(Mono_matrix)

Mono_pd<-new("AnnotatedDataFrame", data =sample_ann)
#build new cell data set
Mono.cds<-newCellDataSet(Mono_matrix,phenoData =Mono_pd,featureData =Mono_fd,expressionFamily=negbinomial.size())

#查看phenodata、featuredata
head(pData(Mono.cds))
head(fData(Mono.cds))
#预处理
Mono.cds <- estimateSizeFactors(Mono.cds)
Mono.cds <- estimateDispersions(Mono.cds)
#根据作者的建议，即使在Seurat包中已经标准化处理过的数据，在转化到Monocle中时仍然需要再一次进行标准化。
#将表达矩阵中所有值进行log标准化
L <- log(exprs(Mono.cds[expressed_genes,]))

#将每个基因都标准化，melt方便作图
library(reshape)
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
#作图，查看标准化的基因表达值的分布
qplot(value, geom = "density", data = melted_dens_df) +
       stat_function(fun = dnorm, size = 0.5, color = 'red') +
       xlab("Standardized log(FPKM)") +
       ylab("Density")
#筛选基因,这里可以根据自己的需要筛选特定的基因
##1.使用monocle选择的高变基因（本次使用该方法）
disp_table <- dispersionTable(Mono.cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
Mono.cds <- setOrderingFilter(Mono.cds, unsup_clustering_genes$gene_id)
##2.使用seurat选择的高变基因
var.seurat <- VariableFeatures(sce)
Mono.cds <- setOrderingFilter(Mono.cds, var.seurat)
plot_ordering_genes(Mono.cds)
##3.使用clusters差异表达基因
deg.cluster <- FindAllMarkers(sce)
diff.genes <- subset(deg.cluster,p_val_adj<0.05)$gene
Mono.cds <- setOrderingFilter(Mono.cds, diff.genes)
plot_ordering_genes(Mono.cds)
#dpFeature根据伪时间表达pattern聚类基因
marker_genes <- row.names(subset(fData(Mono.cds),
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5",
                                                        "ANPEP", "PDGFRA","MYOG",
                                                        "TPM1",  "TPM2",  "MYH2",
                                                        "MYH3",  "NCAM1", "TNNT1",
                                                        "TNNT2", "TNNC1", "CDK1",
                                                        "CDK2",  "CCNB1", "CCNB2",
                                                        "CCND1", "CCNA1", "ID1")))
#自己选择一些marker基因
diff_test_res <- differentialGeneTest(Mono.cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(Mono.cds[sig_gene_names,],
                        num_clusters = 6,
                        cores = 1,
                        show_rownames = T)
#用DDRtree 进行降维分析
Mono.cds <- reduceDimension(
  Mono.cds,
  max_components = 2,
  method = 'DDRTree')
#计算psudotime值
Mono.cds <- orderCells(Mono.cds)
head(pData(Mono.cds))

plot_cell_trajectory(Mono.cds,cell_size = 1)
plot_cell_trajectory(Mono.cds, color_by = "seurat_clusters",cell_size = 1)
plot_cell_trajectory(Mono.cds, color_by = "Pseudotime")
plot_cell_trajectory(Mono.cds, color_by = "seurat_clusters") +
  facet_wrap(~seurat_clusters, nrow = 1)
