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

rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(gplots)
library(ggplot2)
library(monocle)
library(dplyr)
load(file = 'output_of_monocle.Rdata')
cds=my_cds_subset
phe=pData(cds)
counts = Biobase::exprs(cds)

sce <- SingleCellExperiment(assays = List(counts = counts))
sce

# 官网的示例数据，让你熟悉后面的要求
library(slingshot, quietly = FALSE)
data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl 
dim(rd) # data representing cells in a reduced dimensional space
head(rd)
table(cl)
length(cl) # vector of cluster labels
# 其中RD就是降维后的两个主成分，可以是PCA,UMAP 等等，也可以是  # DiffusionMap class {destiny}

## ----genefilt-----------------------------------------------------------------
# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]
sce # class: SingleCellExperiment 

## ----norm---------------------------------------------------------------------
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)
 
pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)


## ----umap, cache=TRUE---------------------------------------------------------
library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

# 可以看到 pca和umap其实都是可以自己制作，并不一定要 是seurat这样的一站式R包

## ----add_RDs, cache=TRUE------------------------------------------------------
reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)


# 聚类也可以自己选择R 包 
## ----clustering_mclust--------------------------------------------------------
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
table(cl1)
#table(cl1,phe$celltype)
colData(sce)$GMM <- cl1
library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
cl2 <- kmeans(rd1, centers = 4)$cluster
table(cl2)
colData(sce)$kmeans <- cl2
plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1) 

table(cl1,cl2)

# 降维聚类分群细节，但是不同流程理论上是有可比性
# 真正的函数在这里 slingshot
sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA') 
summary(sce$slingPseudotime_1) 

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

## ----plot_curve_2-------------------------------------------------------------
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

## ----tradeseq, eval=FALSE-----------------------------------------------------
 library(tradeSeq) 
 sce <- fitGAM(sce) 
 # test for dynamic expression
 ATres <- associationTest(sce) 
 topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250];topgenes
 pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
 heatdata <- assays(sce)$counts[topgenes, pst.ord]
 heatclus <- sce$GMM[pst.ord]

 heatmap(log1p(heatdata), Colv = NA,
         ColSideColors = brewer.pal(9,"Set1")[heatclus]) 
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

rd=rd2
cl=cl2

lin1 <- getLineages(rd, cl, start.clus = '1')
lin1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')

## ----lines_sup_end------------------------------------------------------------
lin2 <- getLineages(rd, cl, start.clus= '1', end.clus = '3')

plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)

## ----curves-------------------------------------------------------------------
crv1 <- getCurves(lin1)
crv1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')

## ----sling_approxpoints-------------------------------------------------------
sce5 <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA',
                  approx_points = 5)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce5$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce5)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce5), lwd=2, col='black')

## ----sling_omega--------------------------------------------------------------
rd2 <- rbind(rd, cbind(rd[,2]-12, rd[,1]-6))
cl2 <- c(cl, cl + 10)
pto2 <- slingshot(rd2, cl2, omega = TRUE, start.clus = c(1,11))

plot(rd2, pch=16, asp = 1,
     col = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl2])
lines(SlingshotDataSet(pto2), type = 'l', lwd=2, col='black')

## ----sling_multtraj-----------------------------------------------------------
plot(rd2, pch=16, asp = 1,
     col = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl2])
lines(SlingshotDataSet(pto2), lwd=2, col='black')

## ----session------------------------------------------------------------------
sessionInfo()


