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

dim(counts)
library(destiny)
# DiffusionMap class {destiny}
destinyObj <- as.ExpressionSet(as.data.frame(t(counts)))
destinyObj$condition <- factor(phe$celltype)
sigma=50
# 类似于 PCA 分析
dm <- DiffusionMap(destinyObj, sigma, rotate = TRUE)
# 也可以看主成分 
plot(
  eigenvalues(dm), 
  ylim = 0:1, 
  pch = 20, 
  xlab ='Diffusion component (DC)', 
  ylab ='Eigenvalue'
) 

# 使用 Slingshot 来进行 pseudotime 分析，基于 diffusion map的主要成分
# 前面的图看的 DC3到DC4之间有一个拐点，所以我们考虑前面的4个diffusion map的主要成分 
# https://bioconductor.org/packages/release/bioc/html/slingshot.html 

table(dm$condition)
library(slingshot)
crv <- slingshot(
  dm@eigenvectors[,1:15], 
  dm$condition, 
  start.clus = 'FCGR3A+ Mono', 
  end.clus='CD14+ Mono',
  maxit=100000,
  #   shrink.method="cosine" ,
shrink.method="tricube"
) 
crv
sce=crv
sce
# 下面的代码要求 slingshot_2.0.0 以上版本
# https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html
summary(sce$slingPseudotime_1)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.000   8.631  21.121  21.414  34.363  43.185
library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

