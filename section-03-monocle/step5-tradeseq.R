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
colnames(pData(cds))
plot_cell_trajectory(cds)

# https://bioconductor.org/packages/release/bioc/html/tradeSeq.html  

library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(monocle)

# https://bioconductor.org/packages/release/bioc/vignettes/tradeSeq/inst/doc/Monocle.html
# https://rdrr.io/github/statOmics/tradeSeq/src/R/extract_monocle_info.R
source('extract_monocle_info.R')
info <- extract_monocle_info(cds)
info
sce <- fitGAM(counts = Biobase::exprs(cds),
              cellWeights = info$cellWeights,
              pseudotime = info$pseudotime)
 
# https://www.bioconductor.org/packages/release/bioc/vignettes/tradeSeq/inst/doc/tradeSeq.html 
table(rowData(sce)$tradeSeq$converged)
assoRes <- associationTest(sce)
head(assoRes)

# Discovering progenitor marker genes
startRes <- startVsEndTest(sce) 
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]

counts= Biobase::exprs(cds)
plotSmoothers(sce, counts, gene = sigGeneStart)

# Comparing specific pseudotime values within a lineage
customRes <- startVsEndTest(sce, pseudotimeValues = c(0.1, 0.8))
endRes <- diffEndTest(sce)  
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
plotSmoothers(sce, counts, sigGene)
# Discovering genes with different expression patterns
patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat]) 
plotSmoothers(sce, counts, gene = rownames(patternRes)[oPat][4])


