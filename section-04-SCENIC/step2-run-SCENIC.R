rm(list = ls()) 
library(Seurat) 

load(file = 'sce-monocyte.Rdata')
sce 
table(Idents(sce))
phe=sce@meta.data   
mat=sce@assays$RNA@counts

mat[1:4,1:4]
exprMat =as.matrix(mat) 
dim(exprMat)
exprMat[1:4,1:4] 
head(phe)

cellInfo <-  phe[,c('seurat_clusters','nCount_RNA' ,'nFeature_RNA' )]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)
cellInfo$CellType=Idents(sce)
table(cellInfo$CellType)

### Initialize settings
# https://github.com/aertslab/SCENIC
# https://pyscenic.readthedocs.io/en/latest/

library(SCENIC)
# https://resources.aertslab.org/cistarget/

db='~/Downloads/cisTarget_databases'
list.files(db)
# 保证cisTarget_databases 文件夹下面有下载好 的文件
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir=db , nCores=4) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
saveRDS(cellInfo, file="int/cellInfo.Rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
length(genesKept)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
# 最耗费时间的就是这个步骤
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
# 这个步骤也很耗时
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
# 因为莫名其妙的错误，需要把 多线程重新设置成为 1 个线程
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir=db , nCores=1) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
# 
# Binary regulon activity: 29 TF regulons x 93 cells.
# (34 regulons including 'extended' versions)
# 29 regulons are active in more than 1% (0.93) cells.

tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# 运行 







