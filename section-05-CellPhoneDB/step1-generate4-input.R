rm(list = ls())
library(Seurat)
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = '../section-01-cluster/basic.sce.pbmc.Rdata')

DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
ggsave('DimPlot-umap.pdf')
sce=pbmc 
table( Idents(sce ))
table(sce@meta.data$seurat_clusters) 
table(sce@meta.data$orig.ident) 
table(Idents(sce))

sce=subset(sce, downsample = 50)
sce
save(sce,file = 'downsample-sce.Rdata')

ct=sce@assays$RNA@counts
ct[1:4,1:4] 
ct=as.data.frame(ct)
ct=ct[rowSums(ct>1)>10,]
dim(ct)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
ids=AnnotationDbi::select(org.Hs.eg.db,keys = rownames(ct),
                          keytype= 'SYMBOL',columns = c('SYMBOL','ENSEMBL'))
head(ids)
ids= na.omit(ids)
ids=ids[!duplicated(ids$ENSEMBL),]
test_counts=ct[ids$SYMBOL,]
rownames(test_counts)=ids$ENSEMBL
sample_ann <- sce@meta.data
colnames(sample_ann)

test_meta=data.frame(Cell=rownames(sample_ann),
                     cell_type= Idents(sce))
test_counts=as.data.frame(test_counts)
identical(colnames(test_counts),test_meta$Cell)

test_counts=cbind(rownames(test_counts),test_counts)
colnames(test_counts)[1]='Gene'

write.table(test_counts, "test_counts.txt",  row.names=F, sep='\t',quote = F)
write.table(test_meta, "test_meta.txt", row.names=F, sep='\t',quote = F)





