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
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = '../section-01-cluster/basic.sce.pbmc.Rdata')

DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
sce=pbmc 
table( Idents(sce ))
table(sce@meta.data$seurat_clusters) 
table(sce@meta.data$orig.ident) 

dat=GetAssayData(sce,
             slot='counts',assay='RNA')
dat[1:4,1:4]
 
groupinfo=data.frame(v1=colnames(dat),
                     v2= Idents(sce ) )
head(groupinfo)

# gtf 文件 
library(AnnoProbe)
geneInfor=annoGene(rownames(dat),
                   "SYMBOL",'human')
colnames(geneInfor)
head(geneInfor)
geneInfor=geneInfor[with(geneInfor,
                         order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
## 这里可以去除性染色体
# 也可以把染色体排序方式改变
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
head(groupinfo)
dat[1:4,1:4]
table(groupinfo$v2)
head(groupinfo)

# 为了节约计算机资源，我直接抽样 即可
kp=sample(1:nrow(groupinfo),500)
groupinfo=groupinfo[kp,]
dat=dat[,kp]
# 如果是真实项目，而且你计算机资源是足够的
# 请忽略这个抽样的操作

expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
write.table(groupinfo,file = groupFiles,sep = '\t',
            quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',
            quote = F,col.names = F,row.names = F)

 
options(stringsAsFactors = F)
library(Seurat)
library(gplots)
library(ggplot2)

expFile='expFile.txt' 
groupFiles='groupFiles.txt'  
geneFile='geneFile.txt'

# 发现了样本名字导致的小bug，需要Linux和R语言共同进行修改
# head -n 1 expFile.txt |tr '\t' '\n' > barcodes.txt
bd=read.table('barcodes.txt')[,1]
head(bd)
groupinfo=read.table('groupFiles.txt',sep = '\t')
groupFiles='groupFiles.txt'
head(groupinfo)
groupinfo[,1]=bd
table(groupinfo[,2])
write.table(groupinfo,file = groupFiles,sep = '\t',
            quote = F,col.names = F,row.names = F)

 
library(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c("NK",'DC',
                                                      'Platelet')) 
## 这个取决于自己的分组信息里面的

# cutoff=1 works well for Smart-seq2, and 
# cutoff=0.1 works well for 10x Genomics
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir= "infercnv_output",  # dir is auto-created for storing outputs
                              cluster_by_groups=F,   # cluster
                              hclust_method="ward.D2", plot_steps=F)
 








