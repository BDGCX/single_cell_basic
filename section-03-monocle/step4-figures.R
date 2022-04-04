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

## 后面是对前面的结果进行精雕细琢
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
colnames(phe)
library(ggsci)
p1=plot_cell_trajectory(cds, color_by = "Cluster")  + scale_color_nejm() 
p1
ggsave('trajectory_by_cluster.pdf')
plot_cell_trajectory(cds, color_by = "celltype")  

p2=plot_cell_trajectory(cds, color_by = "Pseudotime")  
p2
ggsave('trajectory_by_Pseudotime.pdf')

p3=plot_cell_trajectory(cds, color_by = "State")  + scale_color_npg()
p3
ggsave('trajectory_by_State.pdf')
library(patchwork)
p1+p2/p3

phe=pData(cds)
head(phe)
table(phe$State,phe$Cluster) 

library(dplyr)
my_pseudotime_de %>% arrange(qval) %>% head() 
# save the top 6 genes
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene=my_pseudotime_gene[,1]
my_pseudotime_gene
plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])+ scale_color_npg()
ggsave('monocle_top6_pseudotime_by_state.pdf')

plot_genes_jitter(my_cds_subset[my_pseudotime_gene,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )+ scale_color_nejm()
ggsave('monocle_top6_pseudotime_by_cluster.pdf')


# cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster[,1]
gene_to_cluster
colnames(pData(my_cds_subset))
table(pData(my_cds_subset)$Cluster,pData(my_cds_subset)$State) 
ac=pData(my_cds_subset)[c('celltype','State','Pseudotime')]
head(ac)
# 这个热图绘制的并不是纯粹的细胞基因表达量矩阵，而是被 Pseudotime 好了的100列，50行的矩阵

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],
                                                 # num_clusters = 2, 
                                                 # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)
my_pseudotime_cluster

pdf('monocle_top50_heatmap.pdf')
print(my_pseudotime_cluster)
dev.off()

if(file.exists('BEAM_res.Rdata')){
  
  load(file = 'BEAM_res.Rdata')
}else{
  
  # 这个步骤超级耗费时间
  # 不知道为什么在Mac电脑无法开启并行计算了 
  # 不过我测试了在Windows 电脑设置cores = 4是可以的
  # 如果你是Mac电脑，自己修改 cores = 1 即可 
  BEAM_branch1 <- BEAM(my_cds_subset, branch_point = 1, cores = 4)
  BEAM_branch1 <- BEAM_branch1[order(BEAM_branch1$qval),]
  BEAM_branch1 <- BEAM_branch1[,c("gene_short_name", "pval", "qval")]
  head(BEAM_branch1) 
  
  BEAM_branch2 <- BEAM(my_cds_subset, branch_point = 2, cores = 4)
  BEAM_branch2 <- BEAM_branch2[order(BEAM_branch2$qval),]
  BEAM_branch2 <- BEAM_branch2[,c("gene_short_name", "pval", "qval")]
  head(BEAM_branch2)
   
  save(BEAM_branch1,BEAM_branch2,file = 'BEAM_res.Rdata')
  
}


# 使用全部的基因进行绘图 
BEAM_res = BEAM_branch1
my_branched_heatmap <- plot_genes_branched_heatmap(
  my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
  branch_point = 1,
  num_clusters = 4, 
  use_gene_short_name = TRUE,
  show_rownames = F,
  return_heatmap = TRUE)

pdf('monocle_BEAM_branch1_heatmap.pdf')
print(my_branched_heatmap$ph)
dev.off()

BEAM_res = BEAM_branch2
my_branched_heatmap <- plot_genes_branched_heatmap(
  my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
  branch_point = 1,
  num_clusters = 4, 
  use_gene_short_name = TRUE,
  show_rownames = F,
  return_heatmap = TRUE)

pdf('monocle_BEAM_branch2_heatmap.pdf')
print(my_branched_heatmap$ph)
dev.off()



head(my_branched_heatmap$annotation_row)
table(my_branched_heatmap$annotation_row$Cluster) 
my_row <- my_branched_heatmap$annotation_row
my_row <- data.frame(cluster = my_row$Cluster,
                     gene = row.names(my_row),
                     stringsAsFactors = FALSE)

head(my_row[my_row$cluster == 3,'gene']) 

my_gene <- row.names(subset(fData(my_cds_subset),
                            gene_short_name %in% head(my_row[my_row$cluster == 1,'gene'])))
my_gene
# plot genes that are expressed in a branch dependent manner
plot_genes_branched_pseudotime(my_cds_subset[my_gene,],
                               branch_point = 1,
                               ncol = 1)

plot_genes_branched_pseudotime(my_cds_subset[my_gene,],
                               branch_point = 2,
                               ncol = 1)
# 后面的批量绘图，意义不大 
names(pData(my_cds_subset))
head(pData(my_cds_subset))

plot_genes_jitter(my_cds_subset[gene_to_cluster,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow=  10,
                  ncol = NULL )

ggsave('monocle_top50_subCluster.pdf',height = 42)
plot_genes_in_pseudotime(my_cds_subset[head(gene_to_cluster,25),])
ggsave('monocle_top50_pseudotime.pdf',height = 49)


write.csv(my_pseudotime_de,file = 'my_pseudotime_de.csv')






