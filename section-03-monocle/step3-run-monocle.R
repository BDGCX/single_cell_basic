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
library(monocle)

load('input_cds.Rdata')
### 然后查看monocle ### 
cds 
# 接下来很重要，到底是看哪个性状的轨迹
table(pData(cds)$Cluster)
table(pData(cds)$Cluster,pData(cds)$celltype)
plot_cell_clusters(cds, 1, 2 )

## 我们这里并不能使用 monocle的分群
# 还是依据前面的 seurat分群, 也就是说前面的代码仅仅是流程而已，我们没有使用那些结果哦

# 其实取决于自己真实的生物学意图
pData(cds)$Cluster=pData(cds)$celltype
table(pData(cds)$Cluster)

Sys.time()
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Cluster")
Sys.time()

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name)) 
#  挑选差异最显著的基因可视化
plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
cg2=as.character(tail(sig_genes$gene_short_name)) 
plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )

# 前面是找差异基因，后面是做拟时序分析

# 第一步: 挑选合适的基因. 有多个方法，例如提供已知的基因集，
# 这里选取统计学显著的差异基因列表
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds) 
# 第二步: 降维。降维的目的是为了更好的展示数据。函数里提供了很多种方法,
# 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
# 第三步: 对细胞进行排序
cds <- orderCells(cds)
# 最后两个可视化函数，对结果进行可视化
plot_cell_trajectory(cds, color_by = "Cluster")  
ggsave('monocle_cell_trajectory_for_seurat.pdf')

length(cg)
plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Cluster") 
ggsave('monocle_plot_genes_in_pseudotime_for_seurat.pdf')

phe=pData(cds)
boxplot(phe$Pseudotime,phe$Cluster)

# https://davetang.org/muse/2017/10/01/getting-started-monocle/
# 前面根据差异基因，推断好了拟时序，也就是说把差异基因动态化了

# 后面就可以具体推断哪些基因随着拟时序如何的变化
my_cds_subset=cds
# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(my_cds_subset))
# 这个differentialGeneTest会比较耗费时间
my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 1 )
# 不知道为什么在Mac电脑无法开启并行计算了 ，不过我测试了在Windows 电脑设置cores = 4是可以的
# 如果你是Mac电脑，自己修改 cores = 1 即可 
head(my_pseudotime_de)
save( my_cds_subset,my_pseudotime_de,
      file = 'output_of_monocle.Rdata')




