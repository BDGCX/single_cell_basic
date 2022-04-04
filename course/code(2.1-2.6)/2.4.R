load("../../tmp/2.3.Rdata")
### 4、降维，PCA分析，可视化----
#先进行归一化（正态分布）
scRNA <- ScaleData(scRNA, features = (rownames(scRNA)))
#储存到"scale.data"的slot里
GetAssayData(scRNA,slot="scale.data",assay="RNA")[1:8,1:4]
#对比下原来的count矩阵
GetAssayData(scRNA,slot="counts",assay="RNA")[1:8,1:4]
#scRNA@assays$RNA@
#PCA降维，利用之前挑选的hvg，可提高效率
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
#挑选第一，第二主成分对cell可视化
DimPlot(scRNA, reduction = "pca", group.by="Patient_ID")
#发现与原文献中颠倒了
?DimPlot
?RunPCA
#seed.use	:Set a random seed. By default, sets the seed to 42. 
#Setting NULL will not set a seed.
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA),seed.use=3)
#尝试了seed.use的不同取值发现图形只有四种变化（四个拐角），其中以seed.use=3为代表的一类与原文文献一致
DimPlot(scRNA, reduction = "pca", group.by="Patient_ID")
#与文献一致了。个人觉得颠倒与否如果只是随机种子的差别的话，对后续分析应该没影响
p2_1 <- DimPlot(scRNA, reduction = "pca", group.by="Patient_ID")+
  labs(tag = "D")
p2_1

DimPlot(scRNA, reduction = "pca",  split.by = 'Patient_ID')

#挑选主成分，RunPCA默认保留了前50个
scRNA <- JackStraw(scRNA,reduction = "pca", dims=20)
scRNA <- ScoreJackStraw(scRNA,dims = 1:20)

p2_2 <- JackStrawPlot(scRNA,dims = 1:20, reduction = "pca") +
  theme(legend.position="bottom") +
  labs(tag = "E")
p2_2
p2_3 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
p2_3
#结果显示可挑选前20个pc

p2_1| (p2_2 | p2_3) #中图

