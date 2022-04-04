### 1、下载、探索、整理数据----
## 1.1 下载、探索数据
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465
#点击下方supplementary files的GSE84465_GBM_All_data.csv.gz行中间的http下载
sessionInfo()
# 读取文件耗时比较长，请耐心等待
a <- read.table("../../rawdata/GSE84465_GBM_All_data.csv.gz")
a[1:4,1:4]
#行名为symbol ID
#列名为sample，看上去像是两个元素的组合。
summary(a[,1:4]) 
boxplot(a[,1:4])###可以观察到细胞的分布比较集中，只有个别异常值
head(rownames(a))
tail(rownames(a),10)
# 可以看到原文的counts矩阵来源于htseq这个计数软件，所以有一些不是基因的行需要剔除：
#  "no_feature"           "ambiguous"            "too_low_aQual"        "not_aligned"          "alignment_not_unique"
tail(a[,1:4],10)

a=a[1:(nrow(a)-5),]##去除最后5行不是基因的行

#原始counts数据
#点击下方SRA RUN Selector进去点击total那行的Metadata下载
#分析细胞对应的细胞/组织类型和样本来源
#3,589 cells of 4 human primary GBM samples, accession number GSE84465
#2,343 cells from tumor cores and 1,246 cells from peripheral regions
b <- read.table("../../rawdata/SraRunTable.txt",
                sep = ",", header = T)
b <- format(b,scientific=F)##数字不表示为科学计数
b[1:4,1:4]
table(b$Patient_ID) # 4 human primary GBM samples
table(b$TISSUE) # tumor cores and peripheral regions
table(b$TISSUE,b$Patient_ID)


## 1.2 整理数据 
# tumor and peripheral 分组信息
head(colnames(a))
head(b$plate_id)
head(b$Well)
#a矩阵行名（sample）并非为GSM编号，而主要是由相应的plate_id与Well组合而成

b.group <- b[,c("plate_id","Well","TISSUE","Patient_ID")]
paste0("X",b.group$plate_id[1],".",b.group$Well[1])
b.group$sample <- paste0("X",b.group$plate_id,".",b.group$Well)#为b添加一列sample
head(b.group)
identical(colnames(a),b.group$sample)###看a的列明与b.group的sample这一列是否相同

# 筛选tumor cell
index <- which(b.group$TISSUE=="Tumor")
length(index)
group <- b.group[index,] #筛选的是行
head(group)

a.filt <- a[,index] #筛选的是列
dim(a.filt)
identical(colnames(a.filt),group$sample)

sessionInfo()

