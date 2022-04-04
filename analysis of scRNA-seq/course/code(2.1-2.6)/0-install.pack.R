rm(list = ls()) 
#清空当前工作空间变量  
options()$repos  
#查看当前工作空间默认的下载包路径
options()$BioC_mirror 
#查看使用BioCManager下载包的默认路径
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
# 指定使用BioCManager下载的路径
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
# 指定使用install.packages下载包的路径
options()$repos 
options()$BioC_mirror

getOption("BioC_mirror")
getOption("CRAN")
#CRAN基础包
options(CRAN="https://mirrors.ustc.edu.cn/CRAN/")
cran_packages <- c('tidyverse',
                   'ggplot2'
                   ) 
for (pkg in cran_packages){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)

#Bio分析包
Biocductor_packages <- c("Seurat",
						"scran",
						"scater",
						"monocle",
						"DropletUtils",
						"SingleR"
                         )
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
#options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# use BiocManager to install
for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

#最后检查下成功与否
for (pkg in c(Biocductor_packages,cran_packages)){
  require(pkg,character.only=T) 
}




### GEO 
# 
# GEO Platform (GPL)
# GEO Sample (GSM)
# GEO Series (GSE)
# GEO Dataset (GDS)

#一篇文章可以有一个或者多个GSE数据集，一个GSE里面可以有一个或者多个GSM样本。
#多个研究的GSM样本可以根据研究目的整合为一个GDS，不过GDS本身用的很少。
#而每个数据集都有着自己对应的芯片平台，就是GPL。(芯片名与基因名ID转换)
#https://blog.csdn.net/weixin_43569478/article/details/108079337
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=
###以下包用于芯片数据更合适，单细胞于GEO网页下载数据更方便
library(GEOquery)
#BiocManager::install("GEOquery")
library(GEOquery)
gse1009 <- getGEO('GSE1009', destdir=".")
class(gse1009)
length(gse1009)
a <- gse1009[[1]]
class(gse1009[1])
a
b <- exprs(a)
c <- pData(a)
a$platform_id
