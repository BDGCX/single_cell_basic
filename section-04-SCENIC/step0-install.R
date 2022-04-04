# http://www.bio-info-trainee.com/3727.html
options()$repos 
options()$BioC_mirror
#options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror

BiocManager::install(c("GENIE3", "AUCell", "RcisTarget"))
devtools::install_github("aertslab/SCENIC")

# Downloading GitHub repo aertslab/SCENIC@HEAD
# Error in utils::download.file(url, path, method = method, quiet = quiet,  : 
#   无法打开URL'https://api.github.com/repos/aertslab/SCENIC/tarball/HEAD'
#                

## 教程列表：
# SCENIC转录因子分析结果的解读，https://mp.weixin.qq.com/s/eAfkhX0SJu1lytZeXdsh0Q
# 单细胞转录因子分析之SCENIC流程， https://mp.weixin.qq.com/s/pN4qWdUszuGqr8nOJstn8w 
# 基因集的转录因子富集分析，https://mp.weixin.qq.com/s/mBR3IwWvQDcTOXNwM_YCEg

# 单细胞｜pySCENIC转录因子分析：从表达矩阵到结果可视化, https://mp.weixin.qq.com/s/gTOKdqawzqBPNokc-aEomQ

# 将模块中的基因输入NetworkAnalyst在线工具中，预测可能调控hub基因的转录因子TFs，并使用Cytoscape完成调控网络的绘制





