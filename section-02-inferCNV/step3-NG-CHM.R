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

install.packages("devtools")
devtools::install_github("bmbroom/tsvio")
devtools::install_github("bmbroom/NGCHMR", ref="stable")
devtools::install_github("broadinstitute/infercnvNGCHM")
# https://www.ngchm.net/Downloads/index.html

library(NGCHM)
library(infercnvNGCHM)
run.final.infercnv_obj=readRDS('infercnv_output/run.final.infercnv_obj')
ngchm(infercnv_obj          =   run.final.infercnv_obj ,# Object created by inferCNV,
      out_dir              =   'infercnv_output/' ,# Directory to output the heat map file,
      path_to_shaidyMapGen =  'ShaidyMapGen.jar',#Pathway to the java application ShaidyMapGen.jar
      gene_symbol          =   "bio.gene.hugo" ,# Option to add linkouts
      )
