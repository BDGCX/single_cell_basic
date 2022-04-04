rm(list = ls())  
library(SCENIC)

library(SCopeLoomR)
scenicLoomPath='output/scenic.loom'
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
exprMat_log[1:4,1:4] 

scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicOptions
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
library(tidyverse)
savedSelections <- shiny::runApp(aucellApp)




