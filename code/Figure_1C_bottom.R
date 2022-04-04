rm(list = ls())
library(ggplot2)
library(survival)
library(GGally)

inputClinical <- '../data/Figure 1C bottom input.txt'
data <- read.delim(inputClinical, sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
head(data)
# 临床资料，很容易走一下简单的 KM 生存分析
model <- Surv(data$`DFS (month)`, data$`Recurrence status`) ~ data$`Operational method`
data.surv <- survfit(model, data)
data.diff <- survdiff(model)
pval <- pchisq(data.diff$chisq, length(data.diff$n)-1, lower.tail = FALSE)
pval <- format(signif(pval, digits = 3), scientific = TRUE)
pval <- paste(c('P <', pval), collapse = ' ')

plt <- ggsurv(data.surv, surv.col = c('#ee756d', '#1aa6b8'), plot.cens = FALSE) + 
  labs(title = 'Operational methods', x = 'Months', y = 'Disease-free survival') + 
  guides(linetype = F) + 
  theme_classic(base_size = 7) + 
  theme(axis.text = element_text(colour="black"), 
        axis.ticks = element_line(colour='black'),
        plot.title = element_text(hjust=0.5, colour='black'), 
        legend.title = element_blank(), 
        legend.position=c(0.7,0.9),
        legend.background = element_blank()) + 
  annotate("text", x = 100, y = 0.1, label = pval, size = 2) + 
  ylim(0,1)
plt

ggsave('../results/Figure 1C bottom.pdf', plt, width = 6, height = 6, units = 'cm')

## 加载R包
library(survival)
library(survminer)
library(tidyverse)

head(data)
your_data=data
colnames(your_data)=c('stage','status','OS')
##一、画生存曲线
fit <- survfit(Surv(OS, status) ~ stage, data = your_data)

# log-rank test：pvalue
# This function implements the G-rho family of Harrington and Fleming (1982), with weights on each death of S(t)^rho, where S is the Kaplan-Meier estimate of survival. 
# With rho = 0 this is the log-rank or Mantel-Haenszel test, and with rho = 1 it is equivalent to the Peto & Peto modification of the Gehan-Wilcoxon test.
sdf <- survdiff(Surv(your_data$OS,your_data$status)~your_data$stage,rho=0)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
p.val

## photo 1
ggsurvplot(fit,
           pval = TRUE, #添加log-rank检验的p值
           pval.method = TRUE,#添加p值的检验方法
           conf.int = TRUE,#添加置信区间
           risk.table = TRUE, #在下方添加风险表
           risk.table.col = "strata", #根据数据分组为风险表添加颜色
           linetype = "strata", #改变不同组别的生存曲线的线条类型
           surv.median.line = "hv", #标注出中位生存时间
           xlab = "Time in years", #x轴标题
           xlim = c(0,5), #展示x轴的范围
           break.time.by = 2, #x轴间隔
           size=1, #线条大小
           ggtheme = theme_bw(), #为图形添加网格
           palette = c("#00BFFF","#DAA520")) #图形颜色风格

## photo 2
photo2 = ggsurvplot(fit,
                    legend.title = "Stage",#定义图例的名称
                    #legend.labs = c("Stage I/II","Stage III/IV"),
                    #legend = "top",#图例位置
                   # pval = 5.905137e-09, #在图上添加log rank检验的p值
                    #pval.method = TRUE,#添加p值的检验方法
                    #conf.int = TRUE,#添加置信区间
                    risk.table = TRUE, #在图下方添加风险表
                    #risk.table.col = "strata", #根据数据分组为风险表添加颜色
                    risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
                    #linetype = "strata", #改变不同组别的生存曲线的线型
                    #surv.median.line = "hv", #标注出中位生存时间
                    xlab = "Time in years", #x轴标题
                    xlim = c(0,10), #展示x轴的范围
                    break.time.by = 2, #x轴间隔
                    size = 2, #线条大小
                    #ggtheme = theme_bw(), #为图形添加网格
                    palette = c("#00BFFF","#DAA520")#图形颜色风格
) 


photo2  #看一下图

# 修改图例
# 修改风险表的图例名称 
photo2$table <- photo2$table + labs(
  title = "Number at risk")
# Changing the font size, style and color of photo2
# survival curves, risk table
photo2 <- ggpar(
  photo2,
  font.title= c(16, "bold", "darkblue"),#16号字体，粗体，darkblue色         
  font.x = c(14, "bold.italic", "red"),#14号字体，粗斜体，red          
  font.y = c(14, "bold.italic", "darkred"), #14号字体，粗斜体，darkred     
  font.xtickslab = c(12, "plain", "darkgreen"),#调整X轴上字体大小、颜色和风格
  legend = "top" #图例的位置
)

photo2  #再看一下图


## photo 3 ：还可以添加删失表，不过一般不加
photo3 = ggsurvplot(fit,
                    legend.title = "Stage",#定义图例的名称
                    #legend.labs = c("Stage I/II","Stage III/IV"),
                    ncensor.plot = TRUE,# 增加删失表
                    ncensor.plot.height = 0.25,
                    #legend = "top",#图例位置
                    pval = 5.905137e-09, #在图上添加log rank检验的p值
                    #pval.method = TRUE,#添加p值的检验方法
                    #conf.int = TRUE,#添加置信区间
                    risk.table = TRUE, #在图下方添加风险表
                    #risk.table.col = "strata", #根据数据分组为风险表添加颜色
                    risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
                    #linetype = "strata", #改变不同组别的生存曲线的线型
                    #surv.median.line = "hv", #标注出中位生存时间
                    xlab = "Time in years", #x轴标题
                    xlim = c(0,125), #展示x轴的范围
                    break.time.by = 10, #x轴间隔
                    size = 2, #线条大小
                    #ggtheme = theme_bw(), #为图形添加网格
                    palette = c("#00BFFF","#DAA520")#图形颜色风格
) 

photo3

# 参考资料：http://www.sthda.com/english/rpkgs/survminer/index.html



