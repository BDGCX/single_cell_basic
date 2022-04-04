rm(list = ls())  
library(RColorBrewer)
library(ggplot2)
library(scales)

drPiechart <- function(columnNames, Values, Colors, outputPdf){
  data <- data.frame(
    group = columnNames,
    value = Values
  )
  data$group <- factor(data$group, columnNames)
  
  pie <- ggplot(data, aes(x="", y=value, fill=factor(group))) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0) + 
    scale_fill_manual(values = Colors) + 
    coord_polar(theta = "y", direction = -1) +
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank()) 
  
  ggsave(outputPdf, units = 'cm', height = 8, width = 16)
}

#drPiechart(c(), c(), colors = c(), 'output.pdf')
drPiechart(c('Yes', 'No', 'NA'), 
           c(21, 13, 1), 
           c('#e8137f', '#646262', '#bebbbb'),
           '../results/Figure 5C left.pdf')

 