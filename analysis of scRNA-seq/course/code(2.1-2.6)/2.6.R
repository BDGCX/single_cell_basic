### 6、拼图，比较----
p <- (p1_1 | p1_2 | p1_3 ) /
  ((p2_1| p2_2 | p2_3) /
     (p3_1 | p3_2))
ggsave("../../out/my_try.pdf", plot = p, width = 15, height = 18) 

rm(list=ls())
save(scRNA,file = "scRNA.Rdata")