setwd("D:/myprojects/20201109-DEG")
library("survival")
library("ggplot2")
library("survminer")
library(pheatmap)

#OS_data = read.table("20201109-OS-sum.txt",header = T,row.names = 1,sep = '\t')
#fit <- survfit(Surv(times,Status) ~ EZH2, data = OS_data)
#summary(fit)$table
#ggsurvplot(fit, 
#           pval = TRUE, int = TRUE,
#           palette = c("#FF0000","#0000FF"))



aaa = read.table("20201109-CIBERSORT-CF.txt",header = T,row.names = 1,sep = '\t')
aaa = t(aaa)
bbb = read.table("CF-Group.txt",header = T,row.names = 1,sep = '\t')
bbb$Mast.cells.resting = as.factor(as.character(bbb$Mast.cells.resting))
pheatmap(aaa,
         show_colnames = F ,
         annotation_col = bbb,
         #cluster_rows = F,
         cluster_cols = F,
         #color = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
         #cluster_cols = T,
         #treeheight_row = 0,
         #border_color = 'gray',
         #show_rownames=T,
         #angle_col = 45,
         #fontsize_row = 8,fontsize_col = 8,
         #filename = "test.pdf"
)






dev.off()
