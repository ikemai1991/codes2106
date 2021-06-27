setwd('D:/myprojects/20201106-TOP ana')
source('D:/myprojects/20201106-TOP ana/myvolcano.R')
source('D:/myprojects/20201106-TOP ana/myGO_KEGG.R')
library(ggplot2)
library(limma)

#读表预处理
library(ggplot2)
library(limma)
merge_count = read.table('FPKM_LIHC.txt',header = T,row.names = 1,sep = '\t')
exprSet <- merge_count
colnames(merge_count)
group_list = read.table('Read_Group.txt',header = F,row.names = 1)
group_list = group_list$V2
group_list = as.factor(as.character(group_list))
rownames(group_list)


# exprSet <- log2(exprSet+1)
design <- model.matrix(~0+factor(group_list))
# design <- model.matrix(~factor(group_list$group))
head(design)
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
# colnames(design)[2] = 'T_p'
head(design)


contrast.matrix<-makeContrasts(paste0(c('Low','High'),collapse = "-"),levels = design)
contrast.matrix



#limma分析
library(limma) 
fit <- lmFit(exprSet,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = "Low-High", number = Inf)
summary(tempOutput)
dim(tempOutput)

DEOutput <- topTable(fit2, coef= "Low-High",number=Inf,p.value= 1,adjust.method ="BH")
DEOutput = subset(DEOutput,subset=P.Value<0.05)
dim(DEOutput)

summary(DEOutput)
DEOutput


#write.csv(DEOutput,'DEG_limma.csv')
write.table(DEOutput,'20201107-DEG.txt',sep = '\t')


#------火山图
df = DEOutput
colnames(df)
colnames(df)[1] = 'avg_logFC'
colnames(df)[4] = 'p_val'
dim(df)
df = df[-which(abs(df$avg_logFC)>100),]

head(df)
###火山???
library(ggrepel)
library(dplyr)
myvolcano(df,gene.plot = 10,logfc.cutoff=0.5,p.cutoff=0.05)


#DEG热图
DEOutput.1 = subset(DEOutput,logFC>1)
tmps = exprSet[rownames(DEOutput.1),]
grouplist.2 = read.table('Read_Group.txt',header = F)
head(grouplist.2)
grouplist.2$V1 = gsub('-','\\.',grouplist.2$V1)

tmps = tmps[,subset(grouplist.2,V2=='Low' | V2=='High')$V1]

count_data = tmps
count_data = t(scale(t(count_data)))
dim(count_data)
library(pheatmap)
count_data[1:4,1:3]
aaa = subset(grouplist.2,V2=='Low' | V2=='High')
aaa$V2 = as.factor(as.character(aaa$V2))
rownames(aaa)  = aaa$V1
aaa = data.frame(row.names = rownames(aaa),group = aaa$V2)
#png("heatmaptest-2.png",width = 456, height = 757)
pheatmap(count_data,
         annotation_col = aaa,
         cluster_rows = TRUE,
         color = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
         cluster_cols = T,
         treeheight_row = 0,
         border_color = 'gray',
         show_rownames=T,
         angle_col = 45,
         fontsize_row = 8,fontsize_col = 8
)

dev.off()


#KEGG&GO
upgene = rownames(subset(DEOutput, P.Value < 0.1 & logFC > 0.5))
downgene = rownames(subset(DEOutput, P.Value < 0.1 & logFC<(-0.5)))
myGO_KEGG(species = 'h',up_gene = upgene,down_gene = downgene)

