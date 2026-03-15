library(DESeq2)
library(ggplot2)
library(gprofiler2)
library(dplyr)

countData <- read.csv('~/bulk_yo/data/Courties_G_03_OT_vs_YT_Filtered_Raw_Counts.txt', header = TRUE, sep = "\t")

row.names(countData)<-countData[,1]
countData<-countData[,2:ncol(countData)]
countData<-countData[,2:ncol(countData)]
countData[] <- lapply(countData, as.numeric)






age<-c(rep("old",3), rep("young",4))
metaData<-data.frame(age)


dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~age)

dds <- DESeq(dds)
vsdata <- vst(dds, blind=FALSE)



pcaData <- plotPCA(vsdata, intgroup="age", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
colors <- c("young" = "#173C70", "old" = "#6C6E38")

png(filename = "~/bulk_yo/results/macro/pca.png", width = 5, height = 2, units = "in", res=600)
ggplot(pcaData, aes(PC1, PC2, color=age)) +
  geom_point(size=3) +
  scale_color_manual(values=colors) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal( )+
  ggtitle("PCA Plot")+theme_bw()+
  expand_limits(x = c(-10,10), y=c(-10,10))
dev.off()




countn<-counts(dds, normalized=T)

res <- results(dds)
res$gene<-rownames(res)
res<-as.data.frame(res)
openxlsx::write.xlsx(res%>% filter(padj<=0.05), file = '~/bulk_yo/results/macro/degene.xlsx', rowNames=TRUE) 





tmp<- res 
tmp2<-strsplit(rownames(tmp), split = "\\|")
tmp3<-c()
for(el in tmp2){
  tmp3<-c(tmp3,el[3])
}

tmp$gene<-tmp3
saveRDS(tmp, "~/bulk_yo/results/macro/degenes.rds")

tmp$fcsign <- sign(tmp$log2FoldChange)
tmp$logP=-log10(tmp$pvalue)
tmp$metric= -tmp$logP/tmp$fcsign
y<-tmp[,c("gene","metric")]
filtered <- na.omit(y)

#gsea input
write.table(filtered ,file="~/bulk_yo/results/macro/rank_macro.rnk",quote=F,sep="\t",row.names=F, col.names = F)








library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel)







sc_mac.markers<-readRDS("~/young-oldv2/results/merge/merge-final/mac/markers.rds")

row_plot<-(sc_mac.markers[sc_mac.markers$gene %in% (tmp %>% filter(padj<=0.05)%>% filter(log2FoldChange>0))$gene ,] %>% filter(cluster=="Pf4 Mac"))$gene
row_plot<-c(row_plot,(sc_mac.markers[sc_mac.markers$gene %in% (tmp %>% filter(padj<=0.05)%>% filter(log2FoldChange<0))$gene,] %>% filter(cluster=="Vsig4 Mac"))$gene)


library(pheatmap)
rld <- rlog(dds, blind=F)
matrix <- assay(rld)[ c(rownames(res %>% filter(log2FoldChange>0.25)%>% filter(padj<=0.05)),rownames(res %>% filter(log2FoldChange<0)%>% filter(padj<=0.05))), ]
matrix <- matrix - rowMeans(matrix)
annotation_data <- as.data.frame(colData(rld)[c("age")])


tmp2<-strsplit(rownames(matrix), split = "\\|")
tmp3<-c()
for(el in tmp2){
  tmp3<-c(tmp3,el[3])
}

rownames(matrix)<-tmp3
labs.row <- rownames(matrix)
labs.row[!(labs.row %in% row_plot)]<-""
png(filename = "~/bulk_yo/results/macro/heatmap.png", width = 7, height = 10, units = "in", res=600)

pheatmap(matrix, annotation_col=annotation_data,cluster_rows =F, scale = "row", color = viridis_pal()(50),labels_row = labs.row)

dev.off()



