library(DESeq2)
library(ggplot2)
library(gprofiler2)
library(dplyr)

countData <- read.csv('~/bulk_yo/data/Courties_G_03_OF_vs_YF_Filtered_Raw_Counts.txt', header = TRUE, sep = "\t")

row.names(countData)<-countData[,1]
countData<-countData[,2:ncol(countData)]
countData[] <- lapply(countData, as.numeric)




age<-c(rep("old",4), rep("young",3))
metaData<-data.frame(age)


dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~age)

dds <- DESeq(dds)
vsdata <- vst(dds, blind=FALSE)


pcaData <- plotPCA(vsdata, intgroup="age", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
colors <- c("young" = "#173C70", "old" = "#6C6E38")

png(filename = "~/bulk_yo/results/fibro/pca.png", width = 5, height = 2, units = "in", res=600)
ggplot(pcaData, aes(PC1, PC2, color=age)) +
  geom_point(size=3) +
  scale_color_manual(values=colors) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal( )+
  ggtitle("PCA Plot")+theme_bw()+
   expand_limits(x = c(-10,15), y=c(-10,10))
dev.off()

countn<-counts(dds, normalized=T)
res <- results(dds)
res$gene<-rownames(res)
res<-as.data.frame(res)
openxlsx::write.xlsx(res %>% filter(padj<=0.05), file = '~/bulk_yo/results/fibro/degene.xlsx', rowNames=TRUE) 





tmp<- res
tmp2<-strsplit(rownames(tmp), split = "\\|")
tmp3<-c()
for(el in tmp2){
  tmp3<-c(tmp3,el[3])
}

tmp$gene<-tmp3

tmp$fcsign <- sign(tmp$log2FoldChange)
tmp$logP=-log10(tmp$pvalue)
tmp$metric= tmp$logP/tmp$fcsign
y<-tmp[,c("gene","metric")]
filtered <- na.omit(y)

#input for gsea
write.table(filtered ,file="~/bulk_yo/results/fibro/rank_fibro.rnk",quote=F,sep="\t",row.names=F, col.names = F)













resp_ox_stress<-c(
  "ENSMUSG00000021109|15251|Hif1a",#
  "ENSMUSG00000026822|16819|Lcn2",#
  "ENSMUSG00000027737|26570|Slc7a11",#
  "ENSMUSG00000061878|20698|Sphk1",#
  "ENSMUSG00000032802|76650|Srxn1",#
  "ENSMUSG00000020250|50493|Txnrd1",#
  
  "ENSMUSG00000022548|11815|Apod",#
  "ENSMUSG00000002985|11816|Apoe",#
  "ENSMUSG00000002602|26362|Axl",#
  "ENSMUSG00000078566|12176|Bnip3",#
  "ENSMUSG00000032323|13070|Cyp11a1",#
  "ENSMUSG00000033268|99439|Duox1",#
  "ENSMUSG00000031616|13617|Ednra",#
  "ENSMUSG00000031400|14381|G6pdx",#
  "ENSMUSG00000063856|14775|Gpx1",#
  "ENSMUSG00000060803|14870|Gstp1",#
  "ENSMUSG00000030562|50490|Nox4",#
  "ENSMUSG00000054580|18779|Pla2r1",#
  "ENSMUSG00000029167|19017|Ppargc1a",#
  "ENSMUSG00000024953|54683|Prdx5",#
  "ENSMUSG00000047250|19224|Ptgs1",#
  "ENSMUSG00000038393|56338|Txnip"
)

matrix <- assay(rld)[resp_ox_stress, ]
matrix <- matrix - rowMeans(matrix)
matrix<-matrix[,c("Young_1_Fibro", "Young_2_Fibro", "Young_4_Fibro","Old_1_Fibro","Old_2_Fibro","Old_3_Fibro","Old_4_Fibro")]

annotation_data <- as.data.frame(colData(rld)[c("age")])

tmp2<-strsplit(rownames(matrix), split = "\\|")
tmp3<-c()
for(el in tmp2){
  tmp3<-c(tmp3,el[3])
}

rownames(matrix)<-tmp3

png(filename = "~/bulk_yo/results/fibro/heatmap_final-resp_ox_stress.png", width = 5, height = 20, units = "in", res=600)

pheatmap(matrix, annotation_col=annotation_data,cluster_rows =F, cluster_cols = F,scale = "row", color =  colorRampPalette(c('blue',"white", 'red'))(50),border_color = "black",cellwidth = 10, cellheight = 10)

dev.off()









inflammation<-c(
  "ENSMUSG00000035385|20296|Ccl2",#
  "ENSMUSG00000035373|20306|Ccl7",#
  "ENSMUSG00000056501|12608|Cebpb",
  "ENSMUSG00000029380|14825|Cxcl1",#
  "ENSMUSG00000021508|57266|Cxcl14",#
  "ENSMUSG00000029371|20311|Cxcl5",#
  "ENSMUSG00000021109|15251|Hif1a",#
  "ENSMUSG00000070942|107527|Il1rl2", #
  "ENSMUSG00000074115|20208|Saa1",#
  "ENSMUSG00000040026|20210|Saa3",#
  "ENSMUSG00000040152|21825|Thbs1",#
  "ENSMUSG00000019850|21929|Tnfaip3",#
  "ENSMUSG00000044534|59289|Ackr2",#
  "ENSMUSG00000073418|12268|C4b",#
  "ENSMUSG00000019122|20308|Ccl9",#
  "ENSMUSG00000027221|76969|Chst1",#
  "ENSMUSG00000023078|55985|Cxcl13",#
  "ENSMUSG00000019278|13479|Dpep1",#
  "ENSMUSG00000017493|16010|Igfbp4",#
  "ENSMUSG00000042228|17096|Lyn",#
  "ENSMUSG00000028751|26970|Pla2g2e",#
  "ENSMUSG00000047250|19224|Ptgs1",#
  "ENSMUSG00000042312|20196|S100a13",#
  "ENSMUSG00000038264|20361|Sema7a",#
  "ENSMUSG00000031596|11988|Slc7a2",#
  "ENSMUSG00000028599|21938|Tnfrsf1b"#

)

matrix <- assay(rld)[inflammation, ]
matrix <- matrix - rowMeans(matrix)
matrix<-matrix[,c("Young_1_Fibro", "Young_2_Fibro", "Young_4_Fibro","Old_1_Fibro","Old_2_Fibro","Old_3_Fibro","Old_4_Fibro")]

annotation_data <- as.data.frame(colData(rld)[c("age")])

tmp2<-strsplit(rownames(matrix), split = "\\|")
tmp3<-c()
for(el in tmp2){
  tmp3<-c(tmp3,el[3])
}

rownames(matrix)<-tmp3

png(filename = "~/bulk_yo/results/fibro/heatmap_final-inflam.png", width = 5, height = 20, units = "in", res=600)

pheatmap(matrix, annotation_col=annotation_data,cluster_rows =F, cluster_cols = F,scale = "row", color =  colorRampPalette(c('blue',"white", 'red'))(50),border_color = "black",cellwidth = 10, cellheight = 10)

dev.off()








ecm_org<-c(
  "ENSMUSG00000001506|12842|Col1a1",
  "ENSMUSG00000029661|12843|Col1a2",
  "ENSMUSG00000026043|12825|Col3a1",
  "ENSMUSG00000031502|12826|Col4a1",
  "ENSMUSG00000067158|12829|Col4a4",
  "ENSMUSG00000031274|12830|Col4a5",
  "ENSMUSG00000026837|12831|Col5a1",
  "ENSMUSG00000026042|12832|Col5a2",
  "ENSMUSG00000043719|245026|Col6a6",
  "ENSMUSG00000056174|329941|Col8a2",
  "ENSMUSG00000026147|12839|Col9a1",#
  "ENSMUSG00000027966|12814|Col11a1",#
  "ENSMUSG00000024330|12815|Col11a2",
  "ENSMUSG00000028339|12819|Col15a1",
  "ENSMUSG00000000402|54156|Egfl6",
  "ENSMUSG00000029675|13717|Eln",#
  "ENSMUSG00000041559|14264|Fmod",#
  "ENSMUSG00000024529|16948|Lox",#
  "ENSMUSG00000034205|94352|Loxl2",
  "ENSMUSG00000000693|16950|Loxl3",
  "ENSMUSG00000025185|67573|Loxl4",#
  "ENSMUSG00000060572|17150|Mfap2",#
  "ENSMUSG00000042436|76293|Mfap4",#
  "ENSMUSG00000043613|17392|Mmp3",#
  "ENSMUSG00000028226|17389|Mmp16",#
  "ENSMUSG00000032011|21838|Thy1",
  "ENSMUSG00000052353|80982|Cemip",#
  "ENSMUSG00000063564|237759|Col23a1",
  "ENSMUSG00000020467|216616|Efemp1",#
  "ENSMUSG00000037370|18605|Enpp1",
  "ENSMUSG00000000392|14089|Fap",
  "ENSMUSG00000026193|14268|Fn1",#
  "ENSMUSG00000002020|16997|Ltbp2",
  "ENSMUSG00000024940|16998|Ltbp3",#
  "ENSMUSG00000000901|17385|Mmp11",#
  "ENSMUSG00000002603|21803|Tgfb1"#

)

matrix <- assay(rld)[ecm_org, ]
matrix <- matrix - rowMeans(matrix)
matrix<-matrix[,c("Young_1_Fibro", "Young_2_Fibro", "Young_4_Fibro","Old_1_Fibro","Old_2_Fibro","Old_3_Fibro","Old_4_Fibro")]

annotation_data <- as.data.frame(colData(rld)[c("age")])

tmp2<-strsplit(rownames(matrix), split = "\\|")
tmp3<-c()
for(el in tmp2){
  tmp3<-c(tmp3,el[3])
}

rownames(matrix)<-tmp3

png(filename = "~/bulk_yo/results/fibro/heatmap_final-ecm_org.png", width = 5, height = 20, units = "in", res=600)

pheatmap(matrix, annotation_col=annotation_data,cluster_rows =F, cluster_cols = F,scale = "row", color =  colorRampPalette(c('blue',"white", 'red'))(50),border_color = "black",cellwidth = 10, cellheight = 10)

dev.off()

