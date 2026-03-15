library(Seurat)
library(dplyr)
library(ggplot2)
library(enrichR)
library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)

sc<- readRDS("~/young-oldv2/results/merge/sc_postsctransform_ccr.rds")
ident<-"SCT_snn_res.0.45"
sc<-SetIdent(sc, value = ident)

new.cluster.ids <- c("Myoc Fib", "Ugdh Fib","Cxcl12 Fib", "Pf4 Mac", "Col3a1 Fib", "Vsig4 Mac", "Teno", "Cd74 Mac","Neutrophils","Prg4 Fib","Ccl11 Fib","Osteoblasts","Muscle","Monocytes","Lymphocytes", "EC","Mast")
current.cluster.ids<-c(0:16)
sc@active.ident <- plyr::mapvalues(x = sc@active.ident, from = current.cluster.ids, to = new.cluster.ids)
fib_names<-c("Prg4 Fib","Myoc Fib", "Ugdh Fib","Cxcl12 Fib","Col3a1 Fib","Ccl11 Fib")
mac_names<-c("Pf4 Mac", "Vsig4 Mac", "Cd74 Mac")

levels(x = sc) <- c("Prg4 Fib","Myoc Fib", "Ugdh Fib","Cxcl12 Fib","Col3a1 Fib","Ccl11 Fib","Osteoblasts","Teno","Muscle","EC","Pf4 Mac", "Vsig4 Mac", "Cd74 Mac","Monocytes", "Neutrophils",  "Lymphocytes", "Mast")
sc$condition <- factor(sc$condition, levels = c("young", "old"))
cols<-c("#263238","#B71C1C","#006064","#AA00FF","#F57F17","#1A237E","#adad24","#01579B","#880E4F","#855605","#1B5E20","#827717","#E65100","#6ac46d","#260202","#00BCD4","#FF5722")


sc$ident<-sc@active.ident

regulonAUC <- importAUCfromText(  "~/young-oldv2/data/scenic/auc_mtx.csv")
AUCmat <- AUCell::getAUC(regulonAUC)
sc[['pyscenicAUC']] <- CreateAssayObject(data = AUCmat)



DefaultAssay(sc) <- 'pyscenicAUC'
sc<-ScaleData(sc)




sc_fibro<-subset(sc, ident %in% c("Prg4 Fib","Myoc Fib", "Ugdh Fib","Cxcl12 Fib","Col3a1 Fib","Ccl11 Fib"))

sc.markers <- FindAllMarkers(sc_fibro, only.pos = TRUE,logfc.threshold = 0.005)
openxlsx::write.xlsx(sc.markers, file = '~/young-oldv2/results/merge/merge-final2/scenic2/markers-fibro.xlsx', rowNames=TRUE) 
library(dplyr)
library(pheatmap)

# Your Seurat object
# seurat_obj <- readRDS("your_seurat_object.rds")  # Optional: load your Seurat object

# Set the identity class (e.g., by cluster)

# Genes of interest
genes_of_interest <- (sc.markers %>% filter(p_val_adj<=0.05) )$gene

# Aggregate expression data (pseudobulk) by averaging within clusters
pseudobulk_expr <- AverageExpression(sc_fibro, features = genes_of_interest, return.seurat = FALSE, assays = "pyscenicAUC")$pyscenicAUC

# Optional: scale the expression for visualization
scaled_expr <- t(scale(t(pseudobulk_expr)))  # Z-score per gene

# Draw heatmap
png(filename = '~/young-oldv2/results/merge/merge-final2/scenic2/heatmap_fibro.png',width = 10, height = 15, units = "in",res = 600)
pheatmap(scaled_expr, 
         cluster_rows = F, 
         cluster_cols = F, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Pseudobulk Heatmap by Cluster")

dev.off()


png(filename = '~/young-oldv2/results/merge/merge-final2/scenic2/heatmap_fibrohi.png',width = 20, height = 20, units = "in",res = 600)
pheatmap(scaled_expr, 
         cluster_rows = F, 
         cluster_cols = F, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Pseudobulk Heatmap by Cluster")

dev.off()


