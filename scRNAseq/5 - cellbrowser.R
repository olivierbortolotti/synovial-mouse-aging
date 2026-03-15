library(Seurat)
library(dplyr)


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


sc.markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sc@misc$markers<-sc.markers
sc$orig.ident<-NULL
sc$old.ident<-NULL
sc$age<-sc$condition
sc$condition<-NULL
sc$SCT_snn_res.0<-NULL
sc$SCT_snn_res.0.05<-NULL
sc$SCT_snn_res.0.1<-NULL
sc$SCT_snn_res.0.15<-NULL
sc$SCT_snn_res.0.2<-NULL
sc$SCT_snn_res.0.25<-NULL
sc$SCT_snn_res.0.3<-NULL
sc$SCT_snn_res.0.35<-NULL
sc$SCT_snn_res.0.4<-NULL
sc$SCT_snn_res.0.45<-NULL
sc$SCT_snn_res.0.5<-NULL
sc$SCT_snn_res.0.55<-NULL
sc$SCT_snn_res.0.6<-NULL
sc$SCT_snn_res.0.65<-NULL
sc$SCT_snn_res.0.7<-NULL
sc$seurat_clusters<-NULL

sc$ident_age<-paste(sc$ident,sc$age,sep="_")


saveRDS(sc,"~/young-oldv2/results/merge/merge-final2/sc.rds")
