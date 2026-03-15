library(Seurat)
library(dplyr)
library(ggplot2)
library(enrichR)

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

# library(clustree)
# clustree(sc, node_colour = "sc3_stability")




png(filename = "~/young-oldv2/results/merge/merge-final2/dimplot.png", width = 500, height = 500, units = "px")
DimPlot(sc, cols = cols)
dev.off()
png(filename = "~/young-oldv2/results/merge/merge-final2/dimplot-clean.png", width = 10, height = 10, units = "in",res = 600)
DimPlot(sc, cols = cols)+ NoLegend() + NoAxes() +ggtitle(NULL)
dev.off()
png(filename = "~/young-oldv2/results/merge/merge-final2/dimplot-clean-young.png", width = 10, height = 10, units = "in",res = 600)
DimPlot(subset(sc, condition=="young"), group.by = "condition",cols=c("#173C70"))+ NoLegend() + NoAxes()+ggtitle(NULL)
dev.off()
png(filename = "~/young-oldv2/results/merge/merge-final2/dimplot-clean-old.png", width = 10, height = 10, units = "in",res = 600)
DimPlot(subset(sc, condition=="old"), group.by = "condition",cols=c("#6C6E38"))+ NoLegend() + NoAxes()+ggtitle(NULL)
dev.off()



png(filename = "~/young-oldv2/results/merge/merge-final2/ft-ptprc.png", width = 5, height = 5, units = "in",res = 600)
FeaturePlot(sc, "Ptprc")+ NoLegend() + NoAxes()+ggtitle(NULL)
dev.off()
png(filename = "~/young-oldv2/results/merge/merge-final2/ft-ptprc-legend.png", width = 5, height = 5, units = "in",res = 600)
FeaturePlot(sc, "Ptprc")
dev.off()

png(filename = "~/young-oldv2/results/merge/merge-final2/ft-dcn.png", width = 5, height = 5, units = "in",res = 600)
FeaturePlot(sc, "Dcn")+ NoLegend() + NoAxes()+ggtitle(NULL)
dev.off()
png(filename = "~/young-oldv2/results/merge/merge-final2/ft-dcn-legend.png", width = 5, height = 5, units = "in",res = 600)
FeaturePlot(sc, "Dcn")
dev.off()



gene_dotplot<-c(
  "Pdgfra", "Ly6a",
  "Runx2","Alpl",
  "Scx","Fmod",
  "Mylpf","Acta1",
  "Pecam1","Cdh5",
  "Cd68","Adgre1",
  "Plac8","Ly6c2",
  "S100a9","Retnlg",
  "Cd3g","Trbc2",
  "Cma1","Kit"
)


png(filename = "~/young-oldv2/results/merge/merge-final2/dotplot.png", width = 8.5, height = 5, units = "in",res = 600)
DotPlot(sc,features=gene_dotplot,
cols = c ("lightgrey", "#D60A0A")) + ggtitle(NULL)
dev.off()

png(filename = "~/young-oldv2/results/merge/merge-final2/dotplot-split-clean-legend.png", width = 9, height = 7.5, units = "in",res = 600)
DotPlot(sc,features=gene_dotplot,split.by = "condition",
        cols = c ("black", "black"))
dev.off()









composition<-table(sc@active.ident, sc$condition)
write.csv(composition,"~/young-oldv2/results/merge/merge-final2/composition.csv")
openxlsx::write.xlsx(composition, file = '~/young-oldv2/results/merge/merge-final2/composition.xlsx') 

 





sc.markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sc@misc$markers<-sc.markers

write.csv(sc.markers,"~/young-oldv2/results/merge/merge-final2/markers.csv")
openxlsx::write.xlsx(sc.markers, file = '~/young-oldv2/results/merge/merge-final2/markers.xlsx', rowNames=TRUE) 





library(clusterProfiler)
library(org.Mm.eg.db)

yvso<-list()
yvso2<-list()
for(i in levels(sc@active.ident)){
  sctemp<-subset(sc, idents = i)
  sctemp<-SetIdent(sctemp, value = "condition")
  markers<-FindAllMarkers(sctemp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers,paste0("~/young-oldv2/results/merge/merge-final2/markers_youngvsold_",i,".csv"))
  yvso<-append(yvso, list( el = as.data.frame(markers)))
  names(yvso)[length(yvso)] <- i
}
openxlsx::write.xlsx(yvso, file = '~/young-oldv2/results/merge/merge-final2/markers_youngvsold.xlsx') 

# this part is to create enrichment dotplot between young and aged for fibroblasts
ck <- compareCluster(geneCluster = list(
  Prg4_down=(yvso$`Prg4 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Prg4_up=(yvso$`Prg4 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene,
  Myoc_down=(yvso$`Myoc Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Ugdh_down=(yvso$`Ugdh Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Cxcl12_down=(yvso$`Cxcl12 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Col3a1_down=(yvso$`Col3a1 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Ccl11_down=(yvso$`Ccl11 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Myoc_up=(yvso$`Myoc Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene,
  Ugdh_up=(yvso$`Ugdh Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene,
  Cxcl12_up=(yvso$`Cxcl12 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene,
  Col3a1_up=(yvso$`Col3a1 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene,
  Ccl11_up=(yvso$`Ccl11 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene
  ),
  fun = enrichGO,OrgDb=org.Mm.eg.db,keyType="SYMBOL",ont      = "BP")


ck <- compareCluster(geneCluster = list(
  Prg4_down=(yvso$`Prg4 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Prg4_up=(yvso$`Prg4 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene,
  Myoc_down=(yvso$`Myoc Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Myoc_up=(yvso$`Myoc Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene,
  Ugdh_down=(yvso$`Ugdh Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Ugdh_up=(yvso$`Ugdh Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene,
  Cxcl12_down=(yvso$`Cxcl12 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Cxcl12_up=(yvso$`Cxcl12 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene,
  Col3a1_down=(yvso$`Col3a1 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Col3a1_up=(yvso$`Col3a1 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene,
  Ccl11_down=(yvso$`Ccl11 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="young"))$gene,
  Ccl11_up=(yvso$`Ccl11 Fib` %>% filter(p_val_adj<=0.05)%>% filter(cluster=="old"))$gene
),
fun = enrichGO,OrgDb=org.Mm.eg.db,keyType="SYMBOL",ont      = "BP")
library(enrichplot)

source("~/young-oldv2/scripts/merge/v2/custom_dotplot.R")


png(filename = "~/young-oldv2/results/merge/merge-final2/go_yvso_fib2.png",width = 7, height = 5, units = "in",res = 600)
custom_dotplot(ck, showCategory=rev(unique(c(
  "wound healing",
  "connective tissue development",
  "collagen fibril organization",
  "collagen metabolic process",
  "regeneration",
  "biomineral tissue development",
  "cell chemotaxis",
  "oxidative phosphorylation",
  "complement activation",
  "lipid transport",
  "response to interleukin-6",
  "regulation of metallopeptidase activity"
)
)))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

# this part is to create enrichment dotplot between young and aged for macrophages
ck <- compareCluster(geneCluster = list(
  Pf4_down = (yvso$`Pf4 Mac` %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "young"))$gene,
  Vsig4_down = (yvso$`Vsig4 Mac` %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "young"))$gene,
  Cd74_down = (yvso$`Cd74 Mac` %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "young"))$gene,
  Pf4_up = (yvso$`Pf4 Mac` %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "old"))$gene,
  Vsig4_up = (yvso$`Vsig4 Mac` %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "old"))$gene,
  Cd74_up = (yvso$`Cd74 Mac` %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "old"))$gene
),
fun = enrichGO,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP")

png(filename = "~/young-oldv2/results/merge/merge-final2/go_yvso_mac.png",width = 7, height = 5, units = "in",res = 600)
custom_dotplot(ck,showCategory=rev(c(
  "immune response signaling pathway",
  "regulation of inflammatory response",
  "cellular response to molecule of bacterial origin",
  "response to lipopolysaccharide",
  "leukocyte migration",
  "leukocyte chemotaxis",
  "wound healing",
  "phagocytosis",
  "interleukin-10 production",
  "maintenance of location",
  "iron ion homeostasis",
  "fibroblast proliferation",
 "response to interleukin-4",
  
  "oxidative phosphorylation",
  "ATP biosynthetic process",
  "complement activation",
  "calcium ion homeostasis"
)
))
dev.off()

nb_de_genes<-as.data.frame(matrix(0, ncol=3,nrow=0))
colnames(nb_de_genes)<-c("cluster","nb","age")
degenes_y<-c()
degenes_o<-c()
degenes<-as.data.frame(matrix(0,ncol=length(levels(sc@active.ident)), nrow = 0))
colnames(degenes)<-levels(sc@active.ident)

for(clustern in levels(sc@active.ident)){
  tmp_y<-(yvso[[clustern]] %>% filter(p_val_adj<=0.05) %>% filter(cluster=="young") )$gene
  tmp_o<-(yvso[[clustern]] %>% filter(p_val_adj<=0.05) %>% filter(cluster=="old") )$gene
  nb_de_genes[nrow(nb_de_genes)+1,]=c(clustern,length(tmp_y),"young")
  nb_de_genes[nrow(nb_de_genes)+1,]=c(clustern,length(tmp_o),"old")
  degenes_y<-c(degenes_y,tmp_y)
  degenes_o<-c(degenes_o,tmp_o)
  for(x in tmp_y){
    degenes[x,clustern]<--yvso[[clustern]][x,"avg_log2FC"]
  }
  for(x in tmp_o){
    degenes[x,clustern]<-+yvso[[clustern]][x,"avg_log2FC"]
  }
}

degenes[is.na(degenes)]<-0

png(filename = "~/young-oldv2/results/merge/merge-final2/aging-marks/heatmap_de_youngvsold.png",width = 8, height = 15, units = "in",res = 600)
pheatmap::pheatmap(degenes, angle_col = 90, scale = "none",border = NA, cluster_rows=T,cluster_cols=F, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),breaks = seq(-1,1,0.04), show_rownames=F)
dev.off()




nb_de_genes<-as.data.frame(matrix(0, ncol=3,nrow=0))
colnames(nb_de_genes)<-c("cluster","nb","age")
degenes_y<-c()
degenes_o<-c()
degenes<-as.data.frame(matrix(0,ncol=length(levels(sc@active.ident)), nrow = 0))
colnames(degenes)<-levels(sc@active.ident)

for(clustern in levels(sc@active.ident)){
  tmp_y<-(yvso[[clustern]] %>% filter(p_val_adj<=0.05) %>% filter(cluster=="young") )$gene
  tmp_o<-(yvso[[clustern]] %>% filter(p_val_adj<=0.05) %>% filter(cluster=="old") )$gene
  nb_de_genes[nrow(nb_de_genes)+1,]=c(clustern,length(tmp_y),"young")
  nb_de_genes[nrow(nb_de_genes)+1,]=c(clustern,length(tmp_o),"old")
  degenes_y<-c(degenes_y,tmp_y)
  degenes_o<-c(degenes_o,tmp_o)
}

degenes[is.na(degenes)]<-0


nb_de_genes$nb<-as.integer(nb_de_genes$nb)
png(filename = "~/young-oldv2/results/merge/merge-final2/aging-marks/nb_de_youngvsold.png",width = 16, height = 10, units = "in",res = 600)

ggplot(nb_de_genes,aes(x=cluster, y = ifelse(age == "young", -nb, nb), fill=age))+ 
  geom_bar(stat="identity", position="identity")+
  scale_y_continuous(limits = c(- max(nb_de_genes$nb), max(nb_de_genes$nb))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))+
  scale_x_discrete(limits=unique(nb_de_genes$cluster))+
  scale_fill_manual(values = c("firebrick3","navy") )+
  coord_flip()+ theme(panel.background = element_rect(fill="white", colour="white"))
dev.off()



sc$ident<-sc@active.ident
sc_fib<-subset(sc,ident %in% fib_names)

png(filename = "~/young-oldv2/results/merge/merge-final2/fib/dimplot-fib-young.png", width = 5, height = 5, units = "in",res = 600)
DimPlot(subset(sc_fib, condition=="young"), cols=c("#263238","#B71C1C","#006064","#AA00FF","#F57F17","#1A237E"))+ NoLegend() + NoAxes()+ggtitle(NULL)
dev.off()
png(filename = "~/young-oldv2/results/merge/merge-final2/fib/dimplot-fib-old.png", width = 5, height = 5, units = "in",res = 600)
DimPlot(subset(sc_fib, condition=="old"), cols=c("#263238","#B71C1C","#006064","#AA00FF","#F57F17","#1A237E"))+ NoLegend() + NoAxes()+ggtitle(NULL)
dev.off()


sc_fib.markers <- FindAllMarkers(sc_fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cells_to_remove <- c(
  "young_ACTATCTAGGGTATAT-1", "young_AGGGTTTAGTGGAAGA-1", "young_AGTGCCGGTCGCACGT-1", "young_ATCGTGATCAAAGGTA-1", "young_CTCAAGACATGTGACT-1",
  "young_GAATCACCAGACCAGA-1", "old_AAAGGTATCGACGAGA-1", "old_AAGACTCGTCAGACTT-1", "old_ACTTCCGGTGCACAAG-1", "old_ATCGGATCAGTAGAAT-1",
  "old_ATTCTACTCGCGCTGA-1", "old_CATACCCAGGATTTCC-1", "old_CCATAAGTCCGTAGTA-1", "old_CCCATTGTCACGGAGA-1", "old_CCCGGAAGTTTAGAGA-1",
  "old_CTAACCCGTAAGTCAA-1", "old_CTCCCAACAGATGCGA-1", "old_GCACTAAGTTTCGGCG-1", "old_GTGACGCTCATGGAGG-1", "old_TAAGCGTTCCACAAGT-1",
  "old_TACCCGTGTAGACGGT-1", "old_TATACCTTCATACGAC-1", "old_TCCTCCCGTACTGAGG-1", "old_TCGGGACTCGTACACA-1", "old_TCGGTCTGTGGAAGTC-1",
  "old_TGAATGCTCACGTCCT-1", "old_TGCGATATCCTGCTAC-1", "old_TGTTCTACACCTGCGA-1", "old_TTCTCTCAGAGAGAAC-1","young_AAAGGGCTCGGTAACT-1",
  "young_ACAGAAATCTCATTTG-1", "young_CGTTGGGAGAGCAGTC-1", "young_GACGTTACATCCCACT-1", "young_TCCTTTCAGACTCCGC-1", "young_TGGGCTGAGACCTCCG-1",
  "old_AAGGTAACAGTAGAGC-1",   "old_ATCTCTACAATTGTGC-1",   "old_CACGGGTTCTTGCAGA-1" 
)
sc_fib <- subset(sc_fib, cells = setdiff(Cells(sc_fib), cells_to_remove))








ck_top50 <- compareCluster(geneCluster = list(
  `Prg4 Fib` = (sc_fib.markers %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "Prg4 Fib") %>% slice_min(order_by = p_val_adj,n=50))$gene,
  `Myoc Fib` = (sc_fib.markers %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "Myoc Fib")%>% slice_min(order_by = p_val_adj,n=50))$gene,
  `Ugdh Fib` = (sc_fib.markers %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "Ugdh Fib")%>% slice_min(order_by = p_val_adj,n=50))$gene,
  `Cxcl12 Fib` = (sc_fib.markers %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "Cxcl12 Fib")%>% slice_min(order_by = p_val_adj,n=50))$gene,
  `Col3a1 Fib` = (sc_fib.markers %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "Col3a1 Fib")%>% slice_min(order_by = p_val_adj,n=50))$gene,
  `Ccl11 Fib` = (sc_fib.markers %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "Ccl11 Fib")%>% slice_min(order_by = p_val_adj,n=50))$gene
),
fun = enrichGO,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP")
openxlsx::write.xlsx(as.data.frame(ck_top50), file = paste0('~/young-oldv2/results/merge/merge-final2/fib/enrich_top50.xlsx') )




custom_dotplot(ck_top50, ,showCategory=rev(unique(c(
  "extracellular matrix organization",
  "extracellular structure organization",
  "Notch signaling pathway",
  "collagen metabolic process",
  "cell-substrate adhesion",
  "connective tissue development",
  
  "regulation of leukocyte migration",
  "actin filament polymerization",
  "regulation of inflammatory response",
  "chemokine production",
  "response to endoplasmic reticulum stress",
  "response to interleukin-6",
  "reactive oxygen species metabolic process",
  "positive regulation of signal transduction by p53 class mediator",
  "regulation of leukocyte differentiation",
  "cell chemotaxis",
  "collagen fibril organization",
  "response to acid chemical",
  "cartilage development",
  "cell-substrate adhesion",
  "regulation of angiogenesis"
)
)))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
#taxonomy
fls<-c("F13a1","Col22a1","Clic5","Gchfr","Tspan15","Tmem196","Dlx3","Itga6","Fut9","Rab37")
f1<-c("S100a10","Anxa8","Tspo","Crip1","Ociad2","Anxa2","Slurp1","Txn1","Igfbp6","Crip2")
f2<-c("Cxcl12","Ctsb","Pcolce","Gpnmb","Gas6","Dcn","Serping1","Serpina3n","Rarres2","C4b")
f3<-c("Cxcl1","Lpl","Dpep1","Ccl11","Acvr2a","Ntrk2","Gsn","Col4a1","Fst","Bmper")
f4<-c("Pi16","Car8","Efhd1","Anxa3","Aif1l","Mfap5","Cd248","Tek","Zyx","Sema3c")
f5<-c("C7","Abcc9","Fmo2","Rbp1","F3","Kcnj8","Sparcl1","Gdf10","Kitl","Cygb")
f6<-c("Hsd11b1","Vtn","Aldh1a2","Ret","Atp1a2","Cldn15","Ccl11","Hmcn2","D630033O11Rik","Prss12")

sc_fib<-AddModuleScore(sc_fib, features = list(fls=fls,f1=f1,f2=f2,f3=f3,f4=f4,f5=f5,f6=f6))

png(filename = "~/young-oldv2/results/merge/merge-final2/fib/ft-taxonomy.png", width = 12, height = 10, units = "in",res = 600)
FeaturePlot(sc_fib,c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7"), max.cutoff = "q95", min.cutoff = "q5")
dev.off()


openxlsx::write.xlsx(sc_fib.markers, file = '~/young-oldv2/results/merge/merge-final2/fib/markers.xlsx', rowNames=TRUE) 



source("~/young-oldv2/scripts/merge/v2/DoMultiBarHeatmap.R")
cols.use <- list(ident=c(cols[1:6]), condition=c("#5493E0", "#A5A853"))

png(filename = "~/young-oldv2/results/merge/merge-final2/fib/heatmap_top20pvalgroupyvso.png", width = 10, height = 10, units = "in",res = 600)
DoMultiBarHeatmap(sc_fib, top20_pval$gene,group.by='ident', additional.group.by = 'condition',draw.lines = F,label = F, cols.use = cols.use)
dev.off()









macro_names<-c("Pf4 Mac", "Vsig4 Mac", "Cd74 Mac")
sc_mac<-subset(sc,ident%in%macro_names)
sc_mac.markers <- FindAllMarkers(sc_mac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


openxlsx::write.xlsx(sc_mac.markers, file = '~/young-oldv2/results/merge/merge-final2/mac/markers.xlsx', rowNames=TRUE) 
saveRDS(sc_mac.markers, '~/young-oldv2/results/merge/merge-final2/mac/markers.rds')


k0<-c("Aqp1","Fxyd2","Tppp3","Cd9","Lyve1","Anxa1","Cyb5r3","Ifi27l2a")
k1<-c("H2-Ab1","H2-Eb1","H2-Aa","Cd74","H2-DMa","Cd52", "Mgl2","S100a11","Gm2a","Lsp1")
k2<-c("Ccl8","Ccl2","Ccl7","Cxcl2","Dusp1","Zfp36","Junb","Pf4","Fos","Atf3")
k3<-c("Vsig4","Sparc","Ctsd","Lyz2","Srgn","S100b","Fn1","Apoe","Wfdc17","Hexb")
k4<-c("Stmn1","Ube2c","Birc5","Tubb5","Hmgb2","Tuba1b","Tubb4b","Ptma","H2afz","Vim")
k5<-c("Acp5","Ctsk","Mmp9","Atp6v0d2","Atp6v1g1","Atp6v1b2","S100a4","Rplp0","Atp6v0e","Clec12a")
sc_mac<-AddModuleScore(sc_mac, features = list(k0,k1,k2,k3,k4,k5), name="kronke")

png(filename = "~/young-oldv2/results/merge/merge-final2/mac/ft-culemann.png", width = 10, height = 10, units = "in",res = 600)
FeaturePlot(sc_mac,c("kronke1","kronke2","kronke3","kronke4","kronke5","kronke6"), max.cutoff = "q95", min.cutoff = "q5")
dev.off()



ck_top50 <- compareCluster(geneCluster = list(
  Pf4 = (sc_mac.markers %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "Pf4 Mac")%>% slice_min(order_by = p_val_adj,n=50))$gene,
  Vsig4 = (sc_mac.markers %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "Vsig4 Mac")%>% slice_min(order_by = p_val_adj,n=50))$gene,
  Cd74 = (sc_mac.markers %>% filter(p_val_adj <= 0.05) %>% filter(cluster == "Cd74 Mac")%>% slice_min(order_by = p_val_adj,n=50))$gene
),
fun = enrichGO,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP")
openxlsx::write.xlsx(as.data.frame(ck_top50), file = paste0('~/young-oldv2/results/merge/merge-final2/mac/enrich_top5°.xlsx') )

library(enrichplot)

options(enrichplot.colours = c("red","blue"))




png(filename = "~/young-oldv2/results/merge/merge-final2/mac/gomac.png", width = 5, height = 8, units = "in",res = 600)
custom_dotplot(ck_top50, showCategory=rev(unique(c(
  "receptor-mediated endocytosis",
  "cell chemotaxis",
  "wound healing",
  "response to extracellular stimulus",
  "endothelial cell chemotaxis",
  "fibroblast proliferation",
  "apoptotic cell clearance",
  "positive regulation of angiogenesis",
  "positive regulation of vasculature development",
  "regulation of inflammatory response",
  "regulation of inflammatory response",
  "regulation of angiogenesis",
  "regulation of peptidase activity",
  "glycosaminoglycan catabolic process",
  "Tenoitin sulfate metabolic process",
  "sulfur compound catabolic process",
  "aminoglycan catabolic process",
  "cell junction maintenance",
  "tissue remodeling",
  "positive regulation of T cell activation",
  "antigen processing and presentation of exogenous peptide antigen via MHC class II",
  "antigen processing and presentation of peptide antigen via MHC class II",
  "peptide antigen assembly with MHC class II protein complex",
  "lymphocyte mediated immunity",
  "leukocyte cell-cell adhesion"
))
)
)

dev.off()







msc<-c("Hp","Cxcl12","Adipoq","Lpl","Esm1","Gdpd2","Gas6","Cxcl14","Serping1","Dpep1","Grem1","Sfrp4","Pappa","Kitl","Chrdl1","Wisp2","1500009L16Rik","Fbln5","Vcam1","Serpina3g","Cyp1b1","Lepr","Agt","Ebf3","Il34","Fst","Arrdc4","Kng1","Kng2","Fgf7","Slc26a7","H2-Q10","March1","Ms4a4d","Wdr86","Il1rn","Serpina3c","C2","Igfbp4","Cebpa","Serpina12","Cbln1","Cdh11","Mme","Apoe","Tmem176b","Csf1","Ibsp","H2-K1","Serpine2","Tmem176a","H2-D1","Cldn10","Tnc","Trf","Cdh2","Vegfc","Gpr88","Ggt5","Gm16685","S1pr3","Igfbp5","Gpx3","Pdzrn4","Nnmt","Cd302","Ghr","Rarres2","Gm4951","Zeb2","Angpt1","B2m","Il7","Gpm6b","Runx1","Cd1d1","Csmd1","Ebf1","Cp","Ptx3","Steap4","Plpp3","Fstl1","Igsf3","Ackr4","Shroom3","Mdk","H2-Q7","Ccl2","Fmo2","Gja1","Eef1a1","Mmp23","Ptprd","Ccl9","Ifitm3","Tnfrsf19","Pdgfrb","C1ra","Mylk")
sc<-AddModuleScore(sc,features = list(msc),name = 'msc')






sen_mayo<-c("Acvr1b","Ang","Angpt1","Angptl4","Areg","Axl","Bex3","Bmp2","Bmp6","C3","Ccl1","Ccl2","Ccl20","Ccl24","Ccl26","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Cd55","Cd9","Csf1","Csf2","Csf2rb","Cst10","Ctnnb1","Ctsb","Cxcl1","Cxcl10","Cxcl12","Cxcl16","Cxcl2","Cxcl3","Cxcr2","Dkk1","Edn1","Egf","Egfr","Ereg","Esm1","Ets2","Fas","Fgf1","Fgf2","Fgf7","Gdf15","Gem","Gmfg","Hgf","Hmgb1","Icam1","Icam5","Igf1","Igfbp1","Igfbp2","Igfbp3","Igfbp4","Igfbp5","Igfbp6","Igfbp7","Il10","Il13","Il15","Il18","Il1a","Il1b","Il2","Il6","Il6st","Il7","Inha","Iqgap2","Itga2","Itpka","Jun","Kitl","Lcp1","Mif","Mmp13","Mmp10","Mmp12","Mmp13","Mmp14","Mmp2","Mmp3","Mmp9","Nap1l4","Nrg1","Pappa","Pecam1","Pgf","Pigf","Plat","Plau","Plaur","Ptbp1","Ptger2","Ptges","Rps6ka5","Scamp4","Selplg","Sema3f","Serpinb3a","Serpine1","Serpine2","Spp1","Spx","Timp2","Tnf","Tnfrsf11b","Tnfrsf1a","Tnfrsf1b","Tubgcp2","Vegfa","Vegfc","Vgf","Wnt16","Wnt2")
sc<-AddModuleScore(sc,features = list(sen_mayo),name = 'sen_mayo')

png(filename = "~/young-oldv2/results/merge/merge-final2/aging-marks/senmayo.png", width = 10, height = 10, units = "in",res = 600)
VlnPlot(sc,features="sen_mayo1", group.by = "condition", pt.size=0, cols = c(c("#3C8B38","#DC7314")))
dev.off()



sc$ident<-sc@active.ident

sc_mac<-subset(sc,ident==mac_names)

gene_bulk<-readRDS("~/bulk_yo/results/macro/degenes.rds")
gene_down<-(gene_bulk %>% filter(log2FoldChange>0.25)%>% filter(padj<=0.05))$gene
gene_up<-(gene_bulk %>% filter(log2FoldChange<0.25)%>% filter(padj<=0.05))$gene
sc_mac<-AddModuleScore(sc_mac,features = list(gene_down),name = 'down',nbin=1)
sc_mac<-AddModuleScore(sc_mac,features = list(gene_up),name = 'up', nbin = 1)





sc_stroma<-subset(sc, ident %in% c("Prg4 Fib","Myoc Fib", "Ugdh Fib","Cxcl12 Fib","Col3a1 Fib","Ccl11 Fib","Osteoblasts","Teno"))

png(filename = paste0("~/young-oldv2/results/merge/merge-final2/vln-plot-fib/msc.png"), width = 5, height = 3,5, units = "in",res = 600)
print(VlnPlot(sc_stroma, "msc1", pt.size = 0, cols=cols[1:8]))
dev.off()




vlnplotfib <- function(gene){
  png(filename = paste0("~/young-oldv2/results/merge/merge-final2/vln-plot-fib/",gene,".png"), width = 5, height = 3,5, units = "in",res = 600)
  print(VlnPlot(subset(sc,ident==fib_names), gene, pt.size = 0, cols=cols[1:6]))
  dev.off()
}

markers<-c("Prg4", 
           "Tnfaip6",
           "Myoc",
           "Lum","Ccl11","Thy1","Cd34",
           "Gsn","C3","Cxcl1","Ptgs2","Has1","Fbn1","Dpp4","Cd34","Cd248","Pi16","Myoc",
           "Lum","Cxcl14","Ptn","Emb","Cdh11","Col22a1","Tspan15","Clic5","Htra4","F13a1","Thbs1","Ccl11","Smoc2","Fmo2","Bmp5",
           "Chad","Cilp2","Fmod","Thbs4","Scx","Cdkn1a",
           "Col22a1","Htra4","Clic5","Hbegf","Tspan15",
           "Apod","Clec3b","Cd34","Cemip2","Cemip","Chsy3",
           "Il6","Fbn1","Mfap5","Pi16","Slit3","Efemp1",
           "Ugdh","Tnfaip6","Cxcl1","C3","Ifi205","Has1","Has2","Irak3","Myc","Ptgs2",
           "Cxcl12","Lum","Cxcl14","C4b", "Ccl19", "Lepr",
           "Col3a1","Col1a1","C7","Col5a1","Postn","Adamts2",
           "Ccl11","Thbs1","Smoc2",
           "Txnip","Crip1","Disp1","Sema3c","Sdk1","Pi16","Fosb","Junb","Nfkb1","Itga5","Myoc","Fbn1"
           
)
lapply(markers,vlnplotfib)










customvlnplot <- function(gene){
  png(filename = paste0("~/young-oldv2/results/merge/merge-final2/vlnplot/",gene,".png"), width = 10, height = 3,5, units = "in",res = 600)
  print(VlnPlot(sc, gene, pt.size = 0, cols=cols))
  dev.off()
}

markers<-c("Pdgfra","Ly6a",#fibro
           "Cd68","Adgre1",#macro
            "Tnnt3","Actn3","Tpm2", #muscle
           "Il4","Il13","Mcpt4",#mast cells
           "Emcn","Esam","Cldn5", #EC
           "Itk","Il2rb","Trbc2","Trbc1",#lympho
           "S100a8","S100a9","Csf3r","Ccr1",#neutro
           "Chil3","Ccr2","Csf2rb",#mono
           "Sox9","Clu","Cilp2", #Teno
           "Runx2","Alpl","Ostn", #osteo
           "Plac8","Chil3","Apoc2","Adgre4","F10","Ms4a4c","Ly6c2",
           "Trbc2",
           "Kit","Gata2",
           "Trem2",
           "Arg2","Trem1","Retnlg",
           "Scx","Thbs4","Fmod"
           
           
           
)
lapply(markers,customvlnplot)




chakarov<-c("Folr2","Cd163","F13a1","C4b","Ednrb","Fcgrt","Gas6","Cd38","Ninj1","Lyve1","Ccl24","Ltc4s","Timd4","Pf4","Selenop","Rcn3","Ifitm3","Clec10a","Cfh","Smpdl3a","Cd5l","Tcn2","Cfp","Ifitm2","Colec12","Serpinb6a","Pepd","Ctsb","Msr1","Emilin2","Slc9a3r2","Marco","Gypc","Vsig4","Slc40a1","Timp2","Dmpk","Pmp22","Selenbp1","Fam234a","Ceacam1","Snx2","Cd36","Hgsnat","Cxcl13","Nrp1","Mlxipl","Ly6e","Dab2","Fgfr1","Ap2a2","Wwp1","Slc28a2","Ecm1","C6","Tslp","Thbd","Hbb-bs","C1qc","Sptbn1","Itsn1","Ptafr","Aldh2","Stab1","Stard8","Mrc1","Ifi204","Glul","Adam15","Pltp","Ddx60","Fam43a","Fam213b","Maf","Gprc5b","Igfbp4","Nid2","Trf","Blvrb","Plekhg5","Ap2m1","Pros1","Wnt2","Emp1","Hal","Alox5","S1pr1","Tmem8","Gstp1","Rnase4","Tns1","Sesn1","Frmd4b","Ms4a8a","Mpp1","Rasgrp3","Grap","Acly","C1qb","Smagp")
chakarovmhc<-c("Cd9","Mpeg1","Lyz1","Cx3cr1","Cxcl16","Rgs1","Mmp12","H2-Aa","H2-Ab1","Atp6v0d2","H2-DMa","H2-DMb1","Cd74","Axl","Olfml3","H2-Eb1","Il1b","Tmem119","Plbd1","Col14a1","Pald1","Plet1","Car4","Gm2a","Tgfbr1","Dnase1l3","Hexb","Ccr2","Cd72","B4galnt1","Ccrl2","Clec4n","Sema4d","Slamf7","Lpcat2","Spp1","Krt79","Abcg1","Lmna","Cyp4f18","Cyb561a3","Cd300c2","Tlr2","F11r","Mmp13","Lpl","Cd83","Tmem154","Cadm1","Ramp1","Fgl2","Sh2d1b1","Mmp14","Zmynd15","Adam8","Glipr1","Ctss","Syngr2","Psap","Net1","Lsp1","Il1rn","F7","Fam129a","Sirpa","Amz1","Il1rl2","Lrrfip1","Lcp1","Cd44","Cd4","Gpnmb","Egr2","Btg2","Clic4","Olr1","Cytip","Slamf8","Cd300lf","Furin","Rin3","Parvg","Tnfaip3","Sgk1","Fabp1","Siglecf","Skil","Ptgs2","Cd2","Ezr","Fxyd5","Nceh1","Itgax","Slpi","Dusp5","Ctsk","Krt19","Clec7a","Ccnd2","Inmt")
sc<-AddModuleScore(sc,features = list(chakarov),name = 'chakarov')
sc<-AddModuleScore(sc,features = list(chakarovmhc),name = 'chakarovmhc')


png(filename = "~/young-oldv2/results/merge/merge-final2/macro_chakarovlyve1top50.png", width = 7, height = 7, units = "in",res = 600)

SCpubr::do_BoxPlot(subset(sc,ident %in% mac_names), "chakarov1", boxplot.width = 2,colors.use = c("Pf4 Mac"=cols[11], 
                                                                       "Vsig4 Mac"=cols[12], 
                                                                       "Cd74 Mac"=cols[13]
))
dev.off()

png(filename = "~/young-oldv2/results/merge/merge-final2/macro_chakarovmhctop50.png", width = 7, height = 7, units = "in",res = 600)

SCpubr::do_BoxPlot(subset(sc,ident %in% mac_names), "chakarovmhc1", boxplot.width = 2,colors.use = c("Pf4 Mac"=cols[11], 
                                                                                                  "Vsig4 Mac"=cols[12], 
                                                                                                  "Cd74 Mac"=cols[13]
))
dev.off()

 

vlnplotmac <- function(gene){
  png(filename = paste0("~/young-oldv2/results/merge/merge-final2/vln-plot-mac/",gene,".png"), width = 5, height = 3,5, units = "in",res = 600)
  print(VlnPlot(subset(sc,ident==mac_names), gene, pt.size = 0, cols=cols[11:13]))
  dev.off()
}

markers<-c("Pf4","Aqp1","Fxyd2",
           "Cx3cr1","Vsig4","Sparc",
           "Cd74","Ccr2","H2-Ab1","H2-Eb1",
           "Timd4","Lyve1","Folr2","Mrc1",
           "Cd163","F13a1","Ccl24","Gatm","Trem2","Pmepa1","S100b","Mgl2","Ltb4r1")
lapply(markers,vlnplotmac)
