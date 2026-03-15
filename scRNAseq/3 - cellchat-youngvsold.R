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


library(CellChat)
library(scales)


sc1<-subset(sc, condition=="young")
sc2<-subset(sc, condition=="old")


data1.input <- GetAssayData(sc1, assay = "RNA", slot = "data") # normalized data matrix
labels1 <- Idents(sc1)
meta1 <- data.frame(group = labels1, row.names = names(labels1)) # create a dataframe of the cell labels
cellchat1 <- createCellChat(object = data1.input, meta = meta1, group.by = "group")
cellchat1 <- addMeta(cellchat1, meta = meta1, meta.name = "labels")
cellchat1 <- setIdent(cellchat1, ident.use = "labels") # set "labels" as default cell identity


data2.input <- GetAssayData(sc2, assay = "RNA", slot = "data") # normalized data matrix
labels2 <- Idents(sc2)
meta2 <- data.frame(group = labels2, row.names = names(labels2)) # create a dataframe of the cell labels
cellchat2 <- createCellChat(object = data2.input, meta = meta2, group.by = "group")
cellchat2 <- addMeta(cellchat2, meta = meta2, meta.name = "labels")
cellchat2 <- setIdent(cellchat2, ident.use = "labels") # set "labels" as default cell identity




CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat1@DB <- CellChatDB.use
cellchat2@DB <- CellChatDB.use


# subset the expression data of signaling genes for saving computation cost
cellchat1 <- subsetData(cellchat1) # This step is necessary even if using the whole database
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)


cellchat1 <- computeCommunProb(cellchat1)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)

cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)




# subset the expression data of signaling genes for saving computation cost
cellchat2 <- subsetData(cellchat2) # This step is necessary even if using the whole database
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)


cellchat2 <- computeCommunProb(cellchat2)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)

cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)



cellchat1<-netAnalysis_computeCentrality(cellchat1)
cellchat2<-netAnalysis_computeCentrality(cellchat2)

object.list <- list(young = cellchat1, old = cellchat2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat, "~/young-oldv2/results/merge/merge-final2/cellchat/cellchat.rds")
saveRDS(object.list, "~/young-oldv2/results/merge/merge-final2/cellchat/object-list.rds")

cellchat<-readRDS("~/young-oldv2/results/merge/merge-final2/cellchat/cellchat.rds")

object.list<-readRDS("~/young-oldv2/results/merge/merge-final2/cellchat/object-list.rds")

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
png(filename = "~/young-oldv2/results/merge/merge-final2/cellchat/Nombre_force_interation_totale.png",  width = 5, height = 10, units = "in",res=600)
gg1 + gg2
dev.off()


png(filename = "~/young-oldv2/results/merge/merge-final2/cellchat/Nombre_interation_youngvsold.png", width = 5, height = 5, units = "in",res=600)
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(1,2), color.use = cols, margin=c(0.3,0.3,0.3,0.3), vertex.label.cex  = 0.00000000001, title.name = "")
dev.off()
png(filename = "~/young-oldv2/results/merge/merge-final2/cellchat/Force_interation_youngvsold.png", width = 500, height = 500, units = "px")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(1,2))
dev.off()





png(filename = "~/young-oldv2/results/merge/merge-final2/cellchat/Heatmap_force_interation_youngvsold.png", width = 500, height = 500, units = "px")


netVisual_heatmap(cellchat, measure = "weight",comparison = c(1,2))
dev.off()

png(filename = "~/young-oldv2/results/merge/merge-final2/cellchat/Force_interation_young-old-individuel.png", width = 1000, height = 500, units = "px")
weight.max <- getMaxWeight(object.list,  attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Force of interactions - ", names(object.list)[i]))
}
dev.off()


png(filename = "~/young-oldv2/results/merge/merge-final2/cellchat/rank_network_interation_youngvsold.png", width = 10, height = 15, units = "in",res=600)


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(1,2),color.use = c("#173C70","#6C6E38"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(1,2))
gg1 + gg2
dev.off()


png(filename = "~/young-oldv2/results/merge/merge-final2/cellchat/il6fig.png", width = 6, height = 3, units = "in",res=600)

pathways.show <- c("IL6") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.use = cols,color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

dev.off()

library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.10.0
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> The new InteractiveComplexHeatmap package can directly export static 
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 


png(filename = "~/young-oldv2/results/merge/merge-final2/cellchat/outgoing_signal.png", width = 1000, height = 2000, units = "px")





pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20)

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


png(filename = "~/young-oldv2/results/merge/merge-final2/cellchat/incoming_signal.png", width = 1000, height = 2000, units = "px")


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20, color.heatmap = "GnBu")

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()





png(filename = "~/young-oldv2/results/merge/merge-final2/cellchat/signal_prot_youngvsold.png", width = 1500, height = 2000, units = "px")


netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45,signaling = pathways.show)
dev.off()

library(stringr)

pathways.show <- unique(c(cellchat@netP[["young"]][["pathways"]],cellchat@netP[["old"]][["pathways"]]))
for(path in pathways.show){
  dir.create(paste0("~/young-oldv2/results/merge/merge-final2/cellchat/signal/",path))
  dir.create(paste0("~/young-oldv2/results/merge/merge-final2/cellchat/signal/",path,"/genes"))
  
  try({
    weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = path) # control the edge weights across different datasets
    
  png(filename = paste0("~/young-oldv2/results/merge/merge-final2/cellchat/signal/",path,"/circle_youngvsold.png"), width = 10, height = 5, units = "in", res=600)
  
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = path) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = path, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(path, names(object.list)[i]), color.use = cols)
  }
  dev.off()

  })
  try(expr = {
    netVisual_aggregate(object.list[[1]],signaling =path, color.use = cols)
    png(filename = paste0("~/young-oldv2/results/merge/merge-final2/cellchat/signal/",path,"/young.png"), width = 5, height = 5, units = "in", res=600)
    netVisual_aggregate(object.list[[1]],signaling =path, color.use = cols)
    dev.off()
    png(filename = paste0("~/young-oldv2/results/merge/merge-final2/cellchat/signal/",path,"/contrib_young.png"), width = 5, height = 5, units = "in", res=600)
    
    print(netAnalysis_contribution(object.list[[1]], signaling = path))
    dev.off()
    tmp<-netAnalysis_contribution(object.list[[1]], signaling = path, return.data = T)
    extracted_text <- unique(unlist(str_extract_all(as.character(tmp[["LR.contribution"]][["name"]]), "\\w*\\d+\\w*")))
    for(gene in extracted_text){
      png(filename = paste0("~/young-oldv2/results/merge/merge-final2/cellchat/signal/",path,"/genes/",gene,".png"), width = 5, height = 5, units = "in", res=600)
      print(VlnPlot(sc,gene, split.by = "condition"))
      dev.off()
    }
  }
        )
  try(expr = {
    netVisual_aggregate(object.list[[2]],signaling =path, color.use = cols)
    
    png(filename = paste0("~/young-oldv2/results/merge/merge-final2/cellchat/signal/",path,"/old.png"), width = 5, height = 5, units = "in", res=600)
    netVisual_aggregate(object.list[[2]],signaling =path, color.use = cols)
    dev.off()
    png(filename = paste0("~/young-oldv2/results/merge/merge-final2/cellchat/signal/",path,"/contrib_old.png"), width = 5, height = 5, units = "in", res=600)
    
    print(netAnalysis_contribution(object.list[[2]], signaling = path))
    dev.off()
    tmp<-netAnalysis_contribution(object.list[[2]], signaling = path, return.data = T)
    extracted_text <- unique(unlist(str_extract_all(as.character(tmp[["LR.contribution"]][["name"]]), "\\w*\\d+\\w*")))
    for(gene in extracted_text){
      png(filename = paste0("~/young-oldv2/results/merge/merge-final2/cellchat/signal/",path,"/genes/",gene,".png"), width = 5, height = 5, units = "in", res=600)
      print(VlnPlot(sc,gene, split.by = "condition"))
      dev.off()
    }
  }
  
  )
  
}



path<-"IL6"
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = path) # control the edge weights across different datasets

png(filename = paste0("~/young-oldv2/results/merge/merge-final2/cellchat/circle_il6_youngvsold.png"), width = 10, height = 5, units = "in", res=600)

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = path) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = path, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(path, names(object.list)[i]), color.use = cols,vertex.label.cex=0.0000000000001)
}
dev.off()


path<-"TGFb"
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = path) # control the edge weights across different datasets

png(filename = paste0("~/young-oldv2/results/merge/merge-final2/cellchat/circle_tgfb_youngvsold.png"), width = 10, height = 5, units = "in", res=600)

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = path) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = path, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(path, names(object.list)[i]), color.use = cols,vertex.label.cex=0.0000000000001)
}
dev.off()


path<-"WNT"
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = path) # control the edge weights across different datasets

png(filename = paste0("~/young-oldv2/results/merge/merge-final2/cellchat/circle_wnt_youngvsold.png"), width = 10, height = 5, units = "in", res=600)

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = path) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = path, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(path, names(object.list)[i]), color.use = cols,vertex.label.cex=0.0000000000001)
}
dev.off()