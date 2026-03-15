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




exprMat <- sc@assays$RNA@data
cellInfo <- sc@meta.data

loci1 <- which(rowSums(exprMat) > 1*.01*ncol(exprMat))
exprMat_filter <- exprMat[loci1, ]

add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

loom <- build_loom("~/young-oldv2/results/merge/merge-final2/sc.loom", dgem=exprMat_filter)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)





docker run -it --rm  -v /home/debian/scenic2:/data aertslab/pyscenic:0.12.0 arboreto_with_multiprocessing.py /data/sc.loom /data/mm_mgi_tfs.txt --num_workers 32 -o /data/adj.csv --method grnboost2





docker run -it --rm  -v /home/debian/scenic2:/data aertslab/pyscenic:0.12.0 pyscenic ctx /data/adj.csv /data/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /data/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname /data/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl.txt --expression_mtx_fname /data/sc.loom --mode "custom_multiprocessing" --output /data/regulons.csv --num_workers 32




docker run -it --rm -v /home/debian/scenic2:/data aertslab/pyscenic:0.12.0 pyscenic aucell /data/sc.loom /data/regulons.csv -o /data/auc_mtx.csv --num_workers 32