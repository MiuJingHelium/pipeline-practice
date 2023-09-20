library("Seurat")
library("tidyverse")
library("Matrix")
library("readr")
library("MAST")


args = commandArgs(trailingOnly = T)
indir = args[1]
#output into the same directory


#Pre-requisite: have a fully integrated, processed, and annotated Seurat object
load(paste0(indir, "whole_object.Robj")) 

### Assume output into the working directory ###

#### Creating data for plot ######
dataForPlot <- as.data.frame(whole@reductions$tsne@cell.embeddings)
dataForPlot$CellName <- rownames(dataForPlot)
dataForPlot$Sample <- whole@meta.data$orig.ident
dataForPlot$Cluster <-  Idents(whole)
#Cluster names can be changed in case DE is performed at a different population resolution
dataForPlot$nUmi <- whole@meta.data$nCount_RNA
dataForPlot$nGene <- whole@meta.data$nFeature_RNA
dataForPlot$nUmiLog2 <- log2(whole@meta.data$nCount_RNA)
dataForPlot$nGeneLog2 <- log2(whole@meta.data$nFeature_RNA)
write.table(dataForPlot, paste0(indir,"./data_for_plot.tsv"), sep="\t", quote=F,row.names = F)

#####################################

##### Creating average expression table ########

cluster.averages <- AverageExpression(whole, return.seurat = FALSE)
write.table(cluster.averages$RNA, paste0(indir, "avg.exp.tsv"), sep="\t", quote=F, col.names=NA)

###############################################


## FINDING ANS SAVING MARKERS
whole.markers <- FindAllMarkers(object = whole,
                                only.pos = TRUE,
                                min.pct = 0.10,
                                thresh.use = 0.10)
write.table(whole.markers, "markers.tsv", sep="\t", quote=F, row.names=F)

##################################################

