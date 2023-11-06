
library(Seurat)
library(tidyverse)
#library(harmony)
library(gridExtra)

path = "./SCEPrep/"
load(paste0(path,"CD8_reprocessed.Robj"))

gobs <- list(
  p1 = DimPlot(CD8,reduction = "pca",group.by = "Manually_curated_celltype"),
  p2 = DimPlot(CD8,reduction = "tsne",group.by = "Manually_curated_celltype"),
  p3 = DimPlot(CD8,reduction = "umap",group.by = "Manually_curated_celltype"),
  p4 = DimPlot(CD8,reduction = "pca",group.by = "Donor"),
  p5 = DimPlot(CD8,reduction = "tsne",group.by = "Donor"),
  p6 = DimPlot(CD8,reduction = "umap",group.by = "Donor")
)

g <- grid.arrange(grobs = gobs, ncol = 3)
ggsave(paste0(path,"dimplots.pdf"),plot=g, height = 60, width = 90, units = "cm")

#cluster.averages <- AverageExpression(CD8, return.seurat = FALSE)
#write.table(cluster.averages$originalexp, paste0(path,"CD8_avgexp_harmony.tsv"), sep="\t", quote=F, col.names=NA)

#save(CD8,file = "CD8_reprocessed.Robj")
