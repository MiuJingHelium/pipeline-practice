
library("Seurat")
library("tidyverse")
library("Matrix")
library("readr")
library("MAST")

#args = commandArgs(trailingOnly=T)
#indir=args[1]
#object=args[2]
#load(paste0(indir,"/",object))
load("Pan_CD4_processed.Robj")
cluster.averages <- AverageExpression(Pan_CD4, return.seurat = FALSE)
write.table(cluster.averages$originalexp, paste0("Pan_CD4_avgexp_harmony.tsv"), sep="\t", quote=F, col.names=NA)
load("Pan_CD8_processed.Robj")
cluster.averages <- AverageExpression(Pan_CD8, return.seurat = FALSE)
write.table(cluster.averages$originalexp, paste0("Pan_CD8_avgexp_harmony.tsv"), sep="\t", quote=F, col.names=NA)

