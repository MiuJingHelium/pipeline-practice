
library("Seurat")
library("tidyverse")
library("Matrix")
library("readr")
library("MAST")

args = commandArgs(trailingOnly=T)
indir=args[1]
object=args[2]
load(paste0(indir,"/",object))
cluster.averages <- AverageExpression(whole, return.seurat = FALSE)
write.table(cluster.averages$originalexp, paste0(indir, "/avg.exp.tsv"), sep="\t", quote=F, col.names=NA)
