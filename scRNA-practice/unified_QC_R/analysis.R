library(Seurat)
library(dplyr)
library(Matrix)

## GATHERING DATA TOGETHER
#dataFolder <- "/scratch1/fs1/martyomov/carisa/GzmK_practice/results/align"

#dataFolders[[1]] <- paste0(dataFolder,"/Aged_Liver/outs/filtered_feature_bc_matrix")

######## Version note ########
# the directory for loading  #
# the feature, barcodes, and #
# matrix should be provided  #
# via the command line.      #
##############################

# In this script, each dataset is loaded independently for QC. i.e. no merge step. Thus, essentially only
# one directory is needed --> no loop or apply() needed.

# Consider including an argument that indicate whether merged() is used?
# Assume snakemake pipeline is used, so no need to worry about user behavior
# e.g. analysis.R (merge) (sample) indir
 
args <- commandArgs(trailingOnly=T)
# Assume command goes:
# Rscript example.R (merge) (sample) indir
# Either do merge and indir -> list.files(indir) -> merge
# OR do sample and indir -> paste0(indir, sample, suffix)
print(args)
suffix = "/outs/filtered_feature_bc_matrix"
idx = 1
dataFolders <- as.list(paste0(args[idx+1],"/",args[idx],suffix))
if (args[idx] == "merge"){
        indir = args[idx+1]
        sampleNames <- list.files(path = indir)
        dataFolders <- as.list(paste0(indir,"/",sampleNames,suffix))
}
print(dataFolders)



fdata <- list()
fdata[[1]] <- Read10X(dataFolders[[1]])

whole <- CreateSeuratObject(counts=fdata[[1]], project = "A_L003")
## Number of cells before
cells.before <- length(colnames(x = whole))
## NORMALIZATION
mito.genes <- c(grep("^MT-", rownames(x = whole), value = T),
                grep("^mt-", rownames(x = whole), value = T))
percent.mito <- Matrix::colSums(x = GetAssayData(object = whole, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = whole, slot = 'counts'))
whole[['percent.mito']] <- percent.mito

whole <- subset(x = whole, subset = percent.mito <= 0.05)
whole <- NormalizeData(object = whole, normalization.method = "LogNormalize", scale.factor = 10000)
whole <- FindVariableFeatures(object = whole, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
whole <- ScaleData(object = whole, features = VariableFeatures(object = whole), vars.to.regress = c("nCount_RNA", "percent.mito"))
gc()

## PCA
whole <- RunPCA(object = whole,
                features =  VariableFeatures(object = whole),
                verbose=FALSE)


## TSNE
whole <- RunTSNE(object = whole, dims = 1:20)
whole <- RunUMAP(object = whole, dims = 1:20)

## CLUSTERING
whole <- FindNeighbors(object = whole, dims = 1:20)
whole <- FindClusters(object = whole, resolution = 0.6)

DimPlot(object = whole, reduction="tsne")

## SAVING DATA FOR EXPLORER
expData <- GetAssayData(object = whole, slot = 'data')
save(expData, file="expData.Rda")

dataForPlot <- as.data.frame(whole@reductions$tsne@cell.embeddings)
dataForPlot$Sample <- whole@meta.data$orig.ident
dataForPlot$Cluster <-  Idents(object = whole)
dataForPlot$nUmi <- whole@meta.data$nCount_RNA
dataForPlot$nGene <- whole@meta.data$nFeature_RNA
dataForPlot$nUmiLog2 <- log2(whole@meta.data$nCount_RNA)
dataForPlot$nGeneLog2 <- log2(whole@meta.data$nFeature_RNA)

write.table(dataForPlot, "data_for_plot.tsv", sep="\t", quote=F)

## FINDING ANS SAVING MARKERS
whole.markers <- FindAllMarkers(object = whole,
                                only.pos = TRUE,
                                min.pct = 0.10,
                                logfc.threshold = 0.10)
write.table(whole.markers, "markers.tsv", sep="\t", quote=F, row.names=F)

## SAVING
save(whole, file = "whole_object.Robj")

## Number of cells after
cells.after <- length(colnames(x = whole))
print(paste0("cells.before:",cells.before))
print(paste0("cells.after:",cells.after))
print(paste0("cell.diff:", cells.before-cells.after))
