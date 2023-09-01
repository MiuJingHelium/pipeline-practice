install.packages("SoupX")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DropletUtils")
BiocManager::install("harmony")
install.packages('Seurat')

library("harmony")
library(SoupX)
library(Matrix)
library(data.table)
library(DropletUtils)
library(Seurat)
library(tidyverse)

outdir = "../SoupX_outs/"
indir = "./outs/"

samples = list.files(path = indir)
outs <- list()

ctr = 1
for (sample in samples) {
  ## load sample data
  toc <- Seurat::Read10X(paste(indir, sample, "/filtered_feature_bc_matrix", sep = ""))
  tod <- Seurat::Read10X(paste(indir, sample, "/raw_feature_bc_matrix", sep = ""))
  meta <- fread(paste(indir, sample, "/analysis/clustering/graphclust/clusters.csv", sep = ""))
  #because clustering was not performed by cellranger for these two samples,cluter info are loaded from Seurat runs;In this case, need to make sure the data is unfiltered.
  #meta <- read.table(file = paste0("./analysis/",sample,"/data_for_plot.tsv"),sep = "\t",header = T)
  clusters <- meta$Cluster
  names(clusters) <- meta$Barcode
  
  ## run soupx
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  sc <- setClusters(sc, clusters)
  sc <- estimateSoup(sc)
  sc <- autoEstCont(sc)
  out <- adjustCounts(sc)
  outs[[ctr]] <- out
  ctr = ctr + 1
  ## save result
  #DropletUtils:::write10xCounts(paste(outdir, sample, sep = ""), out)
}


ALAW_3 = CreateSeuratObject(outs[[1]])#This will load one sample;for multiple sample, a loop can be used (or better alternatives); If processed interactively, sample name may be recommended for naming the Seurat Object.
ALAW_4 = CreateSeuratObject(outs[[2]])
whole <- merge(x = ALAW_3, y = ALAW_4, orig.ident=samples, add.cell.ids=samples)#Ideally, the orig.ident should be also defined using the sample name.

rm(outs)
rm(out)
rm(meta)
rm(toc)
rm(tod)
rm(ALAW_3)
rm(ALAW_4)
gc()

mito.genes <- c(grep("^MT-|^mt-", rownames(x = whole), value = T))
percent.mito <- Matrix::colSums(x = GetAssayData(object = whole, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = whole, slot = 'counts'))
whole[['percent.mito']] <- percent.mito

#whole@meta.data$orig.ident = substring(names(whole@active.ident),1,6)

jpeg(paste0(outdir,"mito.percentage.jpeg"))

ggplot(whole@meta.data)+
  geom_histogram(aes(x=percent.mito,y = after_stat(width*density*100)),bins = 100)+
  facet_wrap(vars(orig.ident),ncol=2)+
  geom_vline(xintercept = 0.1,col="red") +
  ylab("proportion")+
  xlab("percentage of mitochondrial RNA")
dev.off()

jpeg(paste0(outdir,"n_RNAfeatures.jpeg"))
ggplot(whole@meta.data)+
  geom_histogram(aes(x = log2(nFeature_RNA),y = after_stat(width*density*100)),bins = 100)+
  facet_wrap(vars(orig.ident),ncol=2)+
  geom_vline(xintercept = 8.5,col="red") +
  ylab("proportion")+
  xlab("number of features")
dev.off()

cells.before <- length(colnames(x = whole))
print(paste0("cells before filtering: ",cells.before))

whole <- subset(x = whole, subset = percent.mito < 0.1 & nFeature_RNA_log > 8.5)

## how many cells are left
cells.after <- length(colnames(x = whole))
print(paste0("cells after filtering: ",cells.after))

whole <- SCTransform(whole, vars.to.regress = "percent.mito", verbose = FALSE)

whole <- FindVariableFeatures(object = whole, selection.method = 'vst', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))

whole@assays$RNA@var.features <- whole@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|MAMU-A", whole@assays$RNA@var.features)]
whole <- ScaleData(object = whole, features = VariableFeatures(object = whole), vars.to.regress = c("nCount_RNA", "percent.mito"))

whole <- RunPCA(object = whole, verbose=FALSE)

whole@meta.data$Subject <- gsub("_.*", "", whole@meta.data$orig.ident)#create the subject info for the samples (whether they are from the same subject or not)
whole <- RunHarmony(whole, "Subject")#Run harmony based on the subject

whole <- whole %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
save(whole, file = paste(outdir, "whole_object.Robj", sep = ""))

## identify markers
whole.markers <- FindAllMarkers(object = whole,
                                only.pos = TRUE,
                                min.pct = 0.15,
                                thresh.use = 0.15)
write.table(whole.markers, paste(outdir, "markers.tsv", sep = ""), sep="\t", quote=F, row.names=F)





