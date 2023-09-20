library(SoupX)#Soup removal
library(DoubletFinder)#Doublet removal
library(miQC)#mito percentage filtering
library(data.table) #For fread
library(SeuratWrappers) #For miQC
library(Seurat)
library(Matrix)
library(tidyverse)
library(cowplot)

args <- commandArgs(trailingOnly = T) #only arguments following the script are listed
#indir, outdir, sample

indir = args[1] #Can be made into a command-line input
outdir = args[2]
sample = args[3] #command-line input

print(args)
#samples_dir = list.files(path = indir, full.names = T)

#############list directories for SoupX ##############
toc <- Seurat::Read10X(paste(indir,"/", sample, "/outs/filtered_feature_bc_matrix", sep = ""))
tod <- Seurat::Read10X(paste(indir,"/", sample, "/outs/raw_feature_bc_matrix", sep = ""))
meta <- fread(paste(indir, "/",sample, "/outs/analysis/clustering/graphclust/clusters.csv", sep = ""))

#prep clusterinfo
clusters <- meta$Cluster
names(clusters) <- meta$Barcode

############# run Soupx ##########
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- setClusters(sc, clusters)
sc <- estimateSoup(sc)
sc <- autoEstCont(sc)
out <- adjustCounts(sc)

rm(toc)
rm(tod)
rm(meta)
gc()

print("SoupX finished.")
############  mitochondrial RNA removal and low nFeatureRNA filtering #########

#create object
obj = CreateSeuratObject(out) #it is possible to assign names based on sample name with some downstream technical challenges
#output of eval is possibly non-writable

obj[["percent.mt"]] = PercentageFeatureSet(object = obj, pattern = "^MT-|^mt-")
jpeg(paste0(outdir,"/percent_mt-nFeature_RNA.jpeg"))
FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")
dev.off()

#obj <- RunMiQC(obj, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.75,
#               model.slot = "flexmix_model")

#jpeg(paste0(outdir,"/MiQCplot.jpeg"))
#PlotMiQC(obj, color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")
#dev.off()

#jpeg(paste0(outdir, "/miQC_keep.jpeg"))
#PlotMiQC(obj, color.by = "miQC.keep")
#dev.off()


obj$nFeature_RNA_log2 = log2(obj$nFeature_RNA)
obj_filtered <- subset(obj, subset = percent.mt < 5 & nFeature_RNA_log2 > 8.5)

cells.before <- length(colnames(x = obj))
print(cells.before)

cells.after <- length(colnames(x = obj_filtered))
print(cells.after)



gc()

############## Visualizing and checking #############
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "mvp", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:30) #use top 30 PC

obj_filtered <- NormalizeData(obj_filtered)
obj_filtered <- FindVariableFeatures(obj_filtered, selection.method = "mvp", nfeatures = 3000)
obj_filtered <- ScaleData(obj_filtered)
obj_filtered <- RunPCA(obj_filtered)
obj_filtered <- RunUMAP(obj_filtered, dims = 1:30)

jpeg(paste0(outdir,"/filtering comparison.jpeg"))
plt1 <- FeaturePlot(obj,features = "percent.mt",reduction = "umap") + ggtitle("percent.mt unfiltered")
plt2 <- FeaturePlot(obj_filtered,features = "percent.mt",reduction = "umap") + ggtitle("percent.mt filtered")
plt3 <- FeaturePlot(obj,features = "nFeature_RNA",reduction = "umap") + ggtitle("nGenes unfiltered")
plt4 <- FeaturePlot(obj_filtered,features =  "nFeature_RNA",reduction = "umap") + ggtitle("nGenes filtered")
plot_grid(plt1,plt2,plt3,plt4,labels = "AUTO",nrow = 2)#need "cowplot"
dev.off()

rm(obj)
gc()
##################### DoubleFinder #############

# Find significant PCs
stdv <- obj_filtered[["pca"]]@stdev
sum.stdv <- sum(obj_filtered[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                     percent.stdv[2:length(percent.stdv)]) > 0.1), 
	             decreasing = T)[1] + 1
min.pc <- min(co1, co2)

print(paste0("minPC is ",min.pc))

obj_filtered <- RunUMAP(obj_filtered, dims = 1:30)#parameter values are updated to be consistent with actual analysis
obj_filtered <- FindNeighbors(object = obj_filtered, dims = 1:30)              
obj_filtered <- FindClusters(object = obj_filtered, resolution = 0.9)

print("Finished pre-processing.")

sweep.res.list <- paramSweep_v3(obj_filtered, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

optimal.pk = as.numeric(levels(bcmvn$pK)[which.max(bcmvn$BCmetric)])

print(paste0("Finished finding optimal pk:",  optimal.pk))

annotations <- obj_filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp.poi <- round(optimal.pk * nrow(obj_filtered@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

print("Finding doublets...")

obj_filtered <- doubletFinder_v3(seu = obj_filtered, 
                                 PCs = 1:30, 
                                 pK = optimal.pk,
                                 nExp = nExp.poi.adj)
print("Finished finding doublets. Start processing...")

metadata <- obj_filtered@meta.data

colnames(metadata)[(length(colnames(metadata)))] <- "doublet_finder"
obj_filtered@meta.data <- metadata
obj_filtered <- subset(obj_filtered, doublet_finder == "Singlet")

print("Finished doublets removal. Printing current information and storing object...")
obj_filtered
save(obj_filtered, file = paste0(outdir,"/",sample, "_filtered.Robj"))

sessionInfo()

