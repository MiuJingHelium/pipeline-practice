

library(Seurat)
library(tidyverse)
library(harmony)
library(gridExtra)
load("NKILC_TissueSubset.Robj")


table(NKILC_TissueSubset@meta.data$Manually_curated_celltype)


#DimPlot(NKILC_TissueSubset,reduction = "umap",group.by = "Organ")


## Re-scale everything


NKILC_TissueSubset <- NormalizeData(object = NKILC_TissueSubset, normalization.method = "LogNormalize", scale.factor = 10000)


NKILC_TissueSubset <- FindVariableFeatures(object = NKILC_TissueSubset, selection.method = 'mean.var.plot', mean.cutoff = c(0.01, Inf), dispersion.cutoff = c(0.5, Inf))
gc()


g <- VariableFeaturePlot(NKILC_TissueSubset)
ggsave(filename="NKILC_TissueSubset_VariableFeaturePlot.pdf",plot = g, height = 20,width = 30, units = "cm")


NKILC_TissueSubset@assays$originalexp@var.features <-  NKILC_TissueSubset@assays$originalexp@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|MAMU-A", NKILC_TissueSubset@assays$originalexp@var.features)]
NKILC_TissueSubset <- ScaleData(object = NKILC_TissueSubset, features = VariableFeatures(object = NKILC_TissueSubset), vars.to.regress = c("nCount_originalexp", "percent.mito"))

g <- ElbowPlot(NKILC_TissueSubset,ndims = 40) + geom_hline(yintercept = 2,col="red") + geom_hline(yintercept = 1.5, col ="blue") + geom_hline(yintercept = 1, col ="purple")+ labs(legend =  c( "elbow plot" = "black", "std = 2" = "red", "std = 1.5" = "blue","std = 1" = "purple"))
ggsave(filename="NKILC_TissueSubset_ElbowPlot.pdf", plot = g, height = 20, width = 30, units = "cm")


stdv <- NKILC_TissueSubset[["pca"]]@stdev
sum.stdv <- sum(NKILC_TissueSubset[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                     percent.stdv[2:length(percent.stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min.pc <- min(co1, co2)
print(min.pc)

### Top 30 PCs

NKILC_TissueSubset <- RunPCA(object = NKILC_TissueSubset,
                features =  VariableFeatures(object = NKILC_TissueSubset),
                dims = 1:30)



p1 <- DimPlot(NKILC_TissueSubset,reduction = "pca",group.by = "Sex")

p2 <- DimPlot(NKILC_TissueSubset,reduction = "pca",group.by = "Donor")

p3 <- DimPlot(NKILC_TissueSubset,reduction = "pca",group.by = "Organ")

p4 <- DimPlot(NKILC_TissueSubset,reduction = "pca",group.by = "Tissue_type")

p5 <- DimPlot(NKILC_TissueSubset,reduction = "pca",group.by = "Manually_curated_celltype")

g <- grid.arrange(as.list(p1,p2,p3,p4,p5),nrow = 2)
ggsave("NKILC_TissueSubset_PCAPlots.pdf",plot = g, height = 30, width = 60, units="cm")


NKILC_TissueSubset <- RunHarmony(object = NKILC_TissueSubset, group.by.vars = c("Donor"), assay.use = "originalexp", max.iter.harmony = 20)

PBMC_NK_markersG = c("IFNG","NCAM1","CCR7","IL2RB","IL7R","IL32","XCL1","GZMK","GZMH","GZMB","KLRB1","KLRC1","KLRC2","KLRG1","PRF1","CXCR1","CX3CR1","CD160","FCGR3A","MKI67")

PBMC_NK_markersTF = c("IFNG","NCAM1","MYC","TCF7","EOMES","BACH2","JUND","FOSB","TBX21","ZEB2","PRDM1","ZBTB38")

##### Top 30 PCs

NKILC_TissueSubset <- RunUMAP(NKILC_TissueSubset, dims = 1:30, reduction = "harmony")
gc()
NKILC_TissueSubset <- RunTSNE(NKILC_TissueSubset, dims = 1:30, reduction = "harmony")
gc()
NKILC_TissueSubset <- FindNeighbors(NKILC_TissueSubset, dims = 1:30, reduction = "harmony")
gc()
NKILC_TissueSubset <- FindClusters(object = NKILC_TissueSubset, reduction = "harmony", dims = 1:30,
    resolution = 0.5, print.output = 0, save.SNN = TRUE)

p1 <- DimPlot(NKILC_TissueSubset,reduction = "tsne",group.by = "Donor")
p2 <- DimPlot(NKILC_TissueSubset,reduction = "umap",group.by = "Donor")
p3 <- DimPlot(NKILC_TissueSubset,reduction = "umap",group.by = "Manually_curated_celltype")
p4 <- DimPlot(NKILC_TissueSubset,reduction = "tsne",group.by = "Manually_curated_celltype")
p5 <- DimPlot(NKILC_TissueSubset,reduction = "tsne")
p6 <- DimPlot(NKILC_TissueSubset,reduction = "umap")
g <- grid.arrange(as.list(p1,p2,p3,p4,p5,p6),nrow = 2)
ggsave("NKILC_TissueSubset_DimPlots_30PC.pdf",plot = g, height = 40, width = 60, units = "cm")

#FeaturePlot(CD4_TissueSubset,features = c("SELL","CCR7"),reduction = "umap")
#FeaturePlot(CD4_TissueSubset,features = c("CCR6","CCR4"),reduction = "umap")
#FeaturePlot(CD4_TissueSubset,features = c("CCR10","CCL5"),reduction = "umap")
#FeaturePlot(CD4_TissueSubset,features = c("CXCR5","CXCR3"),reduction = "umap")

#FeaturePlot(CD4_TissueSubset,features = c("GNLY","PRF1"),reduction = "umap")
#FeaturePlot(CD4_TissueSubset,features = c("GZMB","GZMK"),reduction = "umap")
#FeaturePlot(CD4_TissueSubset,features = c("CTLA4","PDCD1"),reduction = "umap")
#FeaturePlot(CD4_TissueSubset,features = c("LAG3","ISG15"),reduction = "umap")

#FeaturePlot(CD4_TissueSubset,features = c("IL2RA","IL4"),reduction = "umap")
#FeaturePlot(CD4_TissueSubset,features = c("IL10","IL7R"),reduction = "umap")
#FeaturePlot(CD4_TissueSubset,features = c("KLRB1","KLRF1"),reduction = "umap")


save(NKILC_TissueSubset, file = "NKILC_TissueSubset_30PC.Robj")

p1 <- DotPlot(NKILC_TissueSubset, features = PBMC_NK_markersG,cluster.idents = T) + RotatedAxis()
p2 <- DotPlot(NKILC_TissueSubset, features = PBMC_NK_markersTF,cluster.idents = T) + RotatedAxis()
g <- grid.arrange(as.list(p1,p2),nrow = 1)
ggsave("NKILC_TissueSubset_Dotplots_30PC.pdf",plot = g, height = 30, width = 60, units = "cm")

#### Top 20


NKILC_TissueSubset <- RunUMAP(NKILC_TissueSubset, dims = 1:20, reduction = "harmony")
gc()
NKILC_TissueSubset <- RunTSNE(NKILC_TissueSubset, dims = 1:20, reduction = "harmony")
gc()
NKILC_TissueSubset <- FindNeighbors(NKILC_TissueSubset, dims = 1:20, reduction = "harmony")
gc()
NKILC_TissueSubset <- FindClusters(object = NKILC_TissueSubset, reduction = "harmony", dims = 1:20,
    resolution = 0.5, print.output = 0, save.SNN = TRUE)


p1 <- DimPlot(NKILC_TissueSubset,reduction = "tsne",group.by = "Donor")
p2 <- DimPlot(NKILC_TissueSubset,reduction = "umap",group.by = "Donor")
p3 <- DimPlot(NKILC_TissueSubset,reduction = "umap",group.by = "Manually_curated_celltype")
p4 <- DimPlot(NKILC_TissueSubset,reduction = "tsne",group.by = "Manually_curated_celltype")
p5 <- DimPlot(NKILC_TissueSubset,reduction = "tsne")
p6 <- DimPlot(NKILC_TissueSubset,reduction = "umap")
g <- grid.arrange(as.list(p1,p2,p3,p4,p5,p6),nrow = 2)
ggsave("NKILC_TissueSubset_DimPlots_20PC.pdf",plot = g, height = 40, width = 60, units = "cm")

p1 <- DotPlot(NKILC_TissueSubset, features = PBMC_NK_markersG,cluster.idents = T) + RotatedAxis()
p2 <- DotPlot(NKILC_TissueSubset, features = PBMC_NK_markersTF,cluster.idents = T) + RotatedAxis()
g <- grid.arrange(as.list(p1,p2),nrow = 1)
ggsave("NKILC_TissueSubset_Dotplots_20PC.pdf",plot = g, height = 30, width = 60, units = "cm")


save(NKILC_TissueSubset, file = "NKILC_TissueSubset_20PC.Robj")



### Top 10 PCs


NKILC_TissueSubset <- RunUMAP(NKILC_TissueSubset, dims = 1:10, reduction = "harmony")
gc()
NKILC_TissueSubset <- RunTSNE(NKILC_TissueSubset, dims = 1:10, reduction = "harmony")
gc()
NKILC_TissueSubset <- FindNeighbors(NKILC_TissueSubset, dims = 1:10, reduction = "harmony")
gc()
NKILC_TissueSubset <- FindClusters(object = NKILC_TissueSubset, reduction = "harmony", dims = 1:10,
    resolution = 0.5, print.output = 0, save.SNN = TRUE)

p1 <- DimPlot(NKILC_TissueSubset,reduction = "tsne",group.by = "Donor")
p2 <- DimPlot(NKILC_TissueSubset,reduction = "umap",group.by = "Donor")
p3 <- DimPlot(NKILC_TissueSubset,reduction = "umap",group.by = "Manually_curated_celltype")
p4 <- DimPlot(NKILC_TissueSubset,reduction = "tsne",group.by = "Manually_curated_celltype")
p5 <- DimPlot(NKILC_TissueSubset,reduction = "tsne")
p6 <- DimPlot(NKILC_TissueSubset,reduction = "umap")
g <- grid.arrange(as.list(p1,p2,p3,p4,p5,p6),nrow = 2)
ggsave("NKILC_TissueSubset_DimPlots_10PC.pdf",plot = g, height = 40, width = 60, units = "cm")

p1 <- DotPlot(NKILC_TissueSubset, features = PBMC_NK_markersG,cluster.idents = T) + RotatedAxis()
p2 <- DotPlot(NKILC_TissueSubset, features = PBMC_NK_markersTF,cluster.idents = T) + RotatedAxis()
g <- grid.arrange(as.list(p1,p2),nrow = 1)
ggsave("NKILC_TissueSubset_Dotplots_10PC.pdf",plot = g, height = 30, width = 60, units = "cm")


save(NKILC_TissueSubset, file = "NKILC_TissueSubset_10PC.Robj")





