
library(Seurat)
library(tidyverse)
library(harmony)
library(gridExtra)

path = "./SCEPrep/"
load(paste0(path,"CD8.Robj"))

###### Checking info ###########
table(CD8@meta.data$Organ)
table(CD8@meta.data$Manually_curated_celltype)
print(colnames(CD8@meta.data))
table(CD8@meta.data$Sex)
table(CD8@meta.data$Age_range)

###############################

CD8@meta.data$Chemistry <- NULL
CD8@meta.data$Majority_voting_CellTypist <- NULL
CD8@meta.data$Predicted_labels_CellTypist <- NULL

CD8@reductions$X_umap <- NULL
CD8@reductions$umap <- NULL
CD8@reductions$pca <- NULL

#remove RBC and platelets
CD8[["CellName"]] <- colnames(CD8)
sum(CD8@assays$originalexp@counts["PTPRC",] == 0 & ((CD8@assays$originalexp@counts["PF4",]>0) |(CD8@assays$originalexp@counts["PPBP",]>0) ))
cell = colnames(CD8)[!(CD8@assays$originalexp@counts["PTPRC",] == 0 & ((CD8@assays$originalexp@counts["PF4",]>0) |(CD8@assays$originalexp@counts["PPBP",]>0) ))]
length(cell)
CD8 = subset(CD8, subset = CellName %in% cell)
length(colnames(x=CD8))

#remove doublets based on gene expression profiles
sum(CD8@assays$originalexp@counts["CD4",] > 0 & ((CD8@assays$originalexp@counts["CD8A",]>0) |(CD8@assays$originalexp@counts["CD8B",]>0)))
cell = colnames(CD8)[!(CD8@assays$originalexp@counts["CD4",] > 0 & ((CD8@assays$originalexp@counts["CD8A",]>0) |(CD8@assays$originalexp@counts["CD8B",]>0)))]
length(cell)
CD8 = subset(CD8, subset = CellName %in% cell)
length(colnames(x=CD8))

sum(CD8@assays$originalexp@counts["CD19",]> 0 & ((CD8@assays$originalexp@counts["CD3E",]>0) |(CD8@assays$originalexp@counts["CD3D",]>0) |(CD8@assays$originalexp@counts["CD3G",]>0)  ))
cell = colnames(CD8)[!(CD8@assays$originalexp@counts["CD19",]> 0 & ((CD8@assays$originalexp@counts["CD3E",]>0) |(CD8@assays$originalexp@counts["CD3D",]>0) |(CD8@assays$originalexp@counts["CD3G",]>0)  ))]
length(cell)
CD8 = subset(CD8, subset = CellName %in% cell)
length(colnames(x=CD8))

CD8@meta.data$CellName <- NULL
gc()

p <- ggplot(CD8@meta.data)+
  geom_histogram(aes(x=percent.mito,y = after_stat(width*density*100)),bins = 100)+
  facet_wrap(vars(Donor),ncol=3)+
  geom_vline(xintercept = 0.05,col="red") +
  ylab("proportion")+
  xlab("fraction of mitochondrial RNA")
ggsave(paste0(path,"mt_percent_byDonor.pdf"),plot = p)


##################################################################################
CD8 <- NormalizeData(object = CD8, normalization.method = "LogNormalize", scale.factor = 10000)
gc()
CD8 <- FindVariableFeatures(object = CD8, selection.method = 'mean.var.plot', mean.cutoff = c(0.01, Inf), dispersion.cutoff = c(0.55, Inf))
gc()
CD8@assays$originalexp@var.features <-  CD8@assays$originalexp@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|MAMU-A", CD8@assays$originalexp@var.features)]
gc()
CD8 <- ScaleData(object = CD8, features = VariableFeatures(object = CD8), vars.to.regress = c("nCount_originalexp", "percent.mito"))
gc()

CD8 <- RunPCA(object = CD8,
              features =  VariableFeatures(object = CD8),
              dims = 1:30)
gc()
CD8 <- RunHarmony(object = CD8, group.by.vars = c("Donor"), assay.use = "originalexp", max.iter.harmony = 20)
gc()
CD8 <- RunUMAP(CD8, dims = 1:30, reduction = "harmony")
gc()
CD8 <- RunTSNE(CD8, dims = 1:30, reduction = "harmony")
gc()
CD8 <- FindNeighbors(CD8, dims = 1:30, reduction = "harmony", k.param = 15)
gc()
CD8 <- FindClusters(CD8, resolution = 1.0)
gc()

gobs <- list(
  p1 = DimPlot(CD8,reduction = "pca",group.by = "Manually_curated_celltype"),
  p2 = DimPlot(CD8,reduction = "tsne",group.by = "Manually_curated_celltype"),
  p3 = DimPlot(CD8,reduction = "umap",group.by = "Manually_curated_celltype"),
  p4 = DimPlot(CD8,reduction = "pca",group.by = "Donor"),
  p5 = DimPlot(CD8,reduction = "tsne",group.by = "Donor"),
  p6 = DimPlot(CD8,reduction = "umap",group.by = "Donor")
)

g <- grid.arrange(grobs = gobs, ncol = 3)
ggsave(paste0(path,"dimplots.pdf"),plot=g)

cluster.averages <- AverageExpression(CD8, return.seurat = FALSE)
write.table(cluster.averages$originalexp, paste0(path,"CD8_avgexp_harmony.tsv"), sep="\t", quote=F, col.names=NA)

save(CD8,file = "CD8_reprocessed.Robj")
