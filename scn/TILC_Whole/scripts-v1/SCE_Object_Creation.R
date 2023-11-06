library(SingleCellExperiment)
library(zellkonverter)
library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(data.table)
library(harmony)
library(Matrix)
library(SCNPrep)
library(RJSONIO)
library(readr)

options(scipen=999)
getPalette.1 <- colorRampPalette(brewer.pal(9, "Set1"))

args <- commandArgs(trailingOnly = T)
indir <- args[1]
res.fld <- args[2]

sce <- readH5AD(indir)#file needs to be the raw one
lymphoid_tissue <- c("THY","LLN","SPL","MLN","BWA")

Tissue_type <- case_when(
              colData(sce)$Organ %in% lymphoid_tissue ~ "LYN",
              .default = "Non-LYN"
)

colData(sce)$Tissue_type <- Tissue_type

whole <- as.Seurat(sce, counts = "counts",data = "X",slot = "originalexp",project="AllTissueTCells")
rm(sce)
gc()

mito.genes <- c(grep("^MT-|^mt-", rownames(x = whole), value = T))
percent.mito <- Matrix::colSums(x = GetAssayData(object = whole, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = whole, slot = 'counts'))
whole[['percent.mito']] <- percent.mito

whole@meta.data$percent.mito <- percent.mito
xx <- whole@meta.data %>%
  group_by(Donor) %>%
  mutate(percent.mito = round(percent.mito * 100, 2)) %>%
  summarise(mean_mt = mean(percent.mito),
            sd_mt = sd(percent.mito),
            median_mt = median(percent.mito),
            mean_plus_sd = mean_mt + sd_mt,
            median_plus_sd = median_mt + sd_mt)


p <- ggplot(whole@meta.data, aes(x = percent.mito)) +
  geom_histogram(bins = 100, aes(y = stat(width*density*100))) +
  theme_classic() +
  facet_wrap(~Donor, ncol = 5) +
  geom_vline(xintercept = 0.1) + 
  ggtitle("mitochondrial content")
ggsave(file=paste0(res.fld,"percent_mito_plots_byDonor.pdf"), plot = p,width = 40, height = 30, units = "cm")

plts <- lapply(c(0.05,0.1,0.15,0.2,0.3,0.5),function(i){
  FeaturePlot(whole, features = "percent.mito",min.cutoff = i) +
    ggtitle(paste0("min.cutoff = ",i))
})
g <- arrangeGrob(grobs = plts,ncol=3)
ggsave(file=paste0(res.fld,"percent_mito_plots_byCutoff.pdf"), plot = g,width = 40, height = 30, units = "cm")

rm(g)
rm(plts)
gc()

##############
whole@meta.data$log2nFeature <- log2(whole@meta.data$nFeature_originalexp)
whole@meta.data$log2nUMI <- log2(whole@meta.data$nCount_originalexp)



whole <- NormalizeData(object = whole, normalization.method = "LogNormalize", scale.factor = 10000)
whole <- FindVariableFeatures(object = whole, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
## in variable genes you can remove genes that can introduce batch effect,
## such as TR-IG genes for humans
whole@assays$originalexp@var.features <-  whole@assays$originalexp@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|MAMU-A", whole@assays$originalexp@var.features)]
whole <- ScaleData(object = whole, features = VariableFeatures(object = whole), vars.to.regress = c("nCount_originalexp", "percent.mito"))
gc()
whole <- RunPCA(object = whole, verbose=FALSE)
g <- DimPlot(whole,reduction="pca",group.by="Donor")
ggsave(file=paste0(res.fld,"pca.pdf"),plot=g,width = 40, height = 30,units="cm")
## run harmony
## label samples that come from the same subjects (depends on how you named
## your samples and how many samples subject has)
#whole@meta.data$Subject <- whole@meta.data$Donor
#whole <- RunHarmony(whole, "Subject")
whole[['umap']] = CreateDimReducObject(embeddings = as.matrix(whole@reductions$X_umap@cell.embeddings), 
                                            key = "UMAP_", 
                                            assay = DefaultAssay(whole))
Idents(whole) <- whole@meta.data$Manually_curated_celltype

## do UMAP and clustering
#whole <- whole %>% 
#  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
#  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
#  FindClusters(resolution = 0.5) %>% 
#  identity()
save(whole, file = paste(res.fld, "whole_object.Robj", sep = ""))

#find markers based on their cell type
#whole.markers <- FindAllMarkers(object = whole,
#                               only.pos = TRUE,
#                               min.pct = 0.15,
#                               thresh.use = 0.15)
#write.table(whole.markers, paste(res.fld, "markers.tsv", sep = ""), sep="\t", quote=F, row.names=F)

#barcodes <- rownames(whole@meta.data)
#ann <- data.frame(whole@meta.data) %>%
#  mutate(Tissue_type = dplyr::case_when(
#         colData(sce)$Organ %in% lymphoid_tissue ~ "LYM",
#         .default = "Non-LYM"
#          )) %>%
#  dplyr::select(Tissue_type)
#rownames(ann) <- barcodes

## factorize each annotation. You can do custom factors
#for (column in colnames(ann)) {
#  ann[, column] <- as.factor(ann[, column])
#}
#additionalAnnotations <- list(ann)



#migrateSeuratObject(whole, res.fld, additionalAnnotations)
#migrateMarkers(paste(res.fld, "markers.tsv", sep = ""), c("CellTypes"), res.fld)
#migrateSeuratObject(object, 
#                    species="hs", 
#                    outdir = res.fld, 
#                    public = T, 
#                    generateMarkers = T,
#                    token="CTIM-T")
