library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(SCNPrep)
library(RJSONIO)
library(readr)

res.fld = "./CTIM_NKILC/"
load(paste(res.fld, "NK_ILC_processed.Robj", sep = ""))

#Idents(whole) <- whole@meta.data$Manually_curated_celltype
#whole@meta.data$log2nFeature <- log2(whole@meta.data$nFeature_originalexp)
#whole@meta.data$log2nUMI <- log2(whole@meta.data$nCount_originalexp)

#find markers based on their cell type
#barcodes <- rownames(whole@meta.data)
#ann <- data.frame(whole@meta.data) %>%
#  mutate(Tissue_type = dplyr::case_when(
#         colData(whole)$Organ %in% lymphoid_tissue ~ "LYM",
#         .default = "Non-LYM"
#          )) %>%
#  dplyr::select(Tissue_type)
#rownames(ann) <- barcodes

## factorize each annotation. You can do custom factors
#for (column in colnames(ann)) {
#  ann[, column] <- as.factor(ann[, column])
#}
#additionalAnnotations <- list(ann)

#SPL <- subset(whole,Organ=="SPL")
#save(SPL,file="./SCEPrep/SPL.Robj")

#migrateSeuratObject(whole, res.fld, additionalAnnotations)
#migrateMarkers(paste(res.fld, "markers.tsv", sep = ""), c("CellTypes"), res.fld)
g <- DimPlot(NK_ILC,reduction="umap",group.by="Donor")
ggsave(file=paste0(res.fld,"Donor_umap.pdf"),plot=g, width=40,height=30,units="cm")
g <- DimPlot(NK_ILC,reduction="umap")
ggsave(file=paste0(res.fld,"umap.pdf"),plot=g,width=40,height=30,units="cm")



migrateSeuratObject(NK_ILC, 
	            assay="originalexp",
                    species="hs", 
                    outdir = res.fld, 
                    public = T,
		    curated = T, 
                    generateMarkers = T,
		    generateGMTS = F,
    		    name='CTIM_NKILC',
                    token="CTIM_NKILC",
		    link='',
                    description = '')
