library(Matrix)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

result <- jsonlite::fromJSON(file.path("./CTIM_NKILC_TS","exp_data.json"))
total_counts <- result$totalCounts
data_ <- rhdf5::h5read(file.path("./CTIM_NKILC_TS","data.h5"),"X")
dims_ <- rhdf5::h5readAttributes(file = file.path("./CTIM_NKILC_TS","data.h5"), 
                                 name = "X")$shape
expData <- Matrix::sparseMatrix(p = data_$indptr, j = data_$indices, 
                                x = as.numeric(data_$data), 
                                dimnames = list(result$features, result$barcodes),
                                dims=c(dims_[2],dims_[1]), repr="R", index1=FALSE)
plotData <- (jsonlite::fromJSON(file.path("./CTIM_NKILC_TS","plot_data.json"))$data) %>% 
  as.data.frame

dimreduction <- plotData[, paste0("UMAP_",c(1,2))]
seurat_obj = CreateSeuratObject(expData, meta.data = plotData)
head(plotData)
seurat_obj[['umap']] =CreateDimReducObject(embeddings = as.matrix(dimreduction),
			                   key = "UMAP_",
                                           assay = DefaultAssay(seurat_obj))
head(seurat_obj@meta.data)
avg.exp <- AverageExpression(seurat_obj,group.by= "Cluster",return.seurat=T)

markers <- list(
  PBMC_CD8T_markersG = c("CD8A","SELL","CCR7","CCR4","CCR5","CCR9","CCL4","CCL5","CX3CR1","CXCR3","IL2","IL4","IL4R","IL7R","ITGAE","CD69","B3GAT1","KLRB1","KLRF1","KLRK1","KLRC2","GZMK","GZMB","HLA-DRB1","MKI67","TYROBP","MX1","ISG15"),

PBMC_CD8T_markersTF = c("CD8A","TCF7","LEF1","ZNF683","IKZF2","EOMES","TBX21","STAT1"),

PBMC_CD4T_markersG = c("CD4","SELL","CCR7","CCR4","CCR6","CCR10","CXCR3","CXCR5","CCL5","IL2RA","IL4","IL10","KLRB1","KLRF1","GNLY","PRF1","GZMK","GZMB","HLA-DRB1","CTLA4","PDCD1","HAVCR2","LAG3","ISG15"),

PBMC_CD4T_markersTF = c("CD4","TCF7","LEF1","ZNF683","IKZF2","EOMES","TBX21","RORC","STAT1","GATA3","AHR","FOXP3"),

PBMC_NK_markersG = c("IFNG","NCAM1","CCR7","IL2RB","IL7R","IL32","XCL1","GZMK","GZMH","GZMB","KLRB1","KLRC1","KLRC2","KLRG1","PRF1","CXCR1","CX3CR1","CD160","FCGR3A","MKI67"),

PBMC_NK_markersTF = c("IFNG","NCAM1","MYC","TCF7","EOMES","BACH2","JUND","FOSB","TBX21","ZEB2","PRDM1","ZBTB38"),

PBMC_DG_markersG = c("TRDV1","TRDV2","TRGV9","CCR7","IL7R","GZMK","GZMB","CXCR4","CX3CR1","KLRB1","KLRC1","KLRC2","KLRF1","TIGIT"),

PBMC_DG_markersTF = c("TCF7","TBX21","JUN","ZNF638","EOMES")
)

gobs <- lapply(markers, function(i){
  DoHeatmap(avg.exp,features = i,size=3.5,draw.lines = F,angle = 90)+NoLegend()+ scale_fill_gradientn(colors = c("blue", "white", "red"))
})
g <- grid.arrange(grobs = gobs, ncol = 2)
ggsave(filename = "./AVGEXP_NKILC_TS.pdf",plot = g,height = 60,width = 40,units = "cm")


g <- FeaturePlot(seurat_obj,features = c("IFNG","EOMES"),reduction = "umap")
ggsave("./Check_Object_NKILC_TS.pdf",plot = g,height = 30,width = 40,units = "cm")
