library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(SCNPrep)
library(RJSONIO)
library(readr)
library(harmony)

res.fld = "./CTIM_CD4T/"
load(paste(res.fld, "Pan_CD4.Robj", sep = ""))

table(Pan_CD4@meta.data$Manually_curated_celltype)
Pan_CD4[["CellName"]] <- colnames(Pan_CD4)
length(colnames(x=Pan_CD4))
sum(Pan_CD4@assays$originalexp@counts["PTPRC",] == 0 & ((Pan_CD4@assays$originalexp@counts["PF4",]>0) |(Pan_CD4@assays$originalexp@counts["PPBP",]>0) ))
cell = colnames(Pan_CD4)[!(Pan_CD4@assays$originalexp@counts["PTPRC",] == 0 & ((Pan_CD4@assays$originalexp@counts["PF4",]>0) |(Pan_CD4@assays$originalexp@counts["PPBP",]>0) ))]
length(cell)
Pan_CD4 = subset(Pan_CD4, subset = CellName %in% cell)
length(colnames(x=Pan_CD4))

Pan_CD4 <- NormalizeData(object = Pan_CD4, normalization.method = "LogNormalize", scale.factor = 10000)
gc()
Pan_CD4 <- FindVariableFeatures(object = Pan_CD4, selection.method = 'mean.var.plot', mean.cutoff = c(0.01, Inf), dispersion.cutoff = c(0.5, Inf))
gc()
Pan_CD4@assays$originalexp@var.features <-  Pan_CD4@assays$originalexp@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|MAMU-A", Pan_CD4@assays$originalexp@var.features)]
gc()
Pan_CD4 <- ScaleData(object = Pan_CD4, features = VariableFeatures(object = Pan_CD4), vars.to.regress = c("nCount_originalexp", "percent.mito"))
gc()

Pan_CD4 <- RunPCA(object = Pan_CD4,
                features =  VariableFeatures(object = Pan_CD4),
                dims = 1:40)
gc()
Pan_CD4 <- RunHarmony(object = Pan_CD4, group.by.vars = c("Donor"), assay.use = "originalexp", max.iter.harmony = 20)
gc()
Pan_CD4 <- RunUMAP(Pan_CD4, dims = 1:40, reduction = "harmony")
gc()
Pan_CD4 <- RunTSNE(Pan_CD4, dims = 1:40, reduction = "harmony")
gc()
Pan_CD4 <- FindNeighbors(Pan_CD4, dims = 1:40, reduction = "harmony")
gc()

for(res in seq(0.5, 1, 0.1))  {
    Pan_CD4 <- FindClusters(object = Pan_CD4, reduction = "harmony", dims = 1:40,
    resolution = res, print.output = 0, save.SNN = TRUE)
}
gc()



g <- DimPlot(Pan_CD4,reduction="umap",group.by="Donor")
ggsave(file=paste0(res.fld,"Donor_umap.pdf"),plot=g, width=40,height=30,units="cm")
g <- DimPlot(Pan_CD4,reduction="umap")
ggsave(file=paste0(res.fld,"umap.pdf"),plot=g,width=40,height=30,units="cm")

cluster.averages <- AverageExpression(Pan_CD4, return.seurat = FALSE)
write.table(cluster.averages$originalexp, paste0("CD4_avgexp_harmony.tsv"), sep="\t", quote=F, col.names=NA)

save(Pan_CD4,file=paste0(res.fld,"Pan_CD4_processed.Robj"))


migrateSeuratObject(Pan_CD4, 
	            assay="originalexp",
                    species="hs", 
                    outdir = res.fld, 
                    public = T,
		    curated = T, 
                    generateMarkers = T,
		    generateGMTS = F,
    		    name='CTIM_CD4T',
                    token="CTIM_CD4T",
		    link='',
                    description = '')
