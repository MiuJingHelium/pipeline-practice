#library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
#library(SCNPrep)
#library(RJSONIO)
#library(readr)

res.fld = "../data/"
load(paste(res.fld, "T_and_ILC.Robj", sep = ""))

Idents(whole) <- whole@meta.data$Manually_curated_celltype
whole@meta.data$log2nFeature <- log2(whole@meta.data$nFeature_originalexp)
whole@meta.data$log2nUMI <- log2(whole@meta.data$nCount_originalexp)

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

#CD3 <- subset(whole,subset = (CD3D > 0) | (CD3E > 0) | (CD3G > 0))
#CD4 <- subset(CD3, subset = CD4 > 0)
#save(CD3, file = "./SCEPrep/CD3.Robj")
#save(CD4,file="./SCEPrep/CD4.Robj")
#rm(CD4)
#gc()
#CD8 <- subset(CD3, subset = (CD8A > 0) | (CD8B > 0))
#save(CD8, file = "./SCEPrep/CD8.Robj")
#rm(CD8)
#gc()
#CD3_DN <- subset(CD3,subset = (CD4 == 0) & (CD8A == 0) & (CD8B == 0))
#save(CD3_DN, file = "./SCEPrep/CD3_DN.Robj")


Pan_CD4 <- subset(whole,idents = c("Tnaive/CM_CD4_activated","Tnaive/CM_CD4","Tregs","Tfh","Teffector/EM_CD4","Trm_Th1/Th17"))
save(Pan_CD4, file="./Pan_CD4.Robj")
Pan_CD8 <- subset(whole, idents = c("Tem/emra_CD8","Tnaive/CM_CD8","Trm/em_CD8","MAIT","Tgd_CRTAM+","Trm_Tgd","Trm_gut_CD8"))
save(Pan_CD8,file="./Pan_CD8.Robj")
NK_ILC <- subset(whole, idents = c("ILC3","NK_CD16+","NK_CD56bright_CD16-"))
save(NK_ILC,file="./NK_ILC.Robj")



#NK <- subset(whole,subset = (NCAM1 > 0 | FCGR3A > 0) & TYROBP > 0 )
#save(NK, file="./SCEPrep/NK.Robj")

#migrateSeuratObject(whole, res.fld, additionalAnnotations)
#migrateMarkers(paste(res.fld, "markers.tsv", sep = ""), c("CellTypes"), res.fld)
#migrateSeuratObject(whole, 
#	            assay="originalexp",
#                    species="hs", 
#                    outdir = res.fld, 
#                    public = T,
#		    curated = F,
#		    markers = read.table("./SCEPrep/markers.tsv",sep = "\t",header = T), 
#                    generateMarkers = F,
#		    generateGMTS = F,
#    		    name='CTIM-T',
#                    token="CTIM-T",
#		    link='',
#                    description = '')
