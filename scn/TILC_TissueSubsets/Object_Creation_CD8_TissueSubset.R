library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(SCNPrep)
library(RJSONIO)
library(readr)
library(harmony)
library(stringi)
#library(argparse)

#res.fld = "./CTIM_CD8T/"
#load(paste(res.fld, "Pan_CD8.Robj", sep = ""))
load("CD8_TissueSubset_20PC.Robj")

table(CD8_TissueSubset@meta.data$Manually_curated_celltype)
length(colnames(x=CD8_TissueSubset))
#sum(Pan_CD8@assays$originalexp@counts["PTPRC",] == 0 & ((Pan_CD8@assays$originalexp@counts["PF4",]>0) |(Pan_CD8@assays$originalexp@counts["PPBP",]>0) ))
#cell = colnames(Pan_CD8)[!(Pan_CD8@assays$originalexp@counts["PTPRC",] == 0 & ((Pan_CD8@assays$originalexp@counts["PF4",]>0) |(Pan_CD8@assays$originalexp@counts["PPBP",]>0) ))]
#length(cell)
#Pan_CD8 = subset(Pan_CD8, subset = CellName %in% cell)

CD8_TissueSubset@meta.data <- CD8_TissueSubset@meta.data %>%
  mutate_if(is.character, as.factor)

# Add markers as a list with names as resolutions
# Testing codes ....

markers <- lapply(seq(0.1,1,0.1),function(x){
  #Note: you can directly use the column name and skip retrieving the variable
  Idents(CD8_TissueSubset) <- paste0("originalexp","_snn_res.",x)
  #lapply cannot assign names inside the function, but you can either pass a named input or name the output outside lapply
  marker <- FindAllMarkers(CD8_TissueSubset, only.pos=T)
  marker$cluster <- as.character(marker$cluster)
  marker
})
names(markers) <- paste0("res.",seq(0.1,1,0.1))
res.fld = "CTIM_CD8T_TS"

migrateSeuratObject(CD8_TissueSubset, 
	            assay="originalexp",
                    species="hs", 
                    outdir = res.fld, 
                    public = T,
		    curated = T,
		    slot = 'counts',
		    marker = markers, 
                    generateMarkers = F,
		    generateGMTS = F,
    		    name='CTIM_CD8T_TS',
                    token="CTIM_CD8T_TS_eBvfghP",
		    link='https://www.tissueimmunecellatlas.org/',
                    description = 'LLN,MLN,LNG,SPL,and BMA subsets of CTIM_CD8T;top 20 PCs used.')
file.rename(from = file.path(res.fld,"dataset.json"),to = file.path(res.fld,"zzzz.json"))
