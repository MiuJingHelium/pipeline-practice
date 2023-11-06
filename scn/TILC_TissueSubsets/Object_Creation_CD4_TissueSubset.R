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

load("CD4_TissueSubset_30PC.Robj")

table(CD4_TissueSubset@meta.data$Manually_curated_celltype)
length(colnames(x=CD4_TissueSubset))

CD4_TissueSubset@meta.data <- CD4_TissueSubset@meta.data %>%
  mutate_if(is.character, as.factor)

# Add markers as a list with names as resolutions
# Testing codes ....

markers <- lapply(seq(0.1,1,0.1),function(x){
  #Note: you can directly use the column name and skip retrieving the variable
  Idents(CD4_TissueSubset) <- paste0("originalexp","_snn_res.",x)
  #lapply cannot assign names inside the function, but you can either pass a named input or name the output outside lapply
  marker <- FindAllMarkers(CD4_TissueSubset, only.pos=T)
  marker$cluster <- as.character(marker$cluster)
  marker
})
names(markers) <- paste0("res.",seq(0.1,1,0.1))
res.fld = "CTIM_CD4T_TS"

migrateSeuratObject(CD4_TissueSubset, 
	            assay="originalexp",
                    species="hs", 
                    outdir = res.fld, 
                    public = T,
		    curated = T,
		    slot = 'counts',
		    marker = markers, 
                    generateMarkers = F,
		    generateGMTS = F,
    		    name='CTIM_CD4T_TS',
                    token="CTIM_CD4T_TS_uIJNVqweuk",
		    link='https://www.tissueimmunecellatlas.org/',
                    description = 'LLN,MLN,LNG,SPL,and BMA subsets of CTIM_CD4T;top 30 PCs used.')
file.rename(from = file.path(res.fld,"dataset.json"),to = file.path(res.fld,"zzzz.json"))
