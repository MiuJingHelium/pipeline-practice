library(Seurat)
library(tidyverse)
library(Matrix)

### The script can be hard to automate without knowing the actual names
### Can try to merge as a part of the the loop? 
indir = "../outs/"
outdir="./"
condition_path = "../"
samples = list.files(path=indir,full.names = F)
n = length(samples)
#objs = list(n)
for (i in (1:n)){
  load(paste0(indir,samples[i],"/",samples[i],"_filtered.Robj"))
  #objs[[i]] = assign(samples[i],obj_filtered)#assume using the same set of script
  assign(samples[i],obj_filtered)
}
rm(obj_filtered)
gc()


######################### NOT GENERALIZED ###############################

whole = merge(eval(as.symbol(samples[1])),y=list(eval(as.symbol(samples[2]))),add.cell.ids = samples)

rm(Spleen_Aged)
rm(Spleen_Young)

gc()

#Add orig.ident based on range

aged_range = range(grep("Aged",rownames(whole@meta.data)))
young_range = range(grep("Young",rownames(whole@meta.data)))
whole@meta.data$orig.ident[aged_range[1]:aged_range[2]] = "Aged"
whole@meta.data$orig.ident[young_range[1]:young_range[2]] = "Young"

######################### NOT GENERALIZED ###############################


######################### Basic plotting and more filtering #####################

### Use threshold value for vline plotting

jpeg(paste0(outdir,"mt_percent_histogram.jpeg"))
ggplot(whole@meta.data)+
  geom_histogram(aes(x=percent.mt,y = after_stat(width*density*100)),bins = 100)+
  facet_wrap(vars(orig.ident),ncol=1)+
  geom_vline(xintercept = 5,col="red") +
  ylab("proportion")+
  xlab("percentage of mitochondrial RNA")
dev.off()

jpeg(paste0(outdir,"log2nFeature_RNA_histogram_before.jpeg"))
ggplot(whole@meta.data, aes(x = log2(nFeature_RNA))) +
  geom_histogram(bins = 100, aes(y = after_stat(width*density*100))) +
  theme_classic() +
  facet_wrap(~orig.ident, ncol = 1) +
  geom_vline(xintercept = 8.5,col="red") + 
  geom_vline(xintercept = 10.5,col="red") + 
  ggtitle("nGene") +
  ylab("")
dev.off()

cells <- length(colnames(x = whole))
print(cells)

whole = subset(whole, subset = nFeature_RNA_log2 <= 12 ) #mt_percentage have been already filtered
cells <- length(colnames(x = whole))

print(cells)
jpeg(paste0(outdir,"log2nUMI_histogram_before.jpeg"))
ggplot(whole@meta.data, aes(x = log2(nCount_RNA))) +
  geom_histogram(bins = 100, aes(y = after_stat(width*density*100))) +
  theme_classic() +
  facet_wrap(~orig.ident, ncol = 1) +
  geom_vline(xintercept = 9,col="red") + 
  geom_vline(xintercept = 15,col="red") + 
  ggtitle("nUMI") +
  ylab("")
dev.off()

whole@meta.data$log2nUMI = log2(whole@meta.data$nCount_RNA)
whole = subset(whole, subset = log2nUMI > 9)
cells <- length(colnames(x = whole))
print(cells)

#remove platelets
sum(whole[["RNA"]]@counts["Ppbp",]>0 | whole[["RNA"]]@counts["Pf4",]>0)
cell = colnames(whole)[!(whole[["RNA"]]@counts["Ppbp",]>0 | whole[["RNA"]]@counts["Pf4",]>0)]
length(cell)
whole = subset(whole, subset = CellName %in% cell)
length(colnames(x=whole))

###### Expression signature based doublet removal ########
#operate over conditions:
#Each condition involves a set of logical conditions (min = 2)
#Condition would be a list of sets of genes that can't co-express
#Condition will generate boolean vectors
#ORing the conditions will return a boolean vector for final cell selection

#Input would be a table of some rows. Each row is a gene set.

whole[["CellName"]] <- colnames(whole)
conditions <- read.table(paste0(condition_path,"conditions.tsv"))
n = nrow(conditions)

bools <- map(1:n,function(i){
  whole[["RNA"]]@counts[conditions[i,1],]>0 & whole[["RNA"]]@counts[conditions[i,2],]>0
}) #return a list with boolean vectors with the same length as the current number cells


summary <- map(1:n,function(i){
  sum(whole[["RNA"]]@counts[conditions[i,1],]>0 & whole[["RNA"]]@counts[conditions[i,2],]>0)
}) 

summary <- unlist(summary)
conditions = cbind(conditions,summary)
write.table(conditions,paste0(condition_path,"doublet_numbers.tsv"),sep = "\t",row.names = F,col.names = F)

cell_bool <- Reduce("|",bools)
cells <- colnames(whole)[!cell_bool]
whole = subset(whole, subset = CellName %in% cells)
length(colnames(x=whole))

##########################################################

################# START ANALYSIS #########################

whole <- NormalizeData(object = whole, normalization.method = "LogNormalize", scale.factor = 10000)
whole <- FindVariableFeatures(object = whole, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
whole@assays$RNA@var.features <-  whole@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|MAMU-A", whole@assays$RNA@var.features)]
whole <- ScaleData(object = whole, features = VariableFeatures(object = whole), vars.to.regress = c("nCount_RNA", "percent.mito"))
gc()
whole <- RunPCA(object = whole,
                features =  VariableFeatures(object = whole),
                dims = 1:30)

whole <- RunTSNE(object = whole, dims = 1:30)
## CLUSTERING
whole <- FindNeighbors(object = whole, dims = 1:30)
whole <- FindClusters(object = whole, resolution = 1.0)

expr.avg <- AverageExpression(whole, return.seurat = T)#Seurat object with the avg exp

#########################################################

############## EXPORT FILES #############################


save(whole, file = paste0(outdir, "whole_object.Robj"))

expData <- GetAssayData(object = whole, slot = 'data')
save(expData, file=paste0(outdir,"expData.Rda"))
save(expr.avg, paste0(outdir, "avgExp.Robj"))

dataForPlot <- as.data.frame(whole@reductions$tsne@cell.embeddings)
dataForPlot$CellName <- wholewhole@meta.data$CellName
dataForPlot$Sample <- whole@meta.data$orig.ident
dataForPlot$Cluster <-  Idents(object = whole)
dataForPlot$nUmi <- whole@meta.data$nCount_RNA
dataForPlot$nGene <- whole@meta.data$nFeature_RNA
dataForPlot$nUmiLog2 <- log2(whole@meta.data$nCount_RNA)
dataForPlot$nGeneLog2 <- log2(whole@meta.data$nFeature_RNA)

write.table(dataForPlot, paste0(outdir,"data_for_plot.tsv"), sep="\t", quote=F,row.names = F)

## FINDING ANS SAVING MARKERS
whole.markers <- FindAllMarkers(object = whole,
                                only.pos = TRUE,
                                min.pct = 0.10,
                                thresh.use = 0.10)
write.table(whole.markers, "markers.tsv", sep="\t", quote=F, row.names=F)

###############################################################









