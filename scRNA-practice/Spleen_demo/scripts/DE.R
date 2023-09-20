library("Seurat")
library("tidyverse")
library("Matrix")
library("readr")
library("MAST")
library("gridExtra")
args = commandArgs(trailingOnly = T)
data_path = args[1]
padj_thres = 0.05
log2FC_thres = 1 #consider adapt into cmd-line args
# should end with "/"



##### The function for within cluster DE #####

de_clusters <- function(current_cluster) {
  groupA <- dplyr::filter(data, Cluster == current_cluster & Sample %in% c("Young"))$CellName 
  #May need to update searching criteria based on sample name
  groupB <- dplyr::filter(data, Cluster == current_cluster & Sample %in% c("Aged"))$CellName
  if (length(groupA) >= 3 && length(groupB) >= 3) { #Needs to have minimal of 3 cells
    both_markers <- FindMarkers(whole, ident.1=groupA, ident.2=groupB, test.use="MAST", min.pct=0.00, logfc.threshold = 0)
    marker_genes <- row.names(both_markers)
    mutate(both_markers, gene = marker_genes, cluster = current_cluster)
  } else {
    print(paste0('Cluster ', current_cluster, ': not enough cells'))
    tibble(p_val = double(), avg_logFC = double(), pct.1 = double(), pct.2 = double(),
           p_val_adj = double(), cluster = character(), gene = character())
  }
}

##############################################

data = read_tsv(paste0(data_path,"data_for_plot.tsv"))
#problem with the original file format; column name for cell names is missing
#CellName column was inserted and row.names=F instead.
load(paste0(data_path,"whole_object.Robj"))
all_clusters <- names(table(data$Cluster)) #Prep list of clusters
print(all_clusters)

test_clusters <- c("CD4 T cell","CD8 T cell","B cell","Macrophage")
#apply(conditions_comparisons, 2, function(current_comparison) {
for (i in test_clusters) {
  output_file <- paste0(data_path,'DE/', paste0('Spleen_Young_vs_Spleen_Aged_', i, '.tsv'))
  output_file_flt <- paste0(data_path,'DE/', paste0('Spleen_Young_vs_Spleen_Aged_', i, '_filtered.tsv'))
  de_results_tbl <- lapply(i, de_clusters) %>% bind_rows() %>%
			arrange(desc(abs(avg_log2FC))) %>%
			mutate(threshold = case_when(
			(p_val_adj < padj_thres) & (abs(avg_log2FC) > log2FC_thres) ~ "Pass",
			.default = "Fail"
			))
  #make volcano plot
  #jpeg(paste0(data_path,'DE/',paste0('Spleen_Young_vs_Spleen_Aged_', i,'Volcano', '.jpeg')))
  #ggplot(data = de_results_tbl) +
#	geom_point(aes(x = avg_log2FC, y= -log10(p_val_adj),col = threshold))+
#	geom_vline(xintercept = -0.5, col="red")+
#	geom_vline(xintercept = 0.5, col="red")+
#	geom_hline(yintercept = 1.5, col="red")+
#	scale_color_manual(values = c("#00AFBB","#E7B800"))+
#	xlab("Average Log2 Fold Change (Young/Aged)")+
#	ylab("-log10(adjusted P-value)")
# dev.off()
  write.table(de_results_tbl, output_file,sep="\t",row.names=T)
  de_results_tbl <- de_results_tbl %>% 
			filter(p_val_adj < padj_thres) %>%
			filter(abs(avg_log2FC) > log2FC_thres)
			
  write.table(de_results_tbl, output_file_flt,sep="\t",row.names=T)
  gc()
}

#Plot volcano plot
# Path needs to contain one type of table.i.e. entirely unfiltered
#path = paste0(data_path,"DE/")
#files = list.files(path = path,full.names = T)
#n = length(files)
#samples = unlist(map(files,function(i){
#	unlist(strsplit(unlist(strsplit(unlist(strsplit(i,
#	split = "/"))[length(unlist(strsplit(i,split = "/")))], 
#	split = "_"))[6],
#	split = "[.]"))[1]
#	}))

#plts = lapply(1:n, function(i) {
#  datai = read.table(file = files[i],row.names = NULL,sep = "\t",header = T)
#  ggplot(data = datai) +
#  geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj),col = threshold))+
#  geom_vline(xintercept = -(log2FC_thres),col="red") + 
#  geom_vline(xintercept = log2FC_thres,col="red") +
#  geom_hline(yintercept = -log10(padj_thres),col="red") +
#  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
#  xlab("Average Log2 Fold Change (Young/Aged)")+
#  ylab("-log10(adjusted P-value)")+
#  labs(title = samples[i])+
#  theme(legend.position = "none")
#  }
#)

#jpeg(paste0(data_path,"DE/Volcano.jpeg"))
#grid.arrange(grobs = plts,ncol=2)
#dev.off()
