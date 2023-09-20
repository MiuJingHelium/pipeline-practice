library("Seurat")
library("tidyverse")
library("Matrix")
library("readr")
library("MAST")

args = commandArgs(trailingOnly = T)
data_path = args[1]
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

test_clusters <- c("B cell","CD4 T cell","CD8 T cell","Macrophage")
#apply(conditions_comparisons, 2, function(current_comparison) {
for (i in all_clusters) {
  output_file <- paste0(data_path,'DE/', paste0('Spleen_Young_vs_Spleen_Aged_', i, '.tsv'))
  #output_file_flt <- paste0(data_path,'DE/', paste0('Spleen_Young_vs_Spleen_Aged_', i, '_filtered.tsv'))
  de_results_tbl <- lapply(i, de_clusters) %>% bind_rows()
  #de_results_tbl <- de_results_tbl[order(abs(de_results_tbl$avg_log2FC),decreasing=T),]
  write_tsv(de_results_tbl, output_file)
  #de_results_tbl <- de_results_tbl[de_results_tbl$p_val_adj < 0.05,]
  #de_results_tbl <- de_results_tbl[abs(de_results_tbl$avg_log2FC) > 0.5]
  #write_tsv(de_results_tbl, output_file_flt)
  gc()
}



