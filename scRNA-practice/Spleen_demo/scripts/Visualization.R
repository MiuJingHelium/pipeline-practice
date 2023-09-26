library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)#for plotting
library(gridExtra)#for plotting

############################  Notes   #3#################################
# More basic visualization via Seurat are made along the analysis process.
# This script contains code for making more customized plots as well as 
# those requiring additional processed data.
#
# Or this might be simply a code example collection and not written intended
# for running as a typical script

###################### Visualization of DE data ####################


path = "../manualQC/analysis/Major_cell_types/raw_DE_tables/" #path may be incorrect
files = list.files(path = path,full.names = T)




samples <- c("B cell","CD4 T cell","CD8 T cell","Macrophage")
n = length(files)
plts = lapply(1:n, function(i) {
  data = read.table(file = files[i],row.names = NULL,sep = "\t",header = T)
  plts[[n]] <- ggplot(data = test) +
    geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj),col = threshold))+
    geom_vline(xintercept = -0.5,col="red") + 
    geom_vline(xintercept = 0.5,col="red") +
    geom_hline(yintercept = 1.5,col="red") +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))+
    xlab("Average Log2 Fold Change (Young/Aged)")+
    ylab("-log10(adjusted P-value)")+
    labs(title = samples[i])
}
)

grid.arrange(grobs = plts,ncol=2)
#another version of placing plots together is available in QC_preprocessing_Single.R

##########################################################

############## Population ratio plot #####################

path = "../manualQC/analysis/"
load(paste0(path,"whole_object.Robj"))


Aged_id_table = table(Idents(subset(whole, subset = orig.ident == "Aged")))
Aged_id_table



Young_id_table = table(Idents(subset(whole, subset = orig.ident == "Young")))
Young_id_table


id_table = left_join(as.data.frame(Young_id_table),as.data.frame(Aged_id_table),by=join_by(Var1==Var1)) %>%
  rename("Var1" = "Cell.Type","Freq.x" = "Freq.Young","Freq.y" = "Freq.Aged")

id_table = id_table %>% 
  rowwise %>%
  mutate(log2FC = log2(Freq.Aged/Freq.Young))%>%
  arrange(desc(log2FC))



as.data.frame(id_table)

ggplot(data = id_table,aes(x=reorder(Cell.Type, -log2FC),y=log2FC,fill = Cell.Type))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Cell Types")+
  ylab("log2FC in cell number (Aged/Young)")


##############################################################################

