library("Seurat")
library("tidyverse")
library("Matrix")
library("readr")
library("MAST")

### Note: clusters should be annotated by this point
load("whole_object.Robj")
outdir = "./" #remember to define output directory
## Below are examples:

MF_DC = subset(whole, idents = c("Macrophage","Xcr1+ cDC","Sirpa+ cDC","pDC"))
FeaturePlot(MF_DC,features = c("Adgre1","Flt3"),reduction = "tsne",split.by = "orig.ident")
DoHeatmap(whole,features = VariableFeatures(whole),group.by = "ident")

### Visualize population ratio changes ###

#### First create a table with the number of cells in each population split by condition ###
Aged_id_table = table(Idents(subset(whole, subset = orig.ident == "Aged")))
Young_id_table = table(Idents(subset(whole, subset = orig.ident == "Young")))
id_table = left_join(as.data.frame(Young_id_table),as.data.frame(Aged_id_table),by=join_by(Var1==Var1)) %>%
  rename("Var1" = "Cell.Type","Freq.x" = "Freq.Young","Freq.y" = "Freq.Aged")

id_table = id_table %>% 
  rowwise %>%
  mutate(log2FC = log2(Freq.Aged/Freq.Young))%>%
  arrange(desc(log2FC))

as.data.frame(id_table)

#### Visualize using ggplot ####

ggplot(data = id_table,aes(x=reorder(Cell.Type, -log2FC),y=log2FC,fill = Cell.Type))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Cell Types")+
  ylab("log2FC in cell number (Aged/Young)")

#### save plot using ggsave ####

g <- ggplot(data = id_table,aes(x=reorder(Cell.Type, -log2FC),y=log2FC,fill = Cell.Type))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Cell Types")+
  ylab("log2FC in cell number (Aged/Young)")

ggsave(paste0(outdir,"population_changes.pdf"),plot = g,width = 40,height = 30,units = "cm")

