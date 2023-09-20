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

path = paste0(data_path,"DE/")
files = list.files(path = path,full.names = T)
n = length(files)
samples = unlist(map(files,function(i){
        unlist(strsplit(unlist(strsplit(unlist(strsplit(i,
        split = "/"))[length(unlist(strsplit(i,split = "/")))],
        split = "_"))[6],
        split = "[.]"))[1]
        }))

plts = lapply(1:n, function(i) {
  datai = read.table(file = files[i],row.names = NULL,sep = "\t",header = T)
  ggplot(data = datai) +
  geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj),col = threshold))+
  geom_vline(xintercept = -(log2FC_thres),col="red") +
  geom_vline(xintercept = log2FC_thres,col="red") +
  geom_hline(yintercept = -log10(padj_thres),col="red") +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  xlab("Average Log2 Fold Change (Young/Aged)")+
  ylab("-log10(adjusted P-value)")+
  labs(title = samples[i])+
  theme(legend.position = "none")
  }
)

jpeg(paste0(data_path,"DE/Volcano.jpeg"))
grid.arrange(grobs = plts,ncol=2)
dev.off()

