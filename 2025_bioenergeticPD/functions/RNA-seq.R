library(Seurat) 
library(ggplot2) 
library(sctransform)

library(readr)
count<-read_tsv('~/drive/bioenergeticsPD/fromXi/data/RNA-seq/GSE157783_Luxembourg/RAW/GSE157783_IPDCO_hg_midbrain_UMI/IPDCO_hg_midbrain_UMI.tsv',col_names='true')
