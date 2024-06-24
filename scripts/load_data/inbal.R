
source("scripts/utils.R")
suppressPackageStartupMessages({
  library(RcppML)
  library(tidyverse)
  library(cowplot) 
})

data_path=paste0(base_path,"inbal/data/")
figure_path=paste0(base_path,"inbal/figures/")

mat = readMM(paste0(data_path,"GSE199317_ONC-retina.mtx"))
colnames(mat) = fread(paste0(data_path,"GSE199317_ONC-retina_cell.tsv"),sep="\t",header=F)$V1
row.names(mat) = fread(paste0(data_path,"GSE199317_ONC-retina_gene.tsv"),sep="\t",header=F)$V1
mat = load_mat(mat,data_path,FALSE)
obj.inbal <-CreateSeuratObject(counts = mat)
cell_types = fread(paste0(data_path,"GSE199317_ONC-retina_celltype.tsv"),sep="\t",header=TRUE)
obj.inbal$cell_type_major = cell_types$cell_type_major
saveRDS(
  object = obj.inbal,
  file = "obj.inbal.Rds",
  destdir = paste0(data_path)
)
obj.inbal <- readRDS(paste0(data_path,"obj.inbal.Rds"))
obj.inbal <- JoinLayers(obj.inbal)
obj.inbal <- subset(obj.inbal, downsample=500)
obj.inbal =calculate_cell_qc(obj.inbal,dataset="inball",figure_path,resolution=2)
obj.inbal = apply_qc_label(obj.inbal,min_genes=400,mit_cutoff=25,max_count=40000)
obj.inbal = offload_mats_to_disk(obj.inbal,data_path)
obj.inbal = RunUMAP(obj.inbal, dims = 1:30, reduction = "pca", return.model=T, reduction.name = "umap.unintegrated",seed.use = 42)
saveRDS(
  object = obj.inbal,
  file = "obj.inbal.subset.qced.Rds",
  destdir = paste0(data_path)
)
plot_cluster(obj.inbal, figure_path, "inbal_", groups=c("cell_type_major"), "umap.unintegrated")

