
############## QC FUNCTIONS ##############

#Calculate percent mitochondrial reads, number of RNAs / cell, number of features (genes) per cell 
calculate_qc_data = function(obj.qc){
  obj.tmp = JoinLayers(obj.qc)
  obj.tmp[["percent.mt"]] <- PercentageFeatureSet(obj.tmp, pattern = "^mt-")
  obj.tmp[["nCount_RNA"]] = colSums(x = obj.tmp, slot = "counts") 
  obj.tmp[["nFeature_RNA"]] = colSums(x = GetAssayData(object = obj.tmp, slot = "counts") > 0) 
  #obj.tmp[["RNA"]] <- split(obj.tmp[["RNA"]], f = obj.tmp$orig.ident)
  return(obj.tmp)
}

#Calculate QC metrics and plot them
calculate_cell_qc = function(obj, dataset, resolution){
  path = paste0("output/figures/",dataset,"/")
  obj = calculate_qc_data(obj)
  plot_qc(obj,path,dataset)
  obj = cluster_data(obj,resolution,run_umap=TRUE,return.model=TRUE)
  plot_cluster_qc(obj,path,dataset)
  return(obj)
}

apply_qc_label=function(obj,min_genes,mit_cutoff,max_count){
  obj$qc = obj$percent.mt<mit_cutoff & obj$nCount_RNA<max_count & obj$nFeature_RNA>min_genes
  return(obj)
}

#Remove low quality clusters based on average mitochondrial reads and number of genes detected
filter_clusters = function(obj.cluster, path, filename, max_mito, min_gene,resolution){
  #obj.cluster <- cluster_data(obj.cluster,resolution=resolution)
  pl.cluster_qc = VlnPlot(obj.cluster, group.by = "seurat_clusters", features = c("nFeature_RNA", "nCount_RNA",
                                                                                  "percent.mt"), ncol = 1) 
  ggsave(filename=filename,path=path,plot=pl.cluster_qc,
         
         device="jpeg",dpi="retina")
  good_clusters=as.data.table(obj.cluster@meta.data)
  #good_clusters=good_clusters%>%
  #  group_by(seurat_clusters) %>%
  #  summarise(mean_mito=mean(percent.mt),mean_count=mean(nCount_RNA),mean_gene=mean(nFeature_RNA))%>%
  #  filter(mean_mito<max_mito,mean_gene>min_gene)
  good_clusters=good_clusters[, .(mean_mito = mean(percent.mt), mean_count = mean(nCount_RNA), mean_gene = mean(nFeature_RNA)),
                              by = seurat_clusters][mean_mito < max_mito & mean_gene > min_gene]
  good_clusters = good_clusters$seurat_clusters 
  obj.cluster$qc_cluster = obj.cluster$seurat_clusters %in% good_clusters
  return(obj.cluster)
}

qc_subset = function(obj,figure_path){
  print("subset")
  obj <- subset(obj, subset = qc)
  print("filter")
  obj<- filter_clusters(obj, 
                        figure_path,
                        "cluster_qc.jpg", 
                        max_mito=10, 
                        min_gene=400, 
                        resolution=2)
  print("subset again")
  obj <- subset(obj, subset = qc_cluster)
  return(obj)
}

