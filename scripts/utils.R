#Server config
#dyn.load("/root/tnn/hdf5-1.12.1/hdf5/lib/libhdf5_hl.so.200")
#dyn.load("/apps/released/gcc-toolchain/gcc-4.x/gcc/gcc-9.3.0/lib64/libstdc++.so")
#.libPaths(c("/PHShome/jph35/margeta/nmf/renv/library/R-4.1/x86_64-pc-linux-gnu",
#            "/PHShome/jph35/.cache/R/renv/sandbox/R-4.1/x86_64-pc-linux-gnu/06572222"))
#lib_path = "/PHShome/jph35/R/x86_64-pc-linux-gnu-library/4.2"
#lib_path = "/PHShome/jph35/margeta/nmf/renv/library/R-4.1/x86_64-pc-linux-gnu/"
#options(bitmapType='cairo')

# suppressPackageStartupMessages({
#   library(Matrix, lib.loc = lib_path)
#   library(dplyr, lib.loc = lib_path)
#   library(sp, lib.loc = lib_path )
#   library(future, lib.loc = lib_path)
#   library(future.apply, lib.loc = lib_path)
#   library(furrr, lib.loc = lib_path)
#   library(BPCells, lib.loc = lib_path)
#   library(ggplot2, lib.loc = lib_path)
#   library(SeuratObject, lib.loc = lib_path)
#   library(Seurat,lib.loc=lib_path)
#   library(SeuratDisk, lib.loc =lib_path)
#   #library(Azimuth)
#   library(dtplyr, lib.loc = lib_path)
#   library(dittoSeq, lib.loc = lib_path)
#   library(caret, lib.loc = lib_path)
#   library(data.table, lib.loc = lib_path)
#   library(DESeq2, lib.loc=lib_path)
#   library(RcppML, lib.loc=lib_path)
#   library(tidyverse, lib.loc = lib_path)
#   library(cowplot, lib.loc = lib_path) 
#   library(png, lib.loc = lib_path)
#   library(Rfast, lib.loc = lib_path)
#   library(glmGamPoi, lib.loc = lib_path)
# })


suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(sp)
  library(future)
  library(future.apply)
  library(furrr)
  library(BPCells)
  library(ggplot2)
  library(SeuratObject)
  library(Seurat)
  #library(SeuratDisk)
  #library(Azimuth)
  library(dtplyr)
  library(dittoSeq)
  library(caret)
  library(data.table)
  library(DESeq2)
  library(RcppML)
  library(tidyverse)
  library(cowplot) 
  library(png)
  library(Rfast)
  library(glmGamPoi)
})

options(Seurat.object.assay.version = "v5")
# needs to be set for large dataset analysis
options(future.globals.maxSize = 3e10)
set.seed(42)

#Manually currated cell type markers 
markers.rod=c('Rho','Pdc','Nrl','Sag','Gnat1','Gngt1','Nr2e3','Gnb1','Pde6a','Ppef2')
markers.cone=c('Opn1mw','Arr3','Rcvrn','Gnat2','Gngt2','Opn1sw','Opn1lw','Pde6h','Guca1a')
markers.horizontal=c('Lhx1','Onecut1','Onecut2','Calb1','Adgrb1')
markers.bipolar=c('Vsx2','Otx2',"Scgn","Isl1","Gng13",'Grm6',"Cdh8","Cdh9","Tacr3","Vsx1",'Prkca','Trpm1','Grik1','Cabp5',"Tmem215","Camk2b","Pcp2")
markers.amacrine=c('Slc6a9','Tfap2a','Tfap2b','Tfap2c','Gad1','Gad2',"C1ql2")
markers.rgc=c('Slc17a6','Thy1','Nefl','Pou4f1','Pou4f2','Pou4f3','Rbpms','Nefm','Sncg','Rbpms2')
markers.muller_glia=c('Rlbp1','Glul','Apoe','Crabp1','Clu','Dkk3','Crym','Aqp4')
markers.pericyte=c('Kcnj8','Mgp','Myl9','Col4a1','Pdgfrb','Acta2')
markers.microglia=c('Cx3cr1','C1qa','C1qb','C1qc','Hexb','Ctss','P2ry12','Tmem119','B2m',"Aif1","Cd163")
markers.endothelial=c('Pecam1','Flt1','Cldn5','Igfbp7',"Adamts9","Cd34","Cdh5","Rgs5","Ly6c1")
markers.oligodendrocyte=c("Mog","Mag",'Gjc3',"Plp1",'Mbp','Mobp')
markers.astrocytes=c("Gfap",'Rbp1','Slc1a3',"Acsbg1","Sox9","Pdgfra")
markers.fibroblasts=c("Pax6","Fbn1")
markers.retina = list("Rod"=markers.rod,"Cone"=markers.cone,"Horizontal"=markers.horizontal,"Bipolar"=markers.bipolar,"Amacrine"=markers.amacrine,"RGC"=markers.rgc,"Muller Glia"=markers.muller_glia,"Pericyte"=markers.pericyte,"Microglia"=markers.microglia,"Endothelial"=markers.endothelial,"Astrocyte"=markers.astrocytes,"Fibroblast"=markers.fibroblasts)
markers.retina.core =c("Rho","Nr2e3","Opn1mw","Lhx1","Vsx2","Slc6a9","Slc17a6","Apoe","Myl9","Pecam1")

#Brain Cell Type Specific Gene Expression and Co-expression Network Architecture
markers.brain.astrocyte=c('Lrig1','Fgfr3','Sox9','Gja1','Gli3','Atp1b2','Slc1a3','Gfap')
markers.brain.ependymal=c("Foxj1","Tuba1a","Pifo","Dynlrb2")
markers.brain.microglia=c('Ccl3l3','Ccl3','Clec7a','Slc11a1','C3','Cxcr3','Ccl4','Tnf','Tmem119','Casp8','P2ry12','Ccr6')
markers.brain.oligo=c('Mbp','Plp1','Cnpase',"Mog","Mag",'Gjc3','Mobp')
markers.brain.opc=c('Pdgfra','Cspg4','Cnp','A2b5','C1ql1')
markers.brain.neuron=c('Fox3','Tubb3','Map2','Vstm2a',"Eno2",'Nefl',"Nefm",'Nefh')
markers.brain.choroid=c('Ttr','Folr1','Prlr')
markers.brain.endothelial=c('Sele','Epas1','Notch3', 'Adgrf5','Kcnj8','Anpep')
markers.brain=list("astrocyte"=markers.brain.astrocyte,"oligos"=markers.brain.oligo,'opc'=markers.brain.opc,
                   'ependymal'=markers.brain.ependymal,'neurons'=markers.brain.neuron,
                   'endothelial'=markers.brain.endothelial,'microglia'=markers.brain.microglia,
                   'choroid'=markers.brain.choroid)
markers.brain.core=c('Lrig1',"Mbp","Fox3","Ccl3l3",'Sele')

#base_path = "/Users/james/Nextcloud/Lab/Margeta_lab/analyses/"
#source("scripts/utils/utils_data_loading.R")
#source("scripts/utils/utils_gene_name_conversion.R")
#source("scripts/utils/utils_qc.R")

######### SEURAT PIPELINE ##########

cluster_data = function(obj.cluster, resolution,return.model=FALSE,run_umap=FALSE){
  if(!("data" %in% names(obj.cluster@assays[["RNA"]]@layers))){
    obj.cluster <- NormalizeData(obj.cluster)
  }
  for (l in Layers(obj.cluster, search='data')){
    obj.cluster[["RNA"]][l] <- as(obj.cluster[["RNA"]][l], Class = "dgCMatrix")
  }
  obj.cluster <- FindVariableFeatures(obj.cluster)
  obj.cluster <- ScaleData(obj.cluster,features=Features(obj.cluster))
  obj.cluster <- RunPCA(obj.cluster, seed.use = 42, npcs=25)
  obj.cluster <- FindNeighbors(obj.cluster, dims = 1:25)
  obj.cluster <- FindClusters(obj.cluster, resolution = resolution, cluster.name = "unintegrated_clusters", random.seed=42)
  if (run_umap){
    obj.cluster <- RunUMAP(obj.cluster, dims = 1:25, reduction = "pca", return.model=return.model, reduction.name = "umap.unintegrated",seed.use = 42)
  }
  return(obj.cluster)
}

project_into_reference=function(ref,query,file_prefix,figure_path,ref_data_column){
  #MERGE LAYERS
  query = JoinLayers(query)
  ref = JoinLayers(ref)
  common_features <- intersect(Features(ref), Features(query))
  ref = subset(x=ref, features=common_features)
  query = subset(x=query, features=common_features)
  ref = clear_normalizations(ref)
  query = clear_normalizations(query)
  ref = cluster_data(ref,2,return.model=TRUE,run_umap = TRUE)
  query = cluster_data(query,2,return.model=TRUE,run_umap = TRUE)
  
  anchor <- FindTransferAnchors(
    reference = ref,
    query = query,
    features=common_features,
    #reference.neighbors = NULL,
    reference.reduction = "pca",
    #normalization.method = "LogNormalize",
    dims = 1:10
  )
   
  query <- MapQuery(
    anchorset = anchor,
    reference = ref,
    query = query,
    refdata = list(celltype.l1 = ref_data_column),
    reference.reduction = "pca",
    reduction.model = "umap.unintegrated"
  )
  
  plot_cluster(obj.plot=query, figure_path, filename=file_prefix, groups=c("predicted.celltype.l1","seurat_clusters","Study"), reduction="umap.unintegrated")
  return(query)
}

load_brain_reference=function(regenerate_subset=FALSE){
  if(regenerate_subset){
    obj.linnarsson <- readRDS(file = "output/processed_data/linnarsson/obj.linnarsson.Rds")
    obj.linnarsson@assays[["RNA"]]@layers[["counts"]] = open_matrix_dir(dir="output/processed_data/linnarsson/linnarsson_bp_matrix_BP") 
    cns_subtypes=c("Neurons","Vascular","Oligos","Astrocyte","Oligos,Cycling","Ependymal","Immune","Ttr")
    obj.linnarsson <- subset(obj.linnarsson, subset=Subclass %in% cns_subtypes)
    obj.linnarsson <- cluster_data(obj.linnarsson,2,return.model=TRUE,run_umap = TRUE)
    saveRDS(
      object = obj.linnarsson,
      file = "output/processed_data/linnarsson/obj.linnarsson.subset.Rds"
    )
  }else{
    obj.linnarsson <- readRDS(file = "output/processed_data/linnarsson/obj.linnarsson.subset.Rds")
  }
  return(obj.linnarsson)
}

identify_cell_types = function(obj,figure_path,project,markers){
  obj <- cluster_data(obj,resolution=1,run_umap=TRUE,return.model=TRUE)
  if(markers=="brain"){
    markers.core = markers.brain.core
    markers.full = markers.brain
    if(!exists("obj.linnarsson")){
        obj.linnarsson=load_brain_reference()
    }
    reference = obj.linnarsson
    ref_data_column = "Subclass"
  }
  if(markers=="retina"){
    markers.core = markers.retina.core
    markers.full = markers.retina
    if(!exists("obj.mac")){
      obj.mac <- readRDS(file = paste0(base_path,"macosko_et_al/obj.mac_qced.Rds"))
    }
    reference = obj.mac
    ref_data_column = "cell_type"
  }
  if(markers=="inbal"){
    markers.core = markers.retina.core
    markers.full = markers.retina
    if(!exists("obj.inbal")){
      obj.inbal <- readRDS(file = paste0(base_path,"inbal/data/obj.inbal.subset.qced.Rds"))
    }
    reference = obj.inbal
    ref_data_column = "cell_type_major"
  }
  plot_cluster(obj,
               figure_path,
               paste0(project,"_clusters_sub_"),
               c("unintegrated_clusters"),
               reduction="umap.unintegrated")
  plot_marker_feature(obj,markers.core,
                      reduction="umap.unintegrated",
                      path=figure_path,
                      filename=paste0(project,"_gene_maps_sub.jpg"))
  plot_marker_dotplot(obj,markers.full,group.by="seurat_clusters",figure_path,paste0(project,"_marker_genes_sub.jpg"))  
  file_name = paste0(project,"_projected_into_",markers,"_")
  obj = project_into_reference(reference,
                               obj,
                               file_name,
                               figure_path,
                               ref_data_column) # This needs to be cell type specific 
  return(obj)
}

label_manual_clusters=function(obj,ids){
  obj <-SetIdent(obj,value=obj$seurat_clusters)
  names(ids) = levels(obj)
  obj <- RenameIdents(obj,ids)
  obj$cell_type=obj@active.ident
  return(obj)
}

refine_clusters=function(obj){
  dt=data.table(predicted=obj$predicted.celltype.l1,seurat=obj$seurat_clusters)
  dt.count = dt %>% 
    group_by(seurat,predicted)%>%
    summarise(N = n())%>% 
    group_by(seurat) %>%
    filter(N == max(N)) %>%
    select(seurat, predicted, max_N = N)%>%
    arrange(seurat)
  return(dt.count)
}

qc_and_normalize=function(obj,dataset,tissue_type){
  obj = calculate_cell_qc(obj,dataset=dataset,resolution=2)
  obj = apply_qc_label(obj,min_genes=400,mit_cutoff=25,max_count=40000)
  #saveRDS(
  #  object = obj,
  #  file =  paste0("output/processed_data/",dataset,"/obj.",dataset,".Rds")
  #)
  obj <- qc_subset(obj,paste0("output/figures/",dataset,"/"))
  #return(obj)
  #obj <- FindVariableFeatures(obj)
  obj <- identify_cell_types(obj,paste0("output/figures/",dataset,"/"),dataset,tissue_type)
  #obj[["condition"]]=get_condition(obj$orig.ident)
  refined_clusters = refine_clusters(obj)
  obj=label_manual_clusters(obj,refined_clusters$predicted)
  plot_cluster(obj,paste0("output/figures/",dataset,"/"),filename=paste0(dataset,"_clusters_"),reduction="umap.unintegrated",groups=c("cell_type","predicted.celltype.l1"))
  #obj = obj %>%
  #  clear_normalizations()%>%
  #  calculate_normalizations(path=paste0("output/processed_data/",dataset,"/"))
  saveRDS(
    object = obj,
    file =  paste0("output/processed_data/",dataset,"/obj.",dataset,".Rds")
  )
  return(obj)
}

######### PLOTTING FUNCTIONS ##########
plot_marker_feature = function(obj,markers,reduction,path,filename){
  plot = FeaturePlot(obj, reduction=reduction, features = markers)
  ggsave(file=filename, path=path, device="jpeg",dpi="retina",
         units="in",height=8,width=11, plot=plot)
}

plot_marker_dotplot = function(obj, markers,group.by, path, filename){
  obj = JoinLayers(obj)
  plot=DotPlot(obj,features=markers,cluster.idents = F,group.by =group.by) + RotatedAxis()
  ggsave(file=filename, path=path, device="jpeg",dpi="retina",
         units="in",height=14,width=25, plot=plot)
  #obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
}

plot_qc = function(qc.obj, path, dataset){
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
  pl.cell_qc = VlnPlot(qc.obj, group.by="orig.ident", features = features, combine=F)
  mapply(ggsave, file=paste0(dataset,"_cell_qc_", features, ".jpg"),path=path, device="jpeg",dpi="retina", plot=pl.cell_qc)
  plot <- FeatureScatter(qc.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  ggsave(file=paste0(dataset,"_cell_qc_", "nCount_vs_perc_mito", ".jpg"), path=path, device="jpeg",dpi="retina", plot=plot)
  plot <- FeatureScatter(qc.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(file=paste0(dataset,"_cell_qc_", "nCount_vs_nFeature", ".jpg"), path=path, device="jpeg",dpi="retina", plot=plot)
  plot <- ggplot(qc.obj[["percent.mt"]], aes(percent.mt)) +
    stat_ecdf(geom = "step")+
    xlab("Mitochondrial genes (%)")+
    ylab("Fraction of all cells")
  ggsave(file=paste0(dataset,"_cell_qc_", "mito_cdf", ".jpg"), path=path, device="jpeg",dpi="retina", plot=plot)
  plot <- ggplot(qc.obj[["nFeature_RNA"]], aes(nFeature_RNA)) +
    stat_ecdf(geom = "step")+
    xlab("Number of genes")+
    ylab("Fraction of all cells")
  ggsave(file=paste0(dataset,"_cell_qc_", "genes_cdf", ".jpg"), path=path, device="jpeg",dpi="retina", plot=plot)
}

plot_cluster_qc = function(obj.plot,path,dataset){
  plot_cluster(obj.plot,path,paste0(dataset,"_umap_cell_qc_"), c("orig.ident","seurat_clusters"), reduction="umap.unintegrated") 
  mit_plot=FeaturePlot(
    obj.plot,
    reduction = "umap.unintegrated",
    features = "percent.mt"
  )
  ggsave(file=paste0(dataset,"_umap_cell_qc_", "mito", ".jpg"), path=path, device="jpeg",dpi="retina", plot=mit_plot)
}

plot_cluster = function(obj.plot, path, filename, groups, reduction){
  plots=DimPlot(obj.plot, group.by = groups, reduction = reduction,label=T,
                combine = FALSE, raster=F)
  mapply(ggsave, file=paste0(filename, groups, ".jpg"),path=path,
         device="jpeg",dpi="retina", units="in",height=8,width=11, plot=plots)
}

plot_cm=function(ref, pred,file_prefix,path){
  cm <- confusionMatrix(factor(ref), factor(pred), dnn = c("Reference", "Prediction"))
  cm_mat = as.matrix(cm)
  cm_mat = melt(cm_mat %*% diag(1/colSums(cm_mat)))
  cm_df = as.data.frame(cm$table)
  cm_df$Freq=cm_mat$value
  plt <- as.data.frame(cm_df)
  cm_df$Prediction <- factor(cm_df$Prediction, levels=rev(levels(cm_df$Prediction)))
  plot = ggplot(cm_df, aes(Reference, Prediction, fill= Freq)) +
    geom_tile() + geom_text(aes(label=round(Freq,digits=2))) +
    scale_fill_gradient(low="blue", high="red") +
    labs(x = "Reference",y = "Prediction")+
    theme(axis.text.x = element_text(angle = -90, vjust = 0, hjust=0))
  ggsave(filename=paste0(file_prefix,"confusion_matrix.jpg"),path=path,plot=plot,device="jpeg",dpi="retina")
}

######## ANALYZE GRF DATA ##########
analyze_grf=function(samples, data_path, figure_path, data_set, projection){
  figure_path=paste0(figure_path,"/",data_set,"/")
  obj = load_grf_data(samples,data_path, data_set)
  obj = calculate_cell_qc(obj,dataset=data_set,resolution=2)
  #plot_cluster_qc(obj.margeta,figure_path,data_set)
  obj = apply_qc_label(obj,min_genes=400,mit_cutoff=25,max_count=40000)
  obj = JoinLayers(obj)
  obj = qc_subset(obj, figure_path)
  saveRDS(
    object = obj,
    file = paste0(data_set,".analyzed.Rds"),
    destdir = data_path)
  return(obj)
  projection_str = paste0(data_set,"_into_",projection)
  obj = identify_cell_types(obj,figure_path,projection_str,projection)
  calculate_celltype_prevalence(obj,paste0(figure_path,projection_str,".csv"))
}

calculate_celltype_prevalence=function(obj,csv_path){
  df = obj@meta.data %>% group_by(orig.ident)%>%
    count(predicted.celltype.l1)
  df = t(dcast(df,orig.ident~predicted.celltype.l1))
  write.csv(df,file=csv_path)
}

get_condition_data=function(obj.list,subgroup){
  obj <- obj.list[[subgroup]]
  DefaultAssay(obj)="RNA"
  obj <- JoinLayers(obj)
  obj[['RNA']]$counts <- as(object = obj[["RNA"]]$counts, Class = "dgCMatrix")  
  return(obj)
}
#########NMF FUNCTIONS ###########
#Some terminology
#w matrix: full matrix of weights from NMF (N genes x M factors)
#gene program: the top 50 gene names for each factor in the w matrix

#Load data
load_dataset_and_subset=function(path,subtype){
  obj = readRDS(file=path)
  obj = subset(obj,subset=cell_type==subtype)
  return(obj)
}

######### MATRIX NORMALIZATIONS ##########

#Take a seurat object and extract the SCT normalized gene expression matrix. Remove genes with no expression in any cells. Convert any negative values to 0
extract_mat=function(obj,normalization){
  DefaultAssay(obj)="RNA"
  rna_row_names=Features(obj)
  variable_features = obj@assays[["SCT"]]@var.features
  if(normalization=="sct_scale_nonneg"){
    mat = obj@assays[["SCT"]]@scale.data
    colnames(mat) = Cells(obj)
    #Row names are the original features order that are subset to just the variable features
    rownames(mat)=Features(obj)[Features(obj) %in% variable_features] 
    mat = pmax(mat,0) #replace all negative values with 0
  }else if(normalization=="sct_scale"){
    mat = obj@assays[["SCT"]]@scale.data
    colnames(mat) = Cells(obj)
    #Row names are the original features order that are subset to just the variable features
    rownames(mat)=Features(obj)[Features(obj) %in% variable_features]
  }else if(normalization=="sct_counts"){
    mat = as.matrix(obj@assays[["SCT"]]@counts)
  }else if(normalization=="sct_data"){
    mat = as.matrix(obj@assays[["SCT"]]@data)
  }else if(normalization=="rna_scale"){
    mat = as.matrix(obj@assays[["RNA"]]@layers[["scale.data"]])
    colnames(mat) = Cells(obj)
    rownames(mat) = obj@assays[["RNA"]]@meta.data[["var.features"]][!is.na(obj@assays[["RNA"]]@meta.data[["var.features"]])]
    #[["RNA"]]@meta.data[["var.features"]]
    #variable_features#rna_row_names
    #obj.lee.exp.sct_overall.immune@misc[["rna_norm"]]@dimnames[[1]]
  }else if(normalization=="rna_data"){
    mat = as.matrix(obj@assays[["RNA"]]@layers[["data"]])
    rownames(mat) = rna_row_names
    colnames(mat) = Cells(obj)
  }else if(normalization=="rna_counts"){
    mat = as.matrix(obj@assays[["RNA"]]@layers[["counts"]])
    rownames(mat) = rna_row_names
    colnames(mat) = Cells(obj)
  }else if(normalization=="rna_norm"){
    mat = as.matrix(obj@misc[["rna_norm"]])
    rownames(mat) = obj@misc[["rna_norm"]]@Dimnames[[1]]
  }else if(normalization=="rna_zscored"){
    mat = as.matrix(obj@misc[["zscore"]])
    rownames(mat) = obj@misc[["zscore"]]@Dimnames[[1]]
  }else if(normalization=="rna_zscored_nonneg"){
    mat = as.matrix(obj@misc[["zscore"]])
    rownames(mat) = obj@misc[["zscore"]]@Dimnames[[1]]
    mat = pmax(mat,0) #replace all negative values with 0
  }else if(normalization=="softmax"){
    mat = as.matrix(obj@misc[["softmax"]])
    rownames(mat) = obj@misc[["softmax"]]@Dimnames[[1]]
  }
  mat = Matrix(mat,sparse=T)
  mat = mat[rownames(mat) %in% variable_features, ] 
  #mat <- mat[rowSums(mat[, -1]) > 0, ]  #remove rows where gene expression is 0
  return(mat)
}

extract_data=function(obj,normalization){
  if(normalization=="sct_scale"){
    d=obj@assays[["SCT"]]@scale.data
    rownames(d)=Features(obj)[Features(obj) %in% obj@assays[["SCT"]]@var.features]
  }else if(normalization=="sct_counts"){
    d=obj@assays[["SCT"]]@counts
    rownames(d)=Features(obj@assays[["SCT"]])
  }else if(normalization=="sct_data"){
    d=obj@assays[["SCT"]]@data
    rownames(d)=Features(obj@assays[["SCT"]])
  }else if(normalization=="rna_data"){
    d=obj@assays[["RNA"]]@layers[["data"]]
    rownames(d)=Features(obj@assays[["RNA"]])
  }else if(normalization=="rna_counts"){
    d=obj@assays[["RNA"]]@layers[["counts"]]
    rownames(d)=Features(obj@assays[["RNA"]])
  }else if(normalization=="rna_norm"){
    d = as.matrix(obj@misc[["rna_norm"]])
    #rownames(d)=Features(obj@assays[["SCT"]])
  }else if(normalization=="rna_zscored"){
    d = as.matrix(obj@misc[["zscore"]])
    #rownames(d)=Features(obj@assays[["SCT"]])
  }else if(normalization=="softmax"){
    d = as.matrix(obj@misc[["softmax"]])
    #rownames(d)=Features(obj@assays[["SCT"]])
  }
  d=as.matrix(d)
  d=Matrix(d,sparse=T)
  return(d)
}

rna_norm=function(obj,path){
  d = as.matrix(obj@assays[["RNA"]]@layers[["counts"]])
  #d=convert_matrix_type(d,"uint32_t")
  d = Matrix(d,sparse=T)
  colnames(d) = Cells(obj)
  rownames(d) = Features(obj@assays[["RNA"]])
  d <- d[rowSums(d[, -1]) > 0, ]  #remove rows where gene exp is 0
  d = sweep(d,2,colSums(d),`/`) #normalize for sequencing depth by dividing each column by the sum of gene expression for the cell. 
  d[is.nan(d)]=0
  Misc(object = obj, slot = "rna_norm") = as(d, Class = "dgCMatrix")
  return(obj)
  write_matrix_dir(
    mat=as(d, Class = "dgCMatrix"),
    dir=path,
    overwrite=TRUE,
    compress=TRUE,
    buffer_size=8192L) 
  Misc(object = obj, slot = "rna_norm") <- open_matrix_dir(path)  
  return(obj)
}

zscore=function(obj,path){
  mat = as.matrix(obj@misc[["rna_norm"]])
  mat=t(scale(t(mat), center = TRUE, scale = TRUE))
  #mat=Rfast::transpose(Rfast::standardise(Rfast::transpose(mat), center = TRUE, scale = TRUE)) #same as above, but with Rfast
  mat[is.nan(mat)]=0
  Misc(object = obj, slot = "zscore") = as(mat, Class="dgCMatrix") 
  return(obj)
  #Dont use BP cells 
  write_matrix_dir(
    mat=as(mat, Class="dgCMatrix"),
    dir=path,
    overwrite=TRUE,
    compress=TRUE,
    buffer_size=8192L) 
  Misc(object = obj, slot = "zscore") <- open_matrix_dir(path)  
  return(obj)
}
softmax = function(par){
  n.par = length(par)
  par1 = sort(par, decreasing = TRUE)
  Lk = par1[1]
  for (k in 1:(n.par-1)) {
    Lk = max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val = exp(par - Lk)
  return(val)
}
softmax_mat = function(obj,path){
  #mat=as.matrix(obj@assays[["RNA"]]@layers[["counts"]])
  #colnames(mat) = Cells(obj)
  #rownames(mat)=Features(obj@assays[["RNA"]])
  mat = as.matrix(obj@misc[["rna_norm"]])
  #mat <- mat[rowSums(mat[, -1]) > 0, ]  #remove rows where gene exp is 0
  mat <- t(apply(mat, MARGIN = 1, FUN = softmax))
  Misc(object = obj, slot = "softmax") = as(mat, Class="dgCMatrix")
  return(obj)
  write_matrix_dir(
    mat=as(mat, Class="dgCMatrix"),
    dir=path,
    overwrite=TRUE,
    compress=TRUE,
    buffer_size=8192L) 
  Misc(object = obj, slot = "softmax") <- open_matrix_dir(path)  
}
clear_normalizations=function(obj){
  DefaultAssay(obj)="RNA"
  obj@assays[["RNA"]]@layers[["scale.data"]]=NULL
  obj@assays[["RNA"]]@meta.data[["var.features"]]=NULL
  obj@assays[["SCT"]]=NULL  
  obj@misc=list()
  return(obj)
}
calculate_normalizations=function(obj,path){
  DefaultAssay(obj)="RNA"
  obj=NormalizeData(obj)
  obj=FindVariableFeatures(obj,nfeatures = 7000)
  obj=ScaleData(obj)
  obj=rna_norm(obj,paste0(path,"rna_norm/"))
  obj=zscore(obj,paste0(path,"zscore/"))
  obj=softmax_mat(obj,paste0(path,"softmax/"))
  if (class(obj[["RNA"]]@layers[["counts"]])!="dgCMatrix"){
    obj[["RNA"]]@layers[["counts"]] <- as(obj[["RNA"]]@layers[["counts"]], Class = "dgCMatrix")
  }
  obj <- SCTransform(obj, method ="glmGamPoi",variable.features.n=7000, do.scale=T, do.center=T, vst.flavor = "v2", seed.use=42) 
  #obj=offload_mats_to_disk(obj,path)
  return(obj)
}
############# NMF and basis calculation #############

#Calculate NMF on a matrix for a given k 
calculate_nmf=function(mat,k){
  set.seed(42)
  #mat <- mat[rowSums(mat[, -1]) > 0, ]  #Remove genes with no expression in any cells
  nmf = nmf(mat,k)
  rownames(nmf[["w"]])=rownames(mat)
  colnames(nmf[["w"]])=1:25
  colnames(nmf[["h"]])=colnames(mat)
  return(nmf)
}

#Calculate the gene programs (top 50 genes, by weight in the "w" matrix) 
calculate_nmf_basis=function(nmf_list,output_path=NULL){
  basis=names(sort(nmf_list[["w"]], decreasing=T))[1:50]
  return(basis)
}

#Read in gene program CSV file
read_gp_csv=function(f){
  basis=read.csv(f)
  if("X"%in%colnames(basis)){
    basis$X=NULL
  }
  return(basis)
}

#Calculate the overlap between 2 gene programs
calculate_gene_program_overlap=function(gp1,gp2,method="jaccard"){
  if(method=="jaccard"){
    overlap = length(intersect(gp1,gp2))/(length(gp1)+length(gp2)-length(intersect(gp1,gp2)))
  }else if(method=="overlap"){
    overlap = length(intersect(gp1,gp2))/length(gp1) 
  } 
  return(overlap)
}

#Create a MxN df of all pairwise gene identity overlaps between 2 dfs that contain M and N gene factors respectively with df1 corresponding to the rows and df2 corresponding to the columns
calculate_gene_program_overlap_matrix=function(df1,df2){
  df = data.frame()
  #iterate across the columns of the first df 
  for(i in 1:dim(df1)[2]){
    row = c()
    #iterate across the columns of the second df
    for(j in 1:dim(df2)[2]){
      #calculate the overlap between these two columns. This value is a single point in the overlap matrix
      overlap=calculate_gene_program_overlap(df1[,i],df2[,j],method="overlap")
      #the first row of values will be comprised of the overlap between the first column of df1 and all the columns of df2
      row = c(row,overlap)
    }
    #repeat for all columns of df1
    df=rbind(df,row)
  }
  #columns represent the gene program from dataframe 2
  colnames(df) <- seq(1,dim(df)[1])
  #rows represent the gene program from dataframe 1
  rownames(df) <- seq(1,dim(df)[2])
  return(df)
}

#Calculate a parameter that describes the amount of gene program overlap between unreleated programs (ie stuff not on the diagonal)
quantify_program_overlap=function(df){
  m=calculate_gene_program_overlap_matrix(df,df)
  o=100*(sum(m)-sum(diag(as.matrix(m))))/((25*25-25)/2)
  return(o)
}

######## NMF PLOTS ###########

#Plot the gene program overlap matrix from "calculate_gene_program_overlap_matrix"
#Take a wide matrix (ie MxN) and convert it to 3 column df. Columns are "values": the matrix value, "row_id"and "programs" are the matrix coordinates i,j corresponding to rows and columns respectively.
plot_gene_program_overlap_matrix=function(plot_df,filename=NULL,overlap=NULL){
  plot_df=plot_df %>%
    rownames_to_column("row_id") %>% #converts the rownames to a column called row_id
    pivot_longer(-c(row_id), names_to = "programs",values_to="values")%>%
    mutate(programs = fct_relevel(programs,colnames(plot_df))) %>%
    mutate(row_id = fct_relevel(row_id,rev(rownames(plot_df)))) 
  pl=ggplot(plot_df,aes(x=programs, y=row_id, fill=values)) + 
    geom_raster() + 
    scale_fill_viridis_c(limits=range(0,1))#+
  #annotate("text",overlap,x=20, y=20,label=round(overlap,digits=2), color="white",size=16)
  if(!(is.null(filename))){
    ggsave(plot=pl,filename=filename,device="jpeg",dpi="retina")
  }
  return(pl)
}

rearrange_rows <- function(df) {
  # Get the indices of the maximum values in each column
  max_indices <- apply(df, 2, which.max)
  # Rearrange rows
  df_rearranged <- df[max_indices, , drop = FALSE]
  return(df_rearranged)
}

#Compare 2 gene programs
compare_gene_programs=function(b1,b2,filename=NULL,sort=TRUE){
  df_overlap=calculate_gene_program_overlap_matrix(b1,b2)
  overlap=NULL
  #if(b1==b2){
  #  overlap=quantify_program_overlap(b1)
  #
  if(sort){
    df_overlap=rearrange_rows(df_overlap)
  }
  if(!(is.null(filename))){
    print(filename)
    p=plot_gene_program_overlap_matrix(df_overlap,filename=filename,overlap=overlap)
  }else{
    p=plot_gene_program_overlap_matrix(df_overlap,overlap=overlap)
  }
  return(p)
}
#Take a list of gene programs and compare all of the overlaps pairwise
compare_gene_programs_pairwise=function(gp_list,filename,figure_prefix){
  pl = list()
  for(i in 1:length(gp_list)){
    for(j in 1:length(gp_list)){
      p=compare_gene_programs(gp_list[[i]],gp_list[[j]])
      #if you want to save each graph individually
      #p=compare_gene_programs(gp_list[[i]],gp_list[[j]],paste0(names(basis_list)[i],"_vs_",names(basis_list[j])))     
      pl[[length(pl)+1]]=p
    } 
  }
  plots=plot_grid(plotlist = pl) 
  ggsave(plot=plots,filename=paste0(figure_prefix,"_gene_overlap_comparision.jpg"),device="jpeg",dpi="retina",width=40,height=30)
  return(plots)
}

#Plot a heatmap of the h matrix
plot_h_mat_heatmap=function(nmf,output_prefix){
  h=nmf[["h"]]
  h=h/colSums(h)
  df=h %>%
    as.data.frame()%>%
    rownames_to_column("factors") %>% #converts the rownames to a column called factors
    pivot_longer(-c(factors), names_to = "cells",values_to="values")  
  max_cell=df %>%
    group_by(cells)%>%
    summarise(max=which.max(values))
  max_cell = max_cell%>%arrange(max)
  max_cell = max_cell$cells
  factor_ordering = as.factor(1:max(as.integer(df$factors)))
  
  p=df%>%
    ggplot(aes(x=cells, y=factors)) + 
    geom_raster(aes(fill=values)) + 
    scale_fill_viridis_c(limits=c(0,1))+
    scale_x_discrete(limits=rev(max_cell))+
    scale_y_discrete(limits=factor_ordering)+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  ggsave(filename=paste0(output_prefix,"_h_matrix_heatmap.jpg"),plot=p,device="jpeg",dpi="retina")
  return(df)
}

plot_h_mat_umap=function(obj,df,output_prefix){
  df = df %>%
    group_by(cells) %>%
    dplyr::slice(which.max(values))
  nmf_factors = df$factors
  names(nmf_factors) = df$cells
  obj$nmf_factor=factor(nmf_factors,levels=as.character(1:max(as.integer(nmf_factors))))
  p=DimPlot(obj,group.by="nmf_factor")
  #file=paste0(output_prefix, "nmf_factor_cluster.jpg")
  ggsave(filename=paste0(output_prefix,"_h_matrix_cluster.jpg"),plot=p,device="jpeg",dpi="retina")
  return(obj)
}

#Example of how to run these functions
# figure_path=paste0(base_path,"lee_et_al_2021/figures/")
# mat.sct_overall.immune.scale = extract_mat(obj.lee.exp.sct_overall.immune,"sct_scale") 
# nmf.sct_overall.immune.scale=calculate_nmf(mat.sct_overall.immune.scale,25)
# #nmf_list = calculate_nmf(mat.sct_overall.immune.scale,3)
# nmf.basis.sct_overall.immune.scale = calculate_nmf_basis(nmf.sct_overall.immune.scale,paste0(figure_path,"lee_et_al_nmf_programs_overall_sct_scale_obj.lee.ctrl.csv"))
# 
# df=plot_h_mat_heatmap(nmf.sct_overall.immune.scale[[1]],"lee_et_all_nmf_sct_overall_immune_scale",figure_path)
# obj=plot_h_mat_umap(obj.lee.exp.sct_overall.immune,df,'lee_et_all_nmf_sct_overall_immune_scale',figure_path)
# plot_h_mat_umap(obj.lee.ctrl.sct_overall.immune.reclustered,df,'lee_et_all_nmf_sct_overall_immune_scale_reclustered',figure_path)


#Sum gene expression for a given nmf factor and plot that summed expression on a umap plot. Do this for all of the different types of RNA expression normalizations. 

#Plot a feature plot (UMAP) of the summed expression of all genes in a gene list
plot_gene_exp_umap_gene_list=function(obj,gene_list,data,title=NULL){
  #gene_list = gene_list[gene_list%in%Features(obj)]
  obj@assays[["SCT"]]=NULL
  obj <- subset(obj,features=gene_list)
  data=data[rownames(data)%in%gene_list,] #Subset the data matrix by gene_list
  exp_sum=colSums(data)
  names(exp_sum)=Cells(obj)
  obj$factor_exp = exp_sum
  p=FeaturePlot(obj,features="factor_exp")
  if(!(is.null(title))){
    p=p+labs(title=title)
  }
  return(p)
}

#Plot the UMAP of summed gene expression for a list of basises for a set of data
plot_gene_exp_umap_all=function(obj,data_list, w_list, dataset){
  for(i in 1:length(data_list)){
    for(j in 1:length(w_list)){
      p = lapply(1:ncol(w_list[[j]]), function(x) plot_gene_exp_umap_gene_list(obj,w_list[[j]][,x],data_list[[i]]))
      p = plot_grid(plotlist = p) 
      ggsave(filename=paste0(dataset,"_",names(data_list[i]),names(w_list[j]),"_w.jpg"),path=paste0(figure_path,"/summed_expression_umaps/"),plot=p,device="jpeg",dpi="retina",width=40,height=40)
    }
  }
}

plot_all_expression_and_comparision_plots=function(obj,normalizations,prefix,figure_prefix,dataset){
  #data_sources=c("rna_counts","rna_data","rna_norm","rna_zscored","sct_counts","sct_data","sct_scale","softmax")  
  #data_sources=c("rna_zscored")
  gp_list=lapply(normalizations, function(x) read_gp_csv(paste0(prefix,x,".csv")))
  names(gp_list)=normalizations
  compare_gene_programs_pairwise(gp_list,prefix,figure_prefix)
  #data_list = lapply(data_sources, function(x) extract_data(obj,x))
  #names(data_list)=data_sources
  #plot_gene_exp_umap_all(obj,data_list, gp_list, dataset)
}

compare_ctrl_to_exp_for_all_noramlizations=function(dataset,cell_type,normalizations){
  p=lapply(normalizations, function(x) compare_gene_programs(
    read_gp_csv(paste0(figure_path,dataset,"_ctrl_",cell_type,"_sct_overall_",x,".csv")),
    read_gp_csv(paste0(figure_path,dataset,"_exp_",cell_type,"_sct_overall_",x,".csv"))))
  plots=plot_grid(plotlist = p, nrow=2) 
  ggsave(plot=plots,filename=paste0(dataset,"_",cell_type,"_exp_ctrl_comparision_by_normalization.jpg"),path=figure_path,dpi="retina",width=40,height=15)
}


############# NMF PIPELINE #################
#PrepSCTFindMarkers(immune.combined.sct)

trim_misc_normalizations=function(obj){
  list.misc = lapply(names(obj@misc), function(x) obj@misc[[x]][,which(colnames(obj@misc[[x]])%in%Cells(obj))])
  names(list.misc)=names(obj@misc)
  obj@misc=list.misc
  return(obj)
}

run_nmf_pipeline=function(obj.nmf,normalization,csv_file,figure_prefix,dataset,obj.recluster=NULL){
  #if(file.exists(csv_file)){
  #  return()
  #}
  mat = extract_mat(obj.nmf,normalization)
  print("Calculating NMF")
  nmf =  calculate_nmf(mat,25)
  #basis = calculate_nmf_basis(nmf)
  #sort the columns from highest to lowest and take the top 50
  basis=apply(nmf[["w"]], 2, function(y) names(sort(y, decreasing = T))[1:50])
  gp_prevalence=rep(0,25)
  names(gp_prevalence)=as.character(c(1:25))
  gp_prevalence[names(table(max.col(t(nmf[["h"]]))))]=table(max.col(t(nmf[["h"]])))
  gp_prevalence=names(sort(gp_prevalence,decreasing = T)) # compute the dominant gene program in each cell and tally the number of occurances of each gene program. sort the gene programs by prevalence as the dominant program in each cell
  write.csv(basis,file=csv_file,col.names=F) #Write the gene programs to file
  #calculate_gene_program_overlap_matrix(basis,basis) %>%
  #plot_gene_program_overlap_matrix(paste0(file_prefix,"_",normalization))
  #Analyze info from h matrix
  df.h=plot_h_mat_heatmap(nmf,figure_prefix) #Plot H heatmap
  plot_h_mat_umap(obj.nmf,df.h,figure_prefix)
  if(!(is.null(obj.recluster))){
    plot_h_mat_cluster(obj.recluster,df.h,figure_prefix)
  }
  return(nmf)
}

run_nmf_parallel <- function(normalization,obj.nmf,csv_file_prefix,figure_prefix,dataset){
  print(normalization)
  nmf = run_nmf_pipeline(
    obj.nmf,
    normalization,
    csv_file = paste0(csv_file_prefix, normalization, ".csv"),
    figure_prefix = paste0(figure_prefix,normalization),
    dataset
  )
  return(nmf)
}


#Calculate NMF for all normalization
nmf_all_normalizations=function(obj.nmf,csv_file_prefix,dataset,normalizations,obj.reclustered=NULL){
  print("Running NMF pipeline")
  figure_prefix = sub("processed_data","figures",csv_file_prefix)
  figure_prefix = sub(paste0("/",dataset,"/"),paste0("/",dataset,"/nmf/"),figure_prefix)
  nmfs = list()
  plan(sequential)#clear old workers that might still be running 
  #plan("future::multisession", workers=6)
  #packages=search()
  #nmfs = future_map(normalizations, 
  #           run_nmf_parallel,
  #           obj.nmf=obj.nmf,
  #           csv_file_prefix=csv_file_prefix,
  #           figure_prefix=figure_prefix, 
  #           dataset=dataset,
             #.options = furrr_options(packages = packages),
  #           .progress=TRUE)
  #plan(sequential)
  for (normalization in normalizations){
    print(normalization)
    if(!(is.null(obj.reclustered))){
      nmf=run_nmf_pipeline(obj.nmf,
                           normalization,
                           csv_file=paste0(csv_file_prefix,normalization,".csv"),
                           figure_prefix=figure_prefix,
                           dataset,
                           obj.reclustered)
    }else{
      nmf=run_nmf_pipeline(obj.nmf,
                           normalization,
                           csv_file=paste0(csv_file_prefix,normalization,".csv"),
                           figure_prefix=figure_prefix,
                           dataset)
    }
    nmfs[paste0(csv_file_prefix,normalization)]=nmf
  }
  plot_all_expression_and_comparision_plots(obj.nmf,normalizations,csv_file_prefix,figure_prefix,dataset)
  return(nmfs)
}

