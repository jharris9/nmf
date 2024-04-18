source("scripts/utils.R")
dataset = "cain"
figure_path = paste0("output/processed_data/figures/",dataset)

get_condition = function(l){
  l=gsub(".*AD0.*","Control",l,perl=T)
  l=gsub(".*AD1.*","AD",l,perl=T)
  l=gsub(".*Path0.*","Control",l,perl=T)
  l=gsub(".*Path1.*","AD",l,perl=T)
  return(l)
}

load_cain = function(){
  genes=read.csv("raw_data/cain/ROSMAP_Brain.snRNAseq_metadata_genes_20230420.csv")$x
  cells=read.csv("raw_data/cain/ROSMAP_Brain.snRNAseq_metadata_cells_20230420.csv",header=T)
  mm=readMM("raw_data/cain/ROSMAP_Brain.snRNAseq_counts_sparse_format_20230420.csv")
  cain_dt=setDT(as.data.frame(mm))[]
  rownames(cain_dt)=genes
  colnames(cain_dt)=cells$cell_name
  cain_dt=convertHumanGeneList(cain_dt)
  obj.cain = CreateSeuratObject(counts = cain_dt)
  obj.cain$Study = "Cain"
  obj.cain$sample = cells$specimenID
  obj.cain$study_cell_type = cells$specimenID
  obj.cain$study_subtype = cells$sub.cluster
  cell_type=cells$Cell.Type
  cell_type=gsub("microglia","Immune",cell_type,perl=T)
  cell_type=gsub("astrocytes","Astrocyte",cell_type,perl=T)
  obj.cain$cell_type=cell_type
  obj.cain[["condition"]] = get_condition(obj.cain@meta.data[["sample"]])
  return(obj.cain)
}
#obj.cain = load_cain()
#qc_and_normalize(obj.cain,dataset,"brain")

obj.cain = readRDS(file = "output/processed_data/cain/obj.cain.Rds")
#obj.cain=FindVariableFeatures(obj.cain,nfeatures = 7000)
#obj.cain=ScaleData(obj.cain)
#obj.cain = offload_mats_to_disk(obj.cain,"output/processed_data/cain/")
#Misc(object = obj.cain, slot = "rna_norm") <- open_matrix_dir("output/processed_data/cain/rna_norm/") 
obj.cain = zscore(obj.cain)

mat = as.matrix(obj.cain@misc[["rna_norm"]])
#mat=t(scale(t(mat), center = TRUE, scale = TRUE))
mat=Rfast::transpose(Rfast::standardise(Rfast::transpose(mat), center = TRUE, scale = TRUE)) #same as above, but with Rfast
mat[is.nan(mat)]=0
#obj@assays[["RNA"]]@layers[["zscore"]]=mat
mat=Matrix(mat,sparse=F)
mat=as(mat, "CsparseMatrix")
write_matrix_dir(
  mat=mat,#as(mat, Class="BPCells::Iterablematrix"),
         dir="output/processed_data/cain/zscore/",
         overwrite=TRUE,
         compress=TRUE,
         buffer_size=8192L)
  Misc(object = obj.cain, slot = "zscore") <- mat
    open_matrix_dir(path)  

obj.cain <- SCTransform(obj.cain, method ="glmGamPoi",variable.features.n=7000, do.scale=T, do.center=T, vst.flavor = "v2", seed.use=42) 
obj.cain = obj.cain %>%
  clear_normalizations()%>%
  calculate_normalizations(path=paste0("output/processed_data/",dataset,"/"))
saveRDS(
  object = obj.overall,
  file =  "output/processed_data/cain/obj.cain.Rds")
#mat=as(obj.cain[["RNA"]]["data"], Class = "dgCMatrix")
#mat=convert_matrix_type(mat,type="uint32_t")@matrix

#write_matrix_dir(
#  mat=mat,#as(obj.cain[["RNA"]]["counts"], Class = "dgCMatrix"),
#  dir="output/processed_data/cain/counts/",
#  overwrite=TRUE,
#  compress=TRUE,
#  buffer_size=8192L) 

#t=open_matrix_dir("output/processed_data/cain/test/") 
old=function(){
  obj.cain=readRDS(paste0(data_path,"obj.cain.Rds"))
  obj.cain <- FindVariableFeatures(obj.cain)
  obj.cain <- identify_cell_types(obj.cain,figure_path,project.cain,"brain")
  
  cain.cluster.ids=c("Oligo","Oligo","Astrocyte","Neurons","Neurons","Neurons","Neurons","Neurons","Astrocyte","Astrocyte","Neurons","Oligos,Cycling","Neurons","Neurons","Neurons","Neurons","Neurons","Neurons","Immune","Neurons","Vascular","Neurons","Neurons","Astrocyte","Neurons","Neurons","Neurons","Vascular","Neurons")
  obj.cain = label_manual_clusters(obj.cain,cain.cluster.ids)
  saveRDS(
    object = obj.cain,
    file = paste0(data_path,"obj.cainRds")
  )
  
  
  normalizations = c("rna_counts","rna_data","rna_norm","rna_scale","rna_zscored","rna_zscored_nonneg","sct_counts","sct_data","sct_scale","sct_scale_nonneg")
  
  plot_cluster(obj.cain,figure_path,filename="cain_et_al_clusters",reduction="umap.unintegrated",groups=c("cell_type","study_subtype"))
  
  expand_cell_type=function(obj){
    df=obj@meta.data
  }
  
  #Calculate NMF after normalization within cell type and condition
  
}