source("scripts/utils.R")
source("scripts/utils/utils_data_loading.R")
dataset = "habib_2020"
figure_path = paste0("output/processed_data/figures/",dataset)
get_condition = function(l){
  l=gsub(".*AD.*","AD",l,perl=T)
  l=gsub(".*Untreated.*","AD",l,perl=T)
  l=gsub(".*_.*","Control",l,perl=T)
  return(l)
} 
load_habib = function(){
  habib_path = "/raw_data/habib_2020/GSE143758_Admouse_Hippocampus_7m_AllNuclei_UMIcounts.txt"
  habib_dt = load_dt(habib_path)
  habib_sample=unique(substr(colnames(habib_dt),1,nchar(colnames(habib_dt))-17))
  habib_sample = habib_sample[-1]
  habib_list <- c()
  for (i in 1:length(habib_sample)) {
    df_sub = habib_dt %>% dplyr::select(starts_with(habib_sample[i]))
    mat=load_mat(df_sub, paste0("output/processed_data/",dataset,"/",habib_sample[i],"_BP"),FALSE)
    habib_list[[i]]<-mat
  }
  names(habib_list) <- habib_sample
  obj <- CreateSeuratObject(counts = habib_list)
  obj$orig.ident<- habib_sample
  obj$Study = dataset
  obj[["condition"]] = get_condition(obj@meta.data[["orig.ident"]])
  return(obj)
}
obj.habib = load_habib()
qc_and_normalize(obj.habib,dataset,"brain")

obj.habib = readRDS(file = "output/processed_data/habib_2020/obj.habib_2020.Rds")

#habib_test=sample(colnames(obj.habib),size=200, replace=F)
#obj=obj.habib[, habib_test]
#obj.linnarsson <- readRDS(file = "output/processed_data/linnarsson/obj.linnarsson.subset.Rds")
#qc_and_normalize(obj,dataset,"brain")


#obj.habib <- readRDS(file = "output/processed_data/habib_2020/obj.habib_2020.Rds")
#obj.habib = calculate_cell_qc(obj.habib,dataset=dataset,resolution=2)
#obj.habib = apply_qc_label(obj.habib,min_genes=400,mit_cutoff=25,max_count=40000)
#saveRDS(
#  object = obj.habib,
#  file =  paste0("output/processed_data/",dataset,"/obj.habib.Rds")
#)
#obj.habib <- qc_subset(obj.habib,figure_path)
#obj.habib <- FindVariableFeatures(obj.habib)
#obj.habib <- identify_cell_types(obj.habib,figure_path,"habib_2020","brain")
#obj.habib[["condition"]] = get_condition(obj.habib@meta.data[["orig.ident"]])
#refined_clusters = refine_clusters(obj.habib)[-18,]
#obj.habib=label_manual_clusters(obj.habib,refined_clusters$predicted)
#plot_cluster(obj.habib,figure_path,filename="habib_et_al_clusters_",reduction="umap.unintegrated",groups=c("cell_type"))
#obj.habib = clear_normalizations(obj.habib) %>%
#  calculate_normalizations(path=paste0(data_path,"global"))
#saveRDS(
#  object = obj.habib,
#  file =paste0(data_path,"obj.habib.qced.Rds")
#)

#saveRDS(
#  object = obj.habib,
#  file =paste0(data_path,"obj.habib.Rds")
#)
#merg_habib =calculate_cell_qc(merg_habib,dataset="habib_et_al_2020",figure_path,resolution=2)
#merg_habib = apply_qc_label(merg_habib,min_genes=400,mit_cutoff=25,max_count=40000)
