suppressPackageStartupMessages({
  library(loomR)
})
source("scripts/utils.R")
#base_path <- "/Users/james/Nextcloud/Lab/Margeta_lab/analyses/"
figure_path <- "output/figures/linnarsson/"

#linnarsson_loom <- connect(filename = "raw_data/linnarsson/l5_all.agg.loom", mode = "r+", skip.validate=TRUE)
#test = Convert(source = linnarsson_loom, to = "loom", filename = #"pbmc_small.loom", 
#    display.progress = FALSE)


remove_duplicate_rows = function(df){
  dups = df[row.names(df) %in% duplicated(row.names(df)),]
  dups = rowsum(dups, as.factor(rownames(dups)))
  df = df[!duplicated(row.names(df)),]
  df = rbind(df,dups)
  return(df)
}

load_linnarsonn = function(){
  loom_linnarsson <- connect(filename = "/raw_data/linnarsson/l5_all.loom", mode="r+")
  #rm(mat_linnarsson)
  mat_linnarsson = loom_linnarsson[["matrix"]][,]
  #mat_linnarsson = convert_matrix_type(mat_linnarsson,type="uint32_t")
  #mat_linnarsson <- as.integer(mat_linnarsson)
  mat_linnarsson = Matrix(mat_linnarsson,sparse=TRUE)
  gene.names <- loom_linnarsson[["row_attrs/Gene"]][]
  cell.names = loom_linnarsson[["col_attrs/CellID"]][]
  cell.names <- paste0(cell.names,1:length(cell.names))
  colnames(mat_linnarsson) = gene.names
  rownames(mat_linnarsson) = cell.names
  mat_linnarsson=t(mat_linnarsson)
  mat_linnarsson=remove_duplicate_rows(mat_linnarsson)
  row.names(mat_linnarsson)[duplicated(row.names(mat_linnarsson))]
  write_matrix_dir(
    mat = mat_linnarsson,
    dir = "output/processed_data/linnarsson/linnarsson_bp_matrix_BP",
    overwrite = TRUE
  )
  metadata.list = names(loom_linnarsson[["col_attrs"]])
  metadata = data.frame(matrix(ncol=0,nrow=length(cell.names)))
  for (i in 1:length(metadata.list)){
    metadata[metadata.list[i]] <- loom_linnarsson[[paste0("col_attrs/",metadata.list[i])]][]
  }
  rownames(metadata) = cell.names
  obj.linnarsson = open_matrix_dir(dir="output/processed_data/linnarsson/linnarsson_bp_matrix_BP")
  obj.linnarsson = CreateSeuratObject(counts = obj.linnarsson, meta.data = metadata)
  obj.linnarsson = calculate_cell_qc(obj.linnarsson,dataset="linnarsson",resolution=2)
  obj.linnarsson = apply_qc_label(obj.linnarsson,min_genes=400,mit_cutoff=25,max_count=40000)
  obj.linnarsson = JoinLayers(obj.linnarsson)
  saveRDS(
    object = obj.linnarsson,
    file = "output/processed_data/linnarsson/obj.linnarsson.Rds",
  )
  obj.linnarsson <- readRDS(file = "output/processed_data/linnarsson/obj.linnarsson.Rds")
  obj.linnarsson@assays[["RNA"]]@layers[["counts"]] = open_matrix_dir(dir="output/processed_data/linnarsson/linnarsson_bp_matrix_BP") 
  p1 = DimPlot(obj.linnarsson, group.by = c("ClusterName","Subclass"), reduction = "umap.unintegrated", combine=T, raster=F, label=T, label.size = 8, repel=T)  
  p1=unlist(p1)
  p1 = p1+ guides(col = guide_legend(ncol=5))
  mapply(ggsave, file=paste0("linnarsson_",c("ClusterName","Subclass"),".jpg"),path=figure_path,
         device="jpeg",dpi="retina", units="in",height=20,width=30, plot=p1)
  
}

load_linnarsonn() 
