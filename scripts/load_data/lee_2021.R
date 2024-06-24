source("scripts/utils.R")
source("scripts/utils/utils_data_loading.R")
dataset = "lee_et_al"
raw_data = paste0("raw_data/", dataset, "/")
figure_path = paste0("output/processed_data/figures/", dataset)

get_condition = function(condition) {
  condition = sub("^NonTg.*", "Control", condition)
  condition = sub("^PS2.*", "AD", condition)
  condition = sub("Tau.*", "AD", condition)
  return(condition)
}

load_lee = function() {
  path.raw_data <- raw_data
  path.output <- "output/processed_data/lee_et_al/"
  files.lee <- c(
    "GSM5511340_SAM24401880",
    "GSM5511341_SAM24401881",
    "GSM5511342_SAM24401883",
    "GSM5511343_SAM24401884",
    "GSM5511344_SAM24401885",
    "GSM5511345_SAM24401886",
    "GSM5511346_SAM24401887",
    "GSM5511347_SAM24401888",
    "GSM5511348_SAM24401889",
    "GSM5511349_SAM24401890",
    "GSM5511350_SAM24401891"
  )

    
  names.lee <-
    c(
      "NonTg1",
      "PS2APP1",
      "NonTg2",
      "PS2APP2",
      "TauPS2APP1",
      "NonTg3",
      "PS2APP3",
      "TauPS2APP2",
      "NonTg4",
      "PS2APP4",
      "TauPS2APP3"
    )
  obj.lee = load_data(path.raw_data, path.output, files.lee, names.lee, convert =
                        "ensemble")
  obj.lee$Study = dataset
  obj.lee[["condition"]] = get_condition(obj.lee$orig.ident)
  
  #Get cell type labels from meta data
  meta_data <-fread(file=paste0(raw_data,"GSE181786_TsneCoordinatesAndClusterAssignments.txt.gz"))

  #Ensure the meta data order matches the order in the Seurat object
  cells=sub("\\..*", "-1", Cells(obj.lee))
  ordered_index <- match(meta_data$cellID, cells)
  meta_data <- meta_data[order(ordered_index),]
  
  obj.lee <- AddMetaData(obj.lee, meta_data)
  
  
  return(obj.lee)
}

obj.lee = load_lee()
obj.lee = JoinLayers(obj.lee)
obj = qc_and_normalize(obj.lee, dataset, "brain")

obj <- FindNeighbors(obj, assay="SCT", dims = 2:25)
obj <- FindClusters(obj, assay="SCT", resolution = 1, cluster.name = "unintegrated_clusters", random.seed=42)
obj <- RunUMAP(obj, dims = 2:25, reduction = "pca", return.model=T, reduction.name = "umap.unintegrated",seed.use = 42)
DimPlot(obj,reduction="pca",group.by = "interpretation")
#obj.lee = readRDS(file = "output/processed_data/lee_2021/obj.lee_2021.Rds")


