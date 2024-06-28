library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(Matrix)


#Load 
obj.lee = readRDS(file = "output/processed_data/lee_et_al/obj.lee_2021.Rds")
obj.lee.microglia <- subset(obj.lee, subset=interpretation=="Microglia")
rm(obj.lee)
obj.lee.microglia.list <- SplitObject(obj.lee.microglia,split.by="Genotype")
rm(obj.lee.microglia)
gc()


#Output the filtered matrix to disk.
convert_seurat_to_files <- function(obj, dir) {
  counts <- obj@assays$RNA$counts
  barcodes <- colnames(counts)
  gene_names <- rownames(counts)
  
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  # Output counts matrix
  writeMM(counts, paste0(dir, 'matrix.mtx'))
  
  # Output cell barcodes
  write.table(as.data.frame(barcodes), paste0(dir, 'barcodes.tsv'),
              col.names = FALSE, row.names = FALSE, sep = "\t")
  
  # Output feature names
  features <- data.frame("gene_id" = gene_names, "gene_name" = gene_names, type = "Gene Expression")
  write.table(as.data.frame(features), sep = "\t", paste0(dir, 'genes.tsv'),
              col.names = FALSE, row.names = FALSE)
}

run_cNMF_analysis <- function(dir, runname, workers, k) {
  #prepare step which normalizes the count matrix and prepares the factorization step
  k_values <- paste(5:k, collapse = " ")
  prepare_cmd <- paste("cnmf prepare --output-dir", dir,
                       "--name", runname,
                       "-c", paste0(dir, 'matrix.mtx'),
                       "--max-nmf-iter 2000",
                       "-k", k_values, "--n-iter 20", sep=" ")
  print(prepare_cmd)
  system(prepare_cmd)

  #Factorization step runs NMF --n-iter times (in this case 20) for each value of K. In this tutorial we run all these jobs sequentially on a single worker but in theory this can be distributed to multiple cores or nodes with separate commands like so:
  for (i in 0:(workers-1)) {
    factorize_cmd <- paste("cnmf factorize --output-dir", dir,
                           "--name", runname,
                           "--worker-index", i, "--total-workers", workers, sep=" ")
    print(factorize_cmd)
    system(factorize_cmd)
  }
  
  #Concatenate the results for each value of K into a single file
  combine_cmd <- paste("cnmf combine --output-dir", dir,
                       "--name", runname, sep=" ")
  print(combine_cmd)
  system(combine_cmd)
  
  # Plot for k selection
  k_selection_cmd <- paste("cnmf k_selection_plot --output-dir", dir,
                           "--name", runname, sep=" ")
  print(k_selection_cmd)
  system(k_selection_cmd)
}


run_consensus_clustering <- function(obj, dir, runname, k, local_density_threshold) {
  # Consensus clustering step
  cmd <- paste("cnmf consensus --output-dir", dir,
               "--name", runname,
               '--components', k,
               '--local-density-threshold', local_density_threshold,
               '--show-clustering', sep=" ")
  print(cmd)
  system(cmd)
  
  # Load resulting files
  usage_file <- file.path(dir, runname, paste(runname, "usages", paste0("k_",k,".dt_0_1"), 'consensus', 'txt', sep="."))
  spectra_score_file <- file.path(dir, runname, paste(runname, "gene_spectra_score", paste0("k_",k,".dt_0_1"), 'txt', sep="."))
  spectra_tpm_file <- file.path(dir, runname, paste(runname, "gene_spectra_tpm", paste0("k_",k,".dt_0_1"), 'txt', sep="."))
  
  usage <- read.table(usage_file, sep='\t', row.names=1, header=TRUE)
  spectra_score <- read.table(spectra_score_file, sep='\t', row.names=1, header=TRUE)
  spectra_tpm <- read.table(spectra_tpm_file, sep='\t', row.names=1, header=TRUE)
  print(head(usage))
  
  # Normalize usage file so that each cell sums to 
  usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))
  
  # Merge usage_norm into Seurat object metadata

  new_metadata <- merge(obj@meta.data, usage_norm, by = "row.names", all.x = TRUE)
  rownames(new_metadata) <- new_metadata$Row.names
  obj@meta.data <- new_metadata
  return(obj)
}

#Run the standard Seurat UMAP pipeline so we can plot the GEP usages over the UMAP
run_seurat <- function(obj) {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  all_genes <- rownames(obj)
  obj <- ScaleData(obj, features = all_genes)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj),npcs=25)
  obj <- FindNeighbors(obj, dims = 1:15)
  obj <- RunUMAP(obj, dims = 1:15)
  
  return(obj)
}

obj=obj.lee.microglia.list[[2]]
dir="output/processed_data/lee_et_al/"
runname=paste0(names(obj.lee.microglia.list[2]),"_cNMF")
convert_seurat_to_files(obj, dir)
run_cNMF_analysis(dir,,6,20)
run_consensus_clustering(obj,dir,runname,,0.1)
obj=run_seurat(obj)
p <- FeaturePlot(obj, features = "X1", combine=F)
p

#To help make sense of the learned GEPs, we extract the top 20 most highly weighted genes for each GEP as below. 

get_top_colnames <- function(row) {
  # Orders the values in descending order and gets the names of the top 20
  print(row[1:5])
  top_indices <- order(row, decreasing = TRUE)[1:20]
  return(colnames(spectra_score)[top_indices])
}

top_colnames <- apply(spectra_score, 1, get_top_colnames)
top_colnames <- as.data.frame(top_colnames)

top_colnames



filtered_dir = './R_Example_Data/filtered/'