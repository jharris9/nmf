
############### DATA LOADING ###############

#Quickly load a datatable from file
load_dt  = function(path){
  dt = fread(path,blank.lines.skip = TRUE)
  #dt = dt[,lapply(.SD, sum), by="Gene name"]
  rows = dt[[1]]
  rownames(dt)=rows
  dt[,1:=NULL]
  return(dt)
}

#Takes a datatable of single cell expression data and writes it to a BPcells MatrixDir
load_mat = function(dt,path,convert){
  if (convert=="ensemble"){
    #mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "mouse")
    ensb=read.csv("raw_data/ensembl_gene_name.txt",sep="\t")
    dt = convert_ensemble_to_gene_name(dt,ensb)
  }else if(convert=="human"){
    #rownames(dt) = unlist(lapply(rownames(dt), function(x) convert_human_to_mouse(x)))
    dt = convertHumanGeneList(dt)
  }
  mat = Matrix(as.matrix(dt),sparse=TRUE)
  rows = rownames(dt)
  mat@x <- as.integer(mat@x)
  mat = convert_matrix_type(mat,type="uint32_t")
  mat = mat@matrix@mat
  rownames(mat)=rows
  write_matrix_dir(
    mat = mat,
    dir = path,
    overwrite = TRUE
  )
  mat <- open_matrix_dir(dir=path)
  return(mat)
}
load_data = function(dir,files,samples,convert){
  data.list <- c()
  metadata <- c()
  for (i in 1:length(files)) {
    path <- paste0(dir,files[i],".txt.gz")
    dt = load_dt(path)
    mat = load_mat(dt,paste0(gsub(".txt.gz", "", path), "_BP"),convert)
    #ensure no duplicate cells 
    colnames(mat)=paste0(colnames(mat),samples[[i]]) 
    data.list[[i]]<-mat
    metadata=append(metadata,rep(samples[[i]],length(colnames(mat))))
  }
  names(data.list) <- samples
  merged.object <- CreateSeuratObject(data.list)
  merged.object$orig.ident<- metadata
  return(merged.object)
}

load_grf_data = function(samples,base_path,study){
  data.list <- c()
  sample_id <- c()
  for (i in 1:length(samples)){
    sample_path <- paste0(base_path,samples[i],"/filtered_feature_bc_matrix.h5")
    print(sample_path)
    mat <- Read10X_h5(sample_path,T,T)
    colnames(mat) = paste0(samples[i],"_",colnames(mat))
    path <- paste0(gsub(".h5", "", sample_path), "_BP")
    write_matrix_dir(
      mat = mat,
      dir = path,
      overwrite = TRUE
    )
    mat <- open_matrix_dir(dir=path)
    data.list[[i]] <- mat
    sample_id = c(sample_id,rep(samples[i],length(colnames(mat))))
  }
  names(data.list) <- samples
  merged.object <- CreateSeuratObject(counts = data.list)
  merged.object$orig.ident<- sample_id
  merged.object$Study = study
  return(merged.object)
}

load_study = function(files,names,project,convert){
  path <- paste0(base_path,project,"/data/")
  obj = load_data(path,files,names,convert=convert)
  obj$Study = project
  obj =calculate_cell_qc(obj,dataset=project,paste0(base_path,project,"/figures/"),resolution=2)
  obj = apply_qc_label(obj,min_genes=400,mit_cutoff=25,max_count=40000)
  return(obj)
}

load_mats_to_mem=function(obj,path){
  for (l in Layers(obj)){
    obj[["RNA"]][l] <- as(obj[["RNA"]][l], Class = "dgCMatrix")
  }
  return(obj)
}
offload_mats_to_disk=function(obj,path){
  for (l in Layers(obj)){
    write_matrix_dir(
      mat=as(obj[["RNA"]][l], Class = "dgCMatrix"),
      dir=paste0(path,l),
      overwrite=TRUE,
      compress=TRUE,
      buffer_size=8192L) 
    obj[["RNA"]][l]=open_matrix_dir(paste0(path,l))  
  }
  return(obj)
}