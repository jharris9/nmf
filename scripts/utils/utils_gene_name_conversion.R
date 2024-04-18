
########## GENE ID CONVERSION ############

#Convert ensemble gene name to all gene id 
convert_ensemble_to_gene_name = function(dt,ensb){
  dt[["ensemble"]]=toupper(rownames(dt))
  #dt = setDT(df)
  dt = merge(x=dt,y=ensb,by.x="ensemble",by.y="Gene.stable.ID")#merge
  dt = dt[,ensemble:=NULL] #drop ensemble
  dt = dt[,lapply(.SD, sum), by="Gene.name"] #sum same gene names
  dt = na.omit(dt,cols="Gene.name") #Filter na
  #dt = dt["Gene.name"!=""]
  #df = setDF(dt) 
  dt = dt %>% na.omit() #%>% filter()
  dt = dt[dt$`Gene.name`!=""]
  row.names(dt)=dt[["Gene.name"]]
  dt = dt[,"Gene.name":=NULL] #drop gene name
  return(dt)
}

convert_synonym_to_gene_name = function(dt,syn){
  dt["ensemble"]=rownames(dt)
}

#For a given gene synonym, return the gene name
get_gene_name_from_syn=function(gene_names,syn){
  DefaultAssay(obj)="RNA"
  lower_names = tolower(gene_names)
  i = match(syn,lower_names)
  if (!is.na(i)){
    return(gene_names[[i]])
  }
  if (syn %in% names(gene_names)){
    return(gene_names[[syn]])
  }else{
    return(syn)
  }
}

#Convert all the gene synonyms in a datatable to gene names and remove duplicates
convert_gene_syn_to_gene_name = function(dt,syn_path){
  syn=fread(syn_path,sep="\t",header=TRUE)
  syn = syn%>%filter(`Gene Synonym`!='')
  syn <- syn[!duplicated(syn$`Gene Synonym`), ]      
  syn_vector = syn[["Gene name"]]
  names(syn_vector) = tolower(syn[["Gene Synonym"]])
  gene_names = unlist(lapply(rownames(obj.mac), function(x) get_gene_name_from_syn(syn_vector,x)))
  dt = dt[,"name":=unlist(gene_names)]
  dt = dt[,lapply(.SD, sum), by=name] 
  row.names(dt)=dt[["name"]]
  dt = dt[,"name":=NULL]
  return(dt)
}


convert_human_to_mouse <- function(x) {
  x=tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x=sub("^Mt-","mt-",x)
  return(x)
}

convertHumanGeneList <- function(dt){
  
  genes.human <- data.frame(HGNC.symbol = rownames(dt),
                            stringsAsFactors = FALSE)
  
  require("biomaRt")
  human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://may2021.archive.ensembl.org")
  mouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "https://may2021.archive.ensembl.org")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genes.human$HGNC.symbol , mart = human, attributesL = c("mgi_symbol"), martL = mouse)
  genes.human=merge(genes.human,genesV2,all.x=T,sort=F)
  genes.human$dup=FALSE
  genes.human$dup[duplicated(genes.human$HGNC.symbol) | duplicated(genes.human$HGNC.symbol,fromLast = T)]=TRUE
  genes.human$MGI.symbol[duplicated(genes.human$HGNC.symbol)]="DUP"
  genes.human = genes.human %>% filter(MGI.symbol!="DUP" | is.na(MGI.symbol))
  genes.human$id=1:nrow(genes.human)
  genes.human = genes.human[order(genes.human$id),]
  rows_to_remove=genes.human$id[is.na(genes.human$MGI.symbol)|genes.human$dup]
  genes.huamn = genes.human %>% filter(!(id %in% rows_to_remove))
  dt$HGNC.symbol=rownames(dt)
  dt=merge(dt,genes.human)
  dt = dt%>%filter(dup==F,!is.na(MGI.symbol))
  dt$dup=NULL
  dt$HGNC.symbol=NULL
  dt = dt[,lapply(.SD, sum), by="MGI.symbol"]
  rownames(dt)=dt$MGI.symbol
  dt$MGI.symbol=NULL
  dt$id=NULL
  humanx <- genesV2[, 2]
  
  no_mouse_genes <- length(genes.human$HGNC.symbol)
  no_human_genes <- length(humanx)
  
  if(no_human_genes != no_mouse_genes){
    print("Some genes could not be translated!")
    genes_not_trans <- setdiff(genes.human$HGNC.symbol,genesV2$HGNC.symbol)
    print("These genes could not be translated:")
    print(genes_not_trans)
    print(paste("A total number of ",length(genes_not_trans),"genes could not be translated!"),sep=" ")
  }else{
    print("All genes were translated successfully!")
  }
  
  # Print all gene names that could not be translated and the number of genes that were not translated
  
  return(dt)
}