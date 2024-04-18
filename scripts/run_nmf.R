print("start")
source("scripts/utils.R")
randomly_sample_cells=function(obj,sample_perc){
  set.seed(111)
  sampled.cells <- sample(colnames(obj), size =length(colnames(obj))*sample_perc, replace=F)
  sampled.cells2 <- setdiff(colnames(obj),sampled.cells)
  obj1 <- obj[, sampled.cells]
  obj1 <- trim_misc_normalizations(obj1)
  obj2 <- obj[, sampled.cells2]
  obj2 <- trim_misc_normalizations(obj2)
  return(c(obj1,obj2))
}

plot_cor=function(obj1,obj2,filepath){
  gene_overlap=intersect(rownames(obj1@misc[["zscore"]]),rownames(obj2@misc[["zscore"]]))
  gene_av1=rowMeans(obj1@misc[["zscore"]][gene_overlap,])
  gene_av2=rowMeans(obj2@misc[["zscore"]][gene_overlap,])
  df=data.frame(g1=gene_av1,g2=gene_av2)
  p1=ggplot(df,aes(x=g1,y=g2))+geom_point()
  #p1=plot(gene_av1,gene_av2)
  gene_cor1=cor(t(as.matrix(obj1@misc[["zscore"]][1:1000,1:1000])))
  gene_cor2=cor(t(as.matrix(obj2@misc[["zscore"]][1:1000,1:1000])))
  df=data.frame(c1=gene_cor1,c2=gene_cor2)
  p2=ggplot(df,aes(x=c1,y=c2))+geom_point()
  #p2=plot(as.numeric(gene_cor1),as.numeric(gene_cor2))
  p=plot_grid(p1, p2, labels = c('A', 'B'))
  ggsave(file=paste0(filepath,"_correlation.jpg"),  device="jpeg",dpi="retina",
         units="in",height=8,width=11, plot=p)
}
compare_gene_programs_all_norms=function(normalizations,f1,f2,file_name){
  p.list=lapply(normalizations, function(x) compare_gene_programs(
    read_gp_csv(paste0(f1,"_",x,".csv")),
    read_gp_csv(paste0(f2,"_",x,".csv"))))
  p = plot_grid(plotlist = p.list) 
  ggsave(file=file_name,plot=p,dpi="retina",width=10,height=10)  
}
compare_random_samples=function(obj,dataset,figure_path,file_prefix,normalizations){
  data_path=sub("figures","processed_data",figure_path)
  data_path=sub("nmf/","",data_path)
  obj.list = randomly_sample_cells(obj,0.5) 
  obj[["randomized"]]=colnames(obj)%in%colnames(obj.list[[1]])
  plot_cluster(obj,figure_path,filename=file_prefix,reduction="umap.unintegrated",groups=c("randomized"))
  lapply(1:2,function(x) nmf_all_normalizations(obj.list[[x]],
                                                csv_file_prefix=paste0(data_path,file_prefix,"_random_sample",x,"_"),
                                                dataset=dataset,
                                                normalizations))
  file_list=paste0(data_path,file_prefix,c("_random_sample1","_random_sample2"))
  compare_gene_programs_all_norms(normalizations,file_list[[1]],file_list[[2]],file_name=paste0(figure_path,file_prefix,"_random_subset.jpg"))
}

compare_nmf=function(obj.overall,dataset,ct,normalization_method,normalizations){
  data_path=paste0("output/processed_data/",dataset,"/")
  figure_path=paste0("output/figures/",dataset,"/nmf/")
  if(normalization_method == "condition"){
    file_prefix= paste0(dataset,"_norm_condition_",ct)
    obj.cell_type = subset(obj.overall,subset=cell_type==ct)
    obj.cell_type.list = SplitObject(obj.cell_type, split.by = "condition")
    obj.cell_type.ctrl = get_condition_data(obj.cell_type.list,"Control") %>%
      clear_normalizations() %>%
      calculate_normalizations(path=paste0(data_path,ct,"/ctrl"))
    obj.cell_type.ad = get_condition_data(obj.cell_type.list,"AD") %>%
      clear_normalizations() %>%
      calculate_normalizations(path=paste0(data_path,ct,"/ad"))
  }
  if(normalization_method=="cell_type"){
    file_prefix= paste0(dataset,"_norm_cell_type_",ct) 
    obj.cell_type = subset(obj.overall,subset=cell_type==ct)
    obj.cell_type = obj.cell_type %>%
      clear_normalizations() %>%
      calculate_normalizations(path=paste0(data_path,ct))
    obj.cell_type.list = SplitObject(obj.cell_type, split.by = "condition")
    obj.cell_type.list = lapply(obj.cell_type.list, function(x) trim_misc_normalizations(x))
    obj.cell_type.ctrl = get_condition_data(obj.cell_type.list,"Control")
    obj.cell_type.ad = get_condition_data(obj.cell_type.list,"AD")
  }
  if(normalization_method=="global"){
    print("global")
    file_prefix= paste0(dataset,"_norm_global_",ct) 
    obj.cell_type = subset(obj.overall,subset=cell_type==ct)
    obj.cell_type.list = SplitObject(obj.cell_type, split.by = "condition")
    obj.cell_type.list = lapply(obj.cell_type.list, function(x) trim_misc_normalizations(x))
    obj.cell_type.ctrl = get_condition_data(obj.cell_type.list,"Control")
    obj.cell_type.ad = get_condition_data(obj.cell_type.list,"AD")
  }
  #return(obj.cell_type.ad)
  nmfs_ad=nmf_all_normalizations(obj.cell_type.ad,
                                 csv_file_prefix=paste0(data_path,file_prefix,"_ad_"),
                                 dataset=dataset,
                                 normalizations) 
  saveRDS(nmfs_ad, file = paste0(data_path,file_prefix,"_ad_nmf.Rds"))
  nmfs_ctrl=nmf_all_normalizations(obj.cell_type.ctrl,
                                   csv_file_prefix=paste0(data_path,file_prefix,"_ctrl_"),
                                   dataset=dataset,
                                   normalizations) 
  saveRDS(nmfs_ctrl, file = paste0(data_path,file_prefix,"_ctrl_nmf.Rds"))
  compare_gene_programs_all_norms(normalizations, 
                                  paste0(data_path,file_prefix,"_ctrl"),
                                  paste0(data_path,file_prefix,"_ad"),
                                  file_name=paste0(figure_path,file_prefix,"_ctrl_v_ad.jpg"))
  #compare_random_samples(obj.cell_type.ctrl,dataset,figure_path=figure_path,paste0(file_prefix,"_ctrl_"),normalizations)
  #plot_cor(obj.cell_type.ctrl,obj.cell_type.ad)
}
#Calculate all NMFs
#Compare ctrl and ad for all samples
datasets = c("cain","lee_2021","habib_2020")
#datasets=c("habib_2020")

datasets=c("cain")
#datasets = c("lee_2021","habib_2020")
cell_types = c("Astrocyte","Immune")
norms = c("cell_type")
cond=c("ctrl","ad")
#d="cain" 
#ct="Immune"
#n="cell_type"
#test=zscore(test,"~/hi")
#mat = as.matrix(test@misc[["rna_norm"]])
#mat_zscored=t(scale(t(mat), center = TRUE, scale = TRUE))
normalizations = c("rna_norm","rna_scale","rna_zscored","rna_zscored_nonneg","sct_scale","sct_scale_nonneg","softmax")
run_nmf=function(){
  for (d in datasets){
    data_path=paste0("output/processed_data/",d,"/")
    obj.overall=readRDS(file = paste0(data_path,"obj.",d,".Rds"))
    for (ct in cell_types){
      for (n in norms){
        compare_nmf(obj.overall,d,ct,n,normalizations)
        #Compare ctrl vs ad
        #b1=read_gp_csv(paste0(data_path,d,"_norm_",n,"_",ct,"_ctrl_",normalizations,".csv"))
        #b2=read_gp_csv(paste0(data_path,d,"_norm_",n,"_",ct,"_ad_",nomalizations,".csv"))
        #compare_gene_programs(b1,b2,
                              #filename=paste0("output/figures/",d,"/nmf/",d,"_norm_",n,"_",ct,"_ctrl_vs_ad.jpg"),
                              #sort=TRUE)
      }
    }
  }
}
run_nmf()
next_steps=function(){
  #Compare human vs mouse 
  human="cain"
  mouse =c("habib_2020","lee_2021")
  for (d in mouse){
    mouse_data_path=paste0("output/processed_data/",d,"/")
    human_data_path=paste0("output/processed_data/",human,"/")
    for (ct in cell_types){
      for (n in norms){
        for (c in cond){
          for (no in normalizations){
          b1=read_gp_csv(paste0(mouse_data_path,d,"_norm_",n,"_",ct,"_",c,"_",no,".csv"))
          b2=read_gp_csv(paste0(human_data_path,human,"_norm_",n,"_",ct,"_",c,"_",no,".csv"))
          compare_gene_programs(b1,b2,
                                filename=paste0("output/figures/human_v_mouse/",d,"_vs_",human,"_norm_",n,"_cond_",c,".jpg"),
                                sort=TRUE)
          }
        }
      }
    }
  }
}
#next_steps()
