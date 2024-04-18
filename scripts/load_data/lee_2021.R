source("scripts/utils.R")
dataset = "lee_2021"
figure_path = paste0("output/processed_data/figures/",dataset)
get_condition = function(condition){
  condition = sub("^NonTg.*","Control",condition)
  condition = sub("^PS2.*","AD",condition)
  condition = sub("Tau.*","AD",condition)
  return(condition)
}
load_lee = function(){
  path.lee <- paste0("raw_data/",dataset,"/")
  files.lee <- c("GSM5511347_SAM24401888",
                 "GSM5511350_SAM24401891",
                 "GSM5511344_SAM24401885",
                 "GSM5511340_SAM24401880",
                 "GSM5511342_SAM24401883",
                 "GSM5511345_SAM24401886",
                 "GSM5511348_SAM24401889",
                 "GSM5511341_SAM24401881",
                 "GSM5511343_SAM24401884",
                 "GSM5511346_SAM24401887",
                 "GSM5511349_SAM24401890")
  names.lee <- c("TauPS2APP1","TauPS2APP2","TauPS2APP3","NonTg1","NonTg2","NonTg3","NonTg4","PS2APP1","PS2APP2","PS2APP3", "PS2APP4")
  obj.lee = load_data(path.lee,files.lee,names.lee,convert="ensemble")
  obj.lee$Study = dataset
  obj.lee[["condition"]]=get_condition(obj.lee$orig.ident)
  return(obj.lee)
}
obj.lee=load_lee()
qc_and_normalize(obj.lee,dataset,"brain")

#obj.lee = readRDS(file = "output/processed_data/lee_2021/obj.lee_2021.Rds")
