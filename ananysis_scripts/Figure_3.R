# get /dragennfs/area8/datnguyen/PRS_arrays_project/KHV_10_folds_imputed/PRS_result/*all_score /media/datn/data/1st_DC_PRS_array_project/PRS_codes/PRS_khv

setwd("/media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores")


# get /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS/PRS_WGS/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores
# get /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008/PRS_result/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores

# get /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR/PRS_WGS/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores
# get /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed/PRS_result_all_pop/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores

require(data.table)
# wgs_height =  wgs =  fread("WGS_VNP_HEIGHT.all_score")
# wgs_bmi = fread("WGS_VNP_BMI.all_score")
# arr = fread("Axiom_JAPONICA_HEIGHT.all_score")
# arr = fread("global-screening-array-v.3_HEIGHT.all_score")
# col = names(wgs)[6]
# cor(wgs[,..col], arr[,..col])

get_cor_array <- function(wgs, array, pop, trait){
  #wgs_path = paste0("WGS_", pop, "_", trait, ".all_score" )
  #wgs = fread()
  cols = names(wgs)[-c(1,2)]
  array_path = paste0(array, "_", pop, "_", trait, ".all_score")
  arr = fread(array_path)
  res = c()
  for(col in cols){
    tem = cor(wgs[,..col], arr[,..col])
    res = c(res, tem)
  }
  res2 = c(array, trait, res, mean(res), sd(res))
  names(res2) = c("array", "trait", cols, "mean", "sd")
  return(res2)
}


#get_cor_array(wgs, "Axiom_JAPONICA", "VNP", "HEIGHT")

array_list = c("Axiom_GW_ASI", "Axiom_GW_CHB", "Axiom_GW_EUR", "Axiom_GW_PanAFR", "Axiom_JAPONICA", "infinium-omnizhonghua-v1.4", "infinium-psycharray-v1.3", "japanese-screening-array-v1.0", "multi-ethnic-eur-eas-sas-v1.0", "multi-ethnic-global-v1.0", "oncoarray-500k", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB_WCSG", "chinese-genotyping-array-v1.0", "cytosnp-850k-v1.2", "global-screening-array-v.3", "human-cytosnp-12-v2.1", "infinium-core-v1.2", "infinium-global-diversity-array-v1.0", "infinium-omni2.5.v1.5", "infinium-omni5-v1.2", "GenomeWideSNP_6.0")





# "Breast_cancer"
get_cor_pop <- function(pop, array_list, trait_list = c("HEIGHT", "BMI", "Type_2_diabetes")){

  res = list()

  for( trait in trait_list){
    wgs_path = paste0("WGS_", pop, "_", trait, ".all_score" )
    wgs = fread(wgs_path)
    all_cor = list()
    for( array in array_list){
    tem = get_cor_array(wgs, array, pop, trait)
    all_cor[[array]] = tem
    }
    res[[trait]] = as.data.frame(do.call(rbind, all_cor))
  }
  df = as.data.frame(do.call(rbind, res))
  df$pop = pop
  return(df)
}





#df = get_cor_pop("VNP", array_list)
pop_list = c("EAS", "EUR", "SAS", "AMR", "AFR", "VNP")
res_all = list()

for(pop in pop_list){
  res_all[[pop]] = get_cor_pop(pop, array_list)
}

df = as.data.frame(do.call(rbind, res_all))


array_size = fread("../array_size.csv")

#
array_size$size = round(array_size$No.Assays / 1000, 0)

array_size$array = array_size$'Short name'

array_size$array_rename = paste0(array_size$array, "_", array_size$size, "k")

df$short_name = array_size$array[match(df$array, array_size$old_name)]

df$array_rename = array_size$array_rename[match(df$array, array_size$old_name)]

#df$array <- factor(df$array, levels = array_size$array)
df$array_rename <- factor(df$array_rename, levels = array_size$array_rename)


df$mean = as.numeric(df$mean)
df$sd = as.numeric(df$sd)

df2 = df
df2$trait[df2$trait == "HEIGHT"] = "Height" 
require(ggplot2)
p = ggplot(data = df2, aes(x = array_rename, y = mean, fill = trait, color = trait)) + geom_point(stat="identity", position = position_dodge(.5)) + geom_errorbar(aes(x = array_rename, ymin = mean - sd, ymax = mean + sd), width = 0, position = position_dodge(.5)) + guides(x = guide_axis(angle = 90)) + scale_y_continuous(breaks=seq(0,1,0.01), limits = c(0.8,1) ) + theme_light() + ylab("PGS correlation") + xlab("SNP array")  + facet_wrap(~ pop, nrow = 2) + theme(legend.position="bottom", axis.title.x = element_blank())

pdf(file= "../../Figure_3.pdf",  width=12, height=9)
p
dev.off()









### get table statistics

get_stat_pop_trait <- function(df, trait, pop_list = c("AFR", "AMR", "EAS", "EUR", "SAS", "VNP")){
  #pop = "EAS"
  #trait = "HEIGHT"
  res = list()
  for(pop in pop_list){
      pick = df$trait == trait & df$pop == pop
      tem = df[pick,]
      tem$mean = round(tem$mean, 4)
      tem$sd = round(tem$sd, 4)
      tem$mean_sd = paste(tem$mean, tem$sd, sep = "±")
      pick_col = c("short_name", "mean_sd")
      tem = tem[,pick_col]
      res[[pop]] = tem
  }
  res_df = as.data.frame(do.call(cbind, res))
  pick = seq(1,ncol(res_df), by = 2)[-1]
  res_df = res_df[,-pick]
  names(res_df)[1] = "short_name"
  return(res_df)
}

res_all = list()
trait_list = c("HEIGHT", "BMI", "Type_2_diabetes")
for(trait in trait_list){
  res_all[[trait]] = get_stat_pop_trait(df, trait) 
}


pop_list = c("EAS", "EUR", "SAS", "AMR", "AFR", "VNP")
for(pop in pop_list){
  res_all[[pop]] = get_stat_pop_trait(df, pop) 
}

# for sorting
require(data.table)
array_size = fread("../array_size.csv")

export_trait <- function(trait, res_all, array_size){
  array_size$size = round(array_size$No.Assays / 1000, 0)
  array_size$array = array_size$'Short name'
  res = res_all[[trait]] 
  res2 = res_all[[trait]] 
  res2$size = array_size$No.Assays[match(res2$short_name, array_size$array)]
  od = order(res2$size)
  res = res[od,]
  out_name = paste0("../../",  trait ,"_PRS_cor_12_cutoff.csv")
  fwrite(res, out_name, sep = ",")
}

for(trait in trait_list){
  export_trait(trait, res_all, array_size)
}
































# # pdf(file= "PRS_acc_height_bmi_6_pop.pdf",  width=12, height=8)
# # p
# # dev.off()

# ### get table statistics

# get_stat_pop_trait <- function(df, pop, trait_list = c("HEIGHT", "BMI", "Type_2_diabetes")){
#   #pop = "EAS"
#   #trait = "HEIGHT"
#   res = list()
#   for(trait in trait_list){
#       pick = df$trait == trait & df$pop == pop
#       tem = df[pick,]
#       tem$mean = round(tem$mean, 4)
#       tem$sd = round(tem$sd, 4)
#       tem$mean_sd = paste(tem$mean, tem$sd, sep = "±")
#       pick_col = c("short_name", "mean_sd")
#       tem = tem[,pick_col]
#       res[[trait]] = tem
#   }
#   res_df = as.data.frame(do.call(cbind, res))
#   pick = seq(1,ncol(res_df), by = 2)[-1]
#   res_df = res_df[,-pick]
#   names(res_df)[1] = "short_name"
#   return(res_df)
# }

# res_all = list()
# pop_list = c("EAS", "EUR", "SAS", "AMR", "AFR", "VNP")
# for(pop in pop_list){
#   res_all[[pop]] = get_stat_pop_trait(df, pop) 
# }

# # for sorting

# array_size = fread("../array_size.csv")

# export_pop <- function(pop_list, res_all, array_size){
#   array_size$size = round(array_size$No.Assays / 1000, 0)
#   array_size$array = array_size$'Short name'

#   require(data.table)
#   #pop_list = c("AFR", "AMR")
#   tem = res_all[pop_list]
#   res = merge(res_all[[pop_list[1]]], res_all[[pop_list[2]]], by = "short_name")
#   res2 = res 
#   res2$size = array_size$No.Assays[match(res2$short_name, array_size$array)]
#   od = order(res2$size)
#   res = res[od,]
#   out_name = paste(pop_list, collapse = "_")
#   out_name = paste0(out_name, "_PRS_cor_12_cutoff.csv")
#   fwrite(res, out_name, sep = ",")
# }


# export_pop(c("AFR", "AMR"), res_all, array_size)
# export_pop(c("EAS", "EUR"), res_all, array_size)
# export_pop(c("SAS", "VNP"), res_all, array_size)