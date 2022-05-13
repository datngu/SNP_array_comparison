setwd("/Users/datn/github/SNP_array_comparison")

setwd("data/PRS_scores")



# get /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS/PRS_WGS/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores
# get /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008/PRS_result/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores

# get /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR/PRS_WGS/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores
# get /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed/PRS_result_all_pop/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores

require(data.table)
require(scales)
require(ggplot2)
# wgs_height =  wgs =  fread("WGS_VNP_HEIGHT.all_score")
# wgs_bmi = fread("WGS_VNP_BMI.all_score")
# arr = fread("Axiom_JAPONICA_HEIGHT.all_score")
# arr = fread("global-screening-array-v.3_HEIGHT.all_score")
# col = names(wgs)[6]
# cor(wgs[,..col], arr[,..col])



get_percentile_dif <- function(wgs, array, pop, trait, cutoff ){
  #wgs_path = paste0("WGS_", pop, "_", trait, ".all_score" )
  #wgs = fread()
  cols = names(wgs)[-c(1,2)]
  array_path = paste0(array, "_", pop, "_", trait, ".all_score")
  arr = fread(array_path)
  wgs2 = wgs[, ..cols]
  arr2 = arr[,..cols]

  x = wgs2[[cutoff]] #rowMeans(wgs2)
  y = arr2[[cutoff]] #rowMeans(arr2)

  x1 = ecdf(x)
  y1 = ecdf(y)

  x = x1(x)
  y = y1(y)
  pct_dif = abs(x - y)*100

  arr_percentile <- ecdf(arr[[3]])
  y = arr_percentile(arr[[3]])

  tem = data.frame(wgs_pct = x, arr_pct = y, pct_dif = pct_dif, trait = trait)
  tem$array = array
  return(tem)
}




get_percentile_trait <- function(pop, array_list, trait, cutoff){

  wgs_path = paste0("WGS_", pop, "_", trait, ".all_score" )
  wgs = fread(wgs_path)
  res = list()

  for(array in array_list){
      res[[array]] = get_percentile_dif(wgs, array, pop, trait, cutoff)
  }
  df = as.data.frame(do.call(rbind, res))
  df$pop = pop
  return(df)
}

get_percentile_pop <- function(pop, array_list, trait_list, cutoff){
  res = list()
  for(trait in trait_list){
      res[[trait]] = get_percentile_trait(pop, array_list, trait, cutoff)
  }
  df = as.data.frame(do.call(rbind, res))
  df$pop = pop
  return(df)
}

rename_array <- function(df, path_csv = "../array_size.csv"){ 
  array_size = fread(path_csv)
  #
  array_size$size = round(array_size$No.Assays / 1000, 0)  
  array_size$array = array_size$'Short name'  
  array_size$array_rename = paste0(array_size$array, "_", array_size$size, "k")  
  df$short_name = array_size$array[match(df$array, array_size$old_name)]  
  df$array_rename = array_size$array_rename[match(df$array, array_size$old_name)]  
  #df$array <- factor(df$array, levels = array_size$array)
  df$array_rename <- factor(df$array_rename, levels = array_size$array_rename)
  return(df)
} 


plot_6_pop <- function(pop_list, array_list, cutoff){
  p = list()
  for(pop in pop_list){
    df = get_percentile_pop(pop, array_list, trait_list, cutoff)
    df = rename_array(df)
    df$trait[df$trait == "HEIGHT"] = "Height" 
    p[[pop]] = ggplot(data=df, aes(x = array_rename, y = pct_dif, fill = array_rename)) + geom_boxplot(outlier.shape = NA) + theme_light() + ylab("Percentile difference") + xlab("SNP array")  + facet_wrap(~ trait, nrow = 1) + theme(legend.position="none", axis.title.x = element_blank()) + guides(x = guide_axis(angle = 90)) + ggtitle(pop) +  scale_y_continuous(breaks=seq(0,100,5), limits = c(0,30) )
  }
    #p[[1]] = p[[1]] + theme(legend.position="none") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) 
    p[[1]] = p[[1]] + theme(legend.position="none") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
    p[[3]] = p[[3]] + theme(legend.position="none") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
    p[[2]] = p[[2]] + theme(legend.position="none") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
    p[[4]] = p[[4]] + theme(legend.position="none") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
    p[[5]] = p[[5]] + theme(legend.position="none") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
    p[[6]] = p[[6]] + theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())
    library(patchwork)
    m = wrap_plots(p, nrow = 6, plot_layout = "auto")
    return(m)
}

# cutoff = "Pt_0.001"

# m = plot_6_pop(pop_list, array_list, cutoff)
# pdf(file= "../../PCT_ABS_dif_all.pdf",  width=12, height=17)
# print(m)
# dev.off()


pop_list = c("AFR", "AMR", "EAS", "EUR", "SAS", "VNP")
#c("EAS", "EUR", "SAS", "AMR", "AFR", "VNP")
trait_list = c("HEIGHT", "BMI", "Type_2_diabetes")

array_list = c("Axiom_GW_ASI", "Axiom_GW_CHB", "Axiom_GW_EUR", "Axiom_GW_PanAFR", "Axiom_JAPONICA", "infinium-omnizhonghua-v1.4", "infinium-psycharray-v1.3", "japanese-screening-array-v1.0", "multi-ethnic-eur-eas-sas-v1.0", "multi-ethnic-global-v1.0", "oncoarray-500k", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB_WCSG", "chinese-genotyping-array-v1.0", "cytosnp-850k-v1.2", "global-screening-array-v.3", "human-cytosnp-12-v2.1", "infinium-core-v1.2", "infinium-global-diversity-array-v1.0", "infinium-omni2.5.v1.5", "infinium-omni5-v1.2", "GenomeWideSNP_6.0")



#cutoff_list = c("Pt_0.001", "Pt_0.05", "Pt_0.1", "Pt_0.15", "Pt_0.2", "Pt_0.25", "Pt_0.3", "Pt_0.35", "Pt_0.4", "Pt_0.45", "Pt_0.5", "Pt_0.6", "Pt_0.7", "Pt_0.8", "Pt_0.9", "Pt_1")

cutoff_list = c("Pt_5e-08", "Pt_1e-07", "Pt_1e-06", "Pt_1e-05", "Pt_0.0001", "Pt_0.001", "Pt_0.01", "Pt_0.1", "Pt_0.2", "Pt_0.3", "Pt_0.5", "Pt_1")



# plotting main fig 4
cutoff = cutoff_list[1]
out_fn = paste0("../../output_paper/Figure_4.pdf")
m = plot_6_pop(pop_list, array_list, cutoff)
pdf(file= out_fn,  width=12, height=15)
print(m)
dev.off()


# ploting supp figs
for(cutoff in cutoff_list[-1]){
  out_fn = paste0("../../output_paper/PCT_ABS_dif_6_pop_3_trait_" , cutoff, ".pdf")
  m = plot_6_pop(pop_list, array_list, cutoff)
  pdf(file= out_fn,  width=12, height=15)
  print(m)
  dev.off()
}





## Get statistics

for(cutoff in cutoff_list){

  res_all = list()  

  for(pop in pop_list){
    res_all[[pop]] = get_percentile_pop(pop, array_list, trait_list= trait_list, cutoff = cutoff)
  } 

  df = as.data.frame(do.call(rbind, res_all))
  df = rename_array(df)
  df$trait[df$trait == "HEIGHT"] = "Height" 
  

  ###################
  # computing statistics  

  trait_l = unique(df$trait)
  pop_l = unique(df$pop)
  tem = as.character(levels(df$array_rename))
  array_l = c()
  for(i in tem){
    array_l = c(array_l, df$short_name[ df$array_rename == i][1])
  } 

  l_sta = list()  

  for(t in trait_l){
    df1 = df[df$trait == t,]
    sta = matrix(0, length(array_l), length(pop_l) )
    for (i in 1:length(array_l)){
      for(j in 1:length(pop_l)){
        pick = df1$short_name == array_l[i] & df1$pop == pop_l[j]
        sta[i,j] = mean(df1$pct_dif[pick])
      }
    }
    sta = as.data.frame(sta)
    colnames(sta) = pop_l
    sta_array = data.frame(Array_name = array_l, Trait = t)
    res = cbind(sta_array, sta)
    l_sta[[t]] = res    
  }
  out_tab = paste0("../../output_paper/tables/PCT_ABS_dif_" , cutoff, ".csv")
  stat_table = as.data.frame(do.call(rbind, l_sta))
  write.table(stat_table, sep = ",", row.names = F, col.names = T, file = out_tab)
}








