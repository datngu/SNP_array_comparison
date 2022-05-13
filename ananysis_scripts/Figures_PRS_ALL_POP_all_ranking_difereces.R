# get /dragennfs/area8/datnguyen/PRS_arrays_project/KHV_10_folds_imputed/PRS_result/*all_score /media/datn/data/1st_DC_PRS_array_project/PRS_codes/PRS_khv

setwd("/media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores")


# get /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS/PRS_WGS/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores
# get /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008/PRS_result/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores

# get /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR/PRS_WGS/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores
# get /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed/PRS_result_all_pop/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores

require(data.table)
require(scales)
# wgs_height =  wgs =  fread("WGS_VNP_HEIGHT.all_score")
# wgs_bmi = fread("WGS_VNP_BMI.all_score")
# arr = fread("Axiom_JAPONICA_HEIGHT.all_score")
# arr = fread("global-screening-array-v.3_HEIGHT.all_score")
# col = names(wgs)[6]
# cor(wgs[,..col], arr[,..col])



get_data_scatter_array <- function(wgs, array, pop, trait){
  #wgs_path = paste0("WGS_", pop, "_", trait, ".all_score" )
  #wgs = fread()
  cols = names(wgs)[-c(1,2)]
  array_path = paste0(array, "_", pop, "_", trait, ".all_score")
  arr = fread(array_path)
  wgs2 = wgs[, ..cols]
  arr2 = arr[,..cols]

  x = wgs2[[3]] #rowMeans(wgs2)
  y = arr2[[3]] #rowMeans(arr2)

  x1 = ecdf(x)
  y1 = ecdf(y)

  x = x1(x)
  y = y1(y)
  pct_dif = abs(x - y)

  arr_percentile <- ecdf(arr[[3]])
  y = arr_percentile(arr[[3]])

  tem = data.frame(wgs_pct = x, arr_pct = y, pct_dif = pct_dif)
  tem$array = array
  return(tem)
}

#df1 = get_data_scatter_array(wgs, array, pop, trait)
#df2 = get_data_scatter_array(wgs, "Axiom_JAPONICA", pop, trait)
#df3 = rbind(df1, df2)


#library(ggplot2)
#ggplot(df3, aes(x= wgs_pct, y= arr_pct)) + geom_point() + geom_abline(intercept = 0, slope = 1) + facet_wrap(~ array, nrow = 2)

#ggplot(data = df3, aes(pct_dif)) + geom_histogram() + facet_wrap(~ array, nrow = 2)

#get_cor_array(wgs, "Axiom_JAPONICA", "VNP", "HEIGHT")

array_list = c("Axiom_GW_ASI", "Axiom_GW_CHB", "Axiom_GW_EUR", "Axiom_GW_PanAFR", "Axiom_JAPONICA", "infinium-omnizhonghua-v1.4", "infinium-psycharray-v1.3", "japanese-screening-array-v1.0", "multi-ethnic-eur-eas-sas-v1.0", "multi-ethnic-global-v1.0", "oncoarray-500k", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB_WCSG", "chinese-genotyping-array-v1.0", "cytosnp-850k-v1.2", "global-screening-array-v.3", "human-cytosnp-12-v2.1", "infinium-core-v1.2", "infinium-global-diversity-array-v1.0", "infinium-omni2.5.v1.5", "infinium-omni5-v1.2", "GenomeWideSNP_6.0")





# "Breast_cancer"
get_data_pop_trait <- function(pop, array_list, trait){

  wgs_path = paste0("WGS_", pop, "_", trait, ".all_score" )
  wgs = fread(wgs_path)
  res = list()

  for(array in array_list){
      res[[array]] = get_data_scatter_array(wgs, array, pop, trait)
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






df = get_data_pop_trait("VNP", array_list, "HEIGHT")
df = rename_array(df)
 ggplot(data = df, aes(pct_dif)) + geom_histogram() + facet_wrap(~ array_rename, nrow = 4)  + theme_light() + ylab("Count") + xlab("Percentile change range")


ggplot(data = df, aes(pct_dif)) + geom_histogram(bins=15) + facet_wrap(~ c, nrow = 4)  + theme_light() + ylab("Count") + xlab("Squere percentile change range")

# geom_histogram(bins=30)

ggplot(df, aes(x=array_rename, y=pct_dif)) + 
  geom_boxplot() + guides(x = guide_axis(angle = 90))


#ggplot(VNP_HEIGHT, aes(x= wgs_pct, y= arr_pct)) + stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) + geom_point(size=0.5) + stat_density2d(geom="tile", aes(fill=..density..^0.25,     alpha=ifelse(..density..^0.25<0.4,0,1)), contour=FALSE) +  scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) + geom_abline(intercept = 0, slope = 1) + facet_wrap(~ array_rename, nrow = 4)

#ggplot(VNP_HEIGHT, aes(x= wgs_pct, y= arr_pct)) + geom_point() + geom_abline(intercept = 0, slope = 1) + facet_wrap(~ array_rename, nrow = 4)

#ggplot(data = VNP_HEIGHT, aes(pct_dif)) + geom_histogram() + facet_wrap(~ array_rename, nrow = 4)
#out_name = paste0(pop, "_", trait, "_23_array_percentile_dif.pdf")
out_name ="all_pop_23_array_percentile_dif.pdf"
pop_list = c("EAS", "EUR", "SAS", "AMR", "AFR", "VNP")

for(pop in pop_list){
 trait = "HEIGHT"
 df = get_data_pop_trait(pop, array_list, trait)
 df = rename_array(df)
 p = ggplot(data = df, aes(pct_dif)) + geom_histogram() + facet_wrap(~ array_rename, nrow = 4)  + theme_light() + ylab("Count") + xlab("Percentile change range") 
 out_name = paste0(pop, "_", trait, "_23_array_percentile_dif.pdf")
 pdf(file = out_name,  width=12, height=9)
 p
 dev.off()
}



plot_pop_trait <- function(pop, trait, array_list){
   df = get_data_pop_trait(pop, array_list, trait)
   df = rename_array(df)
   p = ggplot(data = df, aes(pct_dif)) + geom_histogram() + facet_wrap(~ array_rename, nrow = 4)  + theme_light() + ylab("Count") + xlab("Percentile change range") 
   out_name = paste0(pop, "_", trait, "_23_array_percentile_dif.pdf")
   pdf(file = out_name,  width=12, height=9)
   p
   dev.off()
}
