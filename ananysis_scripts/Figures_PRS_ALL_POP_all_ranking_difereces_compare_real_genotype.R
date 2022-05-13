
get_data_scatter_array <- function(pheno, array, pop, trait){
  #wgs_path = paste0("WGS_", pop, "_", trait, ".all_score" )
  #wgs = fread()
  
  array_path = paste0(array, "_", pop, "_", trait, ".all_score")
  arr = fread(array_path)
  pick = intersect(arr[[1]], pheno[[1]])
  arr = arr[arr[[1]] %in% pick,]
  pheno = pheno[pheno[[1]] %in% pick,]
  cols = names(arr)[-c(1,2)]
  #pheno = pheno[, ..cols]
  arr2 = arr[,..cols]

  x = pheno[[3]] #rowMeans(wgs2)
  y = arr2[[3]] #rowMeans(arr2)

  x1 = ecdf(x)
  y1 = ecdf(y)

  x = x1(x)
  y = y1(y)
  pct_dif = sum(abs(x - y))

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

array_list = c("WGS", "Axiom_GW_ASI", "Axiom_GW_CHB", "Axiom_GW_EUR", "Axiom_GW_PanAFR", "Axiom_JAPONICA", "infinium-omnizhonghua-v1.4", "infinium-psycharray-v1.3", "japanese-screening-array-v1.0", "multi-ethnic-eur-eas-sas-v1.0", "multi-ethnic-global-v1.0", "oncoarray-500k", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB_WCSG", "chinese-genotyping-array-v1.0", "cytosnp-850k-v1.2", "global-screening-array-v.3", "human-cytosnp-12-v2.1", "infinium-core-v1.2", "infinium-global-diversity-array-v1.0", "infinium-omni2.5.v1.5", "infinium-omni5-v1.2", "GenomeWideSNP_6.0")





# "Breast_cancer"
get_data_pop_trait <- function(pheno_path, pop, array_list, trait){

  #wgs_path = paste0("WGS_", pop, "_", trait, ".all_score" )
  pheno = fread(pheno_path)
  res = list()

  for(array in array_list){
      res[[array]] = get_data_scatter_array(pheno, array, pop, trait)
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






df = get_data_pop_trait(pheno_path = "/media/datn/data2gb/GitHub/SNP_array_comparsion/data/Phenotype/Height.txt",  "VNP", array_list, "HEIGHT")
df = rename_array(df)
 ggplot(data = df, aes(pct_dif)) + geom_histogram() + facet_wrap(~ array_rename, nrow = 4)  + theme_light() + ylab("Count") + xlab("Percentile change range")
