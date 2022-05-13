setwd("/media/datn/data/1st_DC_PRS_array_project/array_annotation/annotation_all_hg38")

require(data.table)

array_list = c("Axiom_GW_ASI", "Axiom_GW_CHB", "Axiom_GW_EUR", "Axiom_GW_PanAFR", "Axiom_UKB_WCSG", "GenomeWideSNP_6.0", "Axiom_PMRA", "Axiom_JAPONICA", "Axiom_PMDA", "chinese-genotyping-array-v1.0", "infinium-core-v1.2", "infinium-omni2.5.v1.5", "infinium-omnizhonghua-v1.4", "infinium-psycharray-v1.3", "multi-ethnic-eur-eas-sas-v1.0", "oncoarray-500k", "cytosnp-850k-v1.2", "global-screening-array-v.3", "human-cytosnp-12-v2.1", "infinium-omni5-v1.2", "japanese-screening-array-v1.0", "multi-ethnic-global-v1.0", "infinium-global-diversity-array-v1.0")


array_size = data.frame(array = array_list)
array_size$all_site = 0
array_size$unique_site = 0

# array = array_list[1]
get_stat_array <- function(array){
  autosomes = as.character(c(1:22))
  x = "X"
  y = "Y"
  mt = "MT"
  file = fread(paste0(array, "_hg38.bed"))
  n_assay = nrow(file)
  file$id = paste(file$V1, file$V2, sep = ":")
  dup = duplicated(file$id)
  unique_pos = file[!dup,]
  n_position = length(unique_pos$id)
  n_autosomes = sum(unique_pos$V1 %in% autosomes)
  n_x = sum(unique_pos$V1 %in% x)
  n_y = sum(unique_pos$V1 %in% y)
  n_mt = sum(unique_pos$V1 %in% mt)
  # n_autosomes + n_x + n_y + n_mt
  res = data.frame(No.Assays = n_assay, No.Positions = n_position, No.Autosomal = n_autosomes, No.X = n_x, No.y = n_y, No.MT = n_mt)
}

res_all = list()
for(array in array_list){
  res_all[[array]] = get_stat_array(array)
}

df = do.call("rbind", res_all)
array = rownames(df)
df = cbind(array, df)
fwrite(df, "array_stat.txt", sep = "\t")

df = fread("/media/datn/data/1st_DC_PRS_array_project/array_annotation/array_stat_with_name.csv")
od = order(df$No.Assays)
df = df[od, ]
df = df[,-1]
fwrite(df, "array_stat_ordered_with_name.tsv", sep = "\t")


