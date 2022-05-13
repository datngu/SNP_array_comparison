
require(data.table)
require(ggplot2)

setwd("/media/datn/data2gb/GitHub/SNP_array_comparsion/data/Imputation_correlation")

array_list = c("Axiom_GW_ASI", "Axiom_GW_CHB", "Axiom_GW_EUR", "Axiom_GW_PanAFR", "Axiom_JAPONICA", "infinium-omnizhonghua-v1.4", "infinium-psycharray-v1.3", "japanese-screening-array-v1.0", "multi-ethnic-eur-eas-sas-v1.0", "multi-ethnic-global-v1.0", "oncoarray-500k", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB_WCSG", "chinese-genotyping-array-v1.0", "cytosnp-850k-v1.2", "global-screening-array-v.3", "human-cytosnp-12-v2.1", "infinium-core-v1.2", "infinium-global-diversity-array-v1.0", "infinium-omni2.5.v1.5", "infinium-omni5-v1.2", "GenomeWideSNP_6.0")

# excluded: "infinium-exome-v1.0"

#c("(0-0.01]", "(0.01-0.05]", "(0.01-0.5]", "(0.05-0.5]")
get_acc_pop <- function(pop, array_list, cutoffs = list("(0-0.01]" = c(0, 0.01), "(0.01-0.05]" = c(0.01, 0.05), "(0.01-0.5]" = c(0.01, 0.5),  "(0.05-0.5]" = c(0.05,0.5))){
  require(data.table)
  # pop = "EAS"
  # array_list = c("Axiom_GW_ASI", "Axiom_GW_CHB", "Axiom_GW_EUR", "Axiom_GW_PanAFR", "Axiom_JAPONICA", "infinium-omnizhonghua-v1.4", "infinium-psycharray-v1.3", "japanese-screening-array-v1.0", "multi-ethnic-eur-eas-sas-v1.0", "multi-ethnic-global-v1.0", "oncoarray-500k", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB_WCSG", "chinese-genotyping-array-v1.0", "cytosnp-850k-v1.2", "global-screening-array-v.3", "human-cytosnp-12-v2.1", "infinium-core-v1.2", "infinium-global-diversity-array-v1.0", "infinium-omni2.5.v1.5", "infinium-omni5-v1.2", "GenomeWideSNP_6.0")
  # array = array_list[1]
  # cutoffs = list("(AC>1, MAF<=0.01]" = c(0, 0.01), "(MAF>0.01, MAF<=0.05]" = c(0.01, 0.05))
  #CREATE list to save values
 
  chr_list = as.character(c(1:22))
  tem_vec = rep(0.0, 22)
  names(tem_vec) = chr_list
  tem_list = list()
  sizes = rep(0, length(cutoffs))
  names(sizes) = names(cutoffs)
  for(cut in names(cutoffs)){
    tem_list[[cut]] = tem_vec
  }
  
  res_list = list()
  array_statistis = list(imp_acc = tem_list, imp_cov = tem_list)

  for(array in array_list){
    res_list[[array]] = array_statistis
  }

  sizes = rep(0, length(cutoffs))
  names(sizes) = names(cutoffs)

  for( chr in chr_list){
    path = paste0("chr", chr, "_", pop, "_merged_23_array_correlation.txt.gz")
    tem = fread(path)
    tem = tem[tem$AC >1,]

    for(cut in names(cutoffs)){
        x = cutoffs[[cut]]
        pick = tem$MAF > x[1] & tem$MAF <= x[2]
        sizes[cut] = sizes[cut] + sum(pick)
    }

    for(array in array_list){
      pick_col = c(colnames(tem)[1:4], array)
      tem2 = tem[,..pick_col]
      tem2$r_2 = tem2[[array]]^2
      for(cut in names(cutoffs)){
        x = cutoffs[[cut]]
        pick = tem2$MAF > x[1] & tem2$MAF <= x[2]
        acc = mean(tem2$r_2[pick], na.rm = T)
        cov = sum(tem2$r_2[pick] >= 0.8, na.rm = T) / sum(pick)
        res_list[[array]][["imp_acc"]][[cut]][chr] = acc
        res_list[[array]][["imp_cov"]][[cut]][chr] = cov
      }
    }
  }

  res = list()
  for(array in array_list){
    for(cut in names(cutoffs)){
      mean_r2 = mean(res_list[[array]][["imp_acc"]][[cut]])
      sd_r2 = sd(res_list[[array]][["imp_acc"]][[cut]])
      mean_cov = mean(res_list[[array]][["imp_cov"]][[cut]])
      sd_cov = sd(res_list[[array]][["imp_cov"]][[cut]])
      res_df = data.frame(array = array, number_variant = sizes[cut], mean_r2 = mean_r2, sd_r2 = sd_r2, mean_cov = mean_cov, sd_cov = sd_cov, MAF_range = cut, pop = pop)
      res = c(res, list(res_df))
    }
  }
  df = do.call("rbind", res)
  return(df)
}


#t = get_acc_pop("EAS", array_list)


pop_list = c("EAS", "EUR", "SAS", "AMR", "AFR", "VNP")

res_pop = list()


for( pop in pop_list){
    tem = get_acc_pop(array_list, pop = pop)
    #tem$pop = pop
    # rename KHV to VNP
    #if(pop == "KHV") tem$pop = "VNP"
    res_pop[[pop]] = tem
    #res_pop[[pop]] = get_stat_pop(array_list2, pop = pop)
}

df = as.data.frame(do.call(rbind, res_pop))
#df$pop = df[,7]
#fwrite(df, file = "cor_data_all_pop_rename.csv", sep = "\t")
df = fread("cor_data_all_pop_rename.csv")

array_size = fread("../array_size.csv")


#
array_size$size = round(array_size$No.Assays / 1000, 0)

array_size$array = array_size$'Short name'

array_size$array_rename = paste0(array_size$array, "_", array_size$size, "k")

df$array_rename = array_size$array_rename[match(df$array, array_size$old_name)]

df$short_name = array_size$array[match(df$array, array_size$old_name)]

#df$array <- factor(df$array, levels = array_size$array)
df$array_rename <- factor(df$array_rename, levels = array_size$array_rename)


# ggplot(data=df, aes(x = array_rename, y = mean, fill = trait, color = trait)) + geom_point(stat="identity", position = position_dodge(.5)) + geom_errorbar(aes(x = array_rename, ymin = mean - sd, ymax = mean + sd), width = 0, position = position_dodge(.5)) + guides(x = guide_axis(angle = 90)) + scale_y_continuous(breaks=seq(0,1,0.01), limits = c(0.8,1) ) + theme_light() + ylab("Mean correlation to WGS PRS") + xlab("SNP array")  + facet_wrap(~ pop, nrow = 2)



require(ggplot2)
p1 = ggplot(data=df, aes(x = array_rename, y = mean_r2, fill = MAF_range, color = MAF_range)) +
geom_point(stat="identity", position = position_dodge(.5)) + geom_errorbar(aes(x = array_rename, ymin = mean_r2 - sd_r2, ymax = mean_r2 + sd_r2), width = 0, position = position_dodge(.5)) + guides(x = guide_axis(angle = 90)) + scale_y_continuous(breaks=seq(0,1,0.05)) + theme_light() + ylab("Mean accuracy r2") + xlab("SNP array") + facet_wrap(~ pop, nrow = 2)


p2 = ggplot(data=df, aes(x = array_rename, y = mean_cov, fill = MAF_range, color = MAF_range)) +
geom_point(stat="identity", position = position_dodge(.5)) + geom_errorbar(aes(x = array_rename, ymin = mean_cov - sd_r2, ymax = mean_cov + sd_r2), width = 0, position = position_dodge(.5)) + guides(x = guide_axis(angle = 90)) + scale_y_continuous(breaks=seq(0,1,0.05)) + theme_light() + ylab("Mean coverage (r2s >= 0.8)") + xlab("SNP array") + facet_wrap(~ pop, nrow = 2)


#ggplot(data=df, aes(x = array, y = mean_cov, fill = MAF_range)) + geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(x = array, ymin = mean_cov - sd_r2, ymax = mean_cov + sd_r2), width = 0.5, position=position_dodge(.9)) + guides(x = guide_axis(angle = 45)) + scale_y_continuous(breaks=seq(0,1,0.05)) + theme_light() + ylab("Mean coverage (r2s >= 0.8)") + xlab("SNP array") + facet_wrap(~ pop, nrow = 2)

pdf(file= "Imputation_accuracy_coverage_all_chr_6-pops.pdf",   width=12, height=9)
p1
p2
dev.off()



# #####################################################################
### get table statistics

# list("(AC>1, MAF<=0.01]" = c(0, 0.01), "(MAF>0.01, MAF<=0.05]" = c(0.01, 0.05), "(MAF>0.01, MAF<=0.5]" = c(0.01, 0.5), "(MAF<0.05, MAF<=0.5]" = c(0.05,0.5)))

get_stat_pop_trait_acc <- function(df, pop, cutoff_list = c("(0-0.01]", "(0.01-0.05]", "(0.01-0.5]", "(0.05-0.5]")){
  # pop = "EAS"
  # trait = "HEIGHT"
  # cut = cutoff_list[1]
  df = as.data.frame(df)
  df$mean_r2 = as.numeric(df$mean_r2)
  df$sd_r2 = as.numeric(df$sd_r2)
  df$mean_cov = as.numeric(df$mean_cov)
  df$sd_cov = as.numeric(df$sd_cov)

  res = list()
  for(cut in cutoff_list){
      pick = df$MAF_range == cut & df$pop == pop
      tem = df[pick,]
      tem$mean_r2 = round(tem$mean_r2, 4)
      tem$sd_r2 = round(tem$sd_r2, 4)
      tem$mean_sd_r2 = paste(tem$mean_r2, tem$sd_r2, sep = "±")
      pick_col = c("short_name", "mean_sd_r2")
      tem = tem[,pick_col]
      res[[cut]] = tem
  }
  res_df = as.data.frame(do.call(cbind, res))
  pick = seq(1,ncol(res_df), by = 2)[-1]
  res_df = res_df[,-pick]
  names(res_df)[1] = "short_name"
  return(res_df)
}



res_all = list()
pop_list = c("EAS", "EUR", "SAS", "AMR", "AFR", "VNP")
for(pop in pop_list){
  res_all[[pop]] = get_stat_pop_trait_acc(df, pop) 
}

# for sorting

array_size = fread("../array_size.csv")

export_pop_imp_acc <- function(pop_list, res_all, array_size){
  array_size$size = round(array_size$No.Assays / 1000, 0)
  array_size$array = array_size$'Short name'

  require(data.table)
  #pop_list = c("AFR", "AMR")
  tem = res_all[pop_list]
  res = merge(res_all[[pop_list[1]]], res_all[[pop_list[2]]], by = "short_name")
  res2 = res 
  res2$size = array_size$No.Assays[match(res2$short_name, array_size$array)]
  od = order(res2$size)
  res = res[od,]
  out_name = paste(pop_list, collapse = "_")
  out_name = paste0(out_name, "_imp_acc.csv")
  fwrite(res, out_name, sep = ",")
}


export_pop_imp_acc(c("AFR", "AMR"), res_all, array_size)
export_pop_imp_acc(c("EAS", "EUR"), res_all, array_size)
export_pop_imp_acc(c("SAS", "VNP"), res_all, array_size)



###
# #####################################################################
### get table statistics

# list("(AC>1, MAF<=0.01]" = c(0, 0.01), "(MAF>0.01, MAF<=0.05]" = c(0.01, 0.05), "(MAF>0.01, MAF<=0.5]" = c(0.01, 0.5), "(MAF<0.05, MAF<=0.5]" = c(0.05,0.5)))

get_stat_pop_trait_cov <- function(df, pop, cutoff_list = c("(0-0.01]", "(0.01-0.05]", "(0.01-0.5]", "(0.05-0.5]")){
  # pop = "EAS"
  # trait = "HEIGHT"
  # cut = cutoff_list[1]
  df = as.data.frame(df)
  df$mean_r2 = as.numeric(df$mean_r2)
  df$sd_r2 = as.numeric(df$sd_r2)
  df$mean_cov = as.numeric(df$mean_cov)
  df$sd_cov = as.numeric(df$sd_cov)

  res = list()
  for(cut in cutoff_list){
      pick = df$MAF_range == cut & df$pop == pop
      tem = df[pick,]
      tem$mean_cov = round(tem$mean_cov, 4)
      tem$sd_cov = round(tem$sd_cov, 4)
      tem$mean_sd_cov = paste(tem$mean_cov, tem$sd_cov, sep = "±")
      pick_col = c("short_name", "mean_sd_cov")
      tem = tem[,pick_col]
      res[[cut]] = tem
  }
  res_df = as.data.frame(do.call(cbind, res))
  pick = seq(1,ncol(res_df), by = 2)[-1]
  res_df = res_df[,-pick]
  names(res_df)[1] = "short_name"
  return(res_df)
}



res_all = list()
pop_list = c("EAS", "EUR", "SAS", "AMR", "AFR", "VNP")
for(pop in pop_list){
  res_all[[pop]] = get_stat_pop_trait_cov(df, pop) 
}

# for sorting

array_size = fread("../array_size.csv")

export_pop_imp_cov <- function(pop_list, res_all, array_size){
  array_size$size = round(array_size$No.Assays / 1000, 0)
  array_size$array = array_size$'Short name'

  require(data.table)
  #pop_list = c("AFR", "AMR")
  tem = res_all[pop_list]
  res = merge(res_all[[pop_list[1]]], res_all[[pop_list[2]]], by = "short_name")
  res2 = res 
  res2$size = array_size$No.Assays[match(res2$short_name, array_size$array)]
  od = order(res2$size)
  res = res[od,]
  out_name = paste(pop_list, collapse = "_")
  out_name = paste0(out_name, "_imp_cov.csv")
  fwrite(res, out_name, sep = ",")
}


export_pop_imp_cov(c("AFR", "AMR"), res_all, array_size)
export_pop_imp_cov(c("EAS", "EUR"), res_all, array_size)
export_pop_imp_cov(c("SAS", "VNP"), res_all, array_size)
