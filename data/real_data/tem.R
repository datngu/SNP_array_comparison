setwd("/Users/datn/github/SNP_array_comparison/data/real_data/PRS_result")

get_prs_correlation_array <- function(array_name, wgs_name, flag){
  
  array_name = "gsa_real_array"
  wgs_name = "wgs_24_sample"
  trait = "HEIGHT"
  flag = "GSA_Real"
  
  trait_list = c("HEIGHT", "BMI", "Type_2_diabetes")
  res_trait = list()
  for(trait in trait_list){
    arr_prs_fn = paste0(array_name, "_VNP_", trait, ".all_score")
    wgs_prs_fn = paste0(wgs_name, "_VNP_", trait, ".all_score")
    
    arr = fread(arr_prs_fn)
    wgs = fread(wgs_prs_fn)
    arr = as.data.frame(arr)
    wgs = as.data.frame(wgs)
    
    cutoffs = names(arr)[-c(1:2)]
    res = list()
    for(c in cutoffs){
      res[[c]] =  cor(arr[,c], wgs[,c], method = "pearson")
    }
    df = do.call("rbind", res)
    df = as.data.frame(df)
    df$trait = trait
    res_trait[[trait]] = df
  }
  
  
  df = do.call("rbind", res_trait)
  df = as.data.frame(df)
  df$flag = flag
  names(df) = c("prs_cor", "trait", "Type")
  return(df)
}



p1 <- ggplot(df, aes( y = prs_cor, x = trait, color=  Type)) + geom_boxplot() + theme_light() + ylab("Mean PGS correlation") + xlab("") + guides(x = guide_axis(angle = 90)) 

