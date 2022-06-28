
setwd("/Users/datn/github/SNP_array_comparison")
setwd("data/real_data")

require(scales)
require(ggplot2)
require(data.table)


# p1: imputation
cutoffs = list("(0-0.01]" = c(0, 0.01), "(0.01-0.05]" = c(0.01, 0.05), "(0.01-0.5]" = c(0.01, 0.5),  "(0.05-0.5]" = c(0.05,0.5))


get_mean_imp_acc <- function(df, cutoffs, flag = "flag"){

  res = list()
  for( i in 1: length(cutoffs)){
    cut = cutoffs[[i]]
    cut_name = names(cutoffs)[i]
    df2 = df[df$AF > cut[1] & df$AF <= cut[2],]
    m = mean(df2$r_2, na.rm =T)
    t = paste(cut_name, flag, sep = "_")
    r = c(m,t, flag)
    res[[t]] = r  
  }
  rdf = do.call("rbind", res)
  return(rdf)
}


get_imp_acc_array <- function(array_path, flag){
  res = list()
  cutoffs = list("(0-0.01]" = c(0, 0.01), "(0.01-0.05]" = c(0.01, 0.05), "(0.01-0.5]" = c(0.01, 0.5),  "(0.05-0.5]" = c(0.05,0.5))
  for(i in 1:22){
    x = fread( paste0(array_path, "/chr", i , "_imputed.correlation.txt.gz"))
    res[[i]] =  get_mean_imp_acc(x, cutoffs, flag)
  }
  df = do.call("rbind", res)
  df = as.data.frame(df)
  names(df) = c("imp_acc", "Bin", "Type")
  df$imp_acc = as.numeric(df$imp_acc)
  return(df)
}



# res = list()

# for(i in 1:22){
#   # real

#   x = fread( paste0("95_shared_sample_pmra/chr", i , "_imputed.correlation.txt.gz"))
#   res[[i]] =  get_mean_imp_acc(x, cutoffs, "Real_PMRA")

# }

# df1 = do.call("rbind", res)



# res = list()

# for(i in 1:22){
#   # simulated

#   x = fread( paste0("95_shared_sample_simulated/chr", i , "_imputed.correlation.txt.gz"))
#   res[[i]] =  get_mean_imp_acc(x, cutoffs, "Simulated_PMRA")
# }

# df2 = do.call("rbind", res)

df1 = get_imp_acc_array("pmra_real_array", "PMRA_Real")
df2 = get_imp_acc_array("pmra_simulated_array", "PMRA_Simulated")

df3 = get_imp_acc_array("gsa_real_array", "GSA_Real")
df4 = get_imp_acc_array("gsa_simulated_array", "GSA_Simulated")

df = rbind(df1, df2, df3, df4)

# names(df) = c("imp_acc", "Bin", "Type")

# df$imp_acc = as.numeric(df$imp_acc)

p1 <- ggplot(df, aes( y = imp_acc, x = Bin, color=  Type)) + geom_boxplot() + theme_light() + ylab("Mean imputation r2") + xlab("") + guides(x = guide_axis(angle = 90)) 


pdf(file= "../../output_paper/Figure_5.pdf",  width=8, height=6)
p1
dev.off()



bin = unique(df$Bin)
df$imp_acc = as.numeric(df$imp_acc)
mean_imp_acc = c()
sd_imp_acc = c()
for(b in bin){
  t1 = mean(df$imp_acc[df$Bin == b])
  t2 = sd(df$imp_acc[df$Bin == b])
  mean_imp_acc = c(mean_imp_acc, t1)
  sd_imp_acc = c(sd_imp_acc, t2)
}

mean_imp_acc = round(mean_imp_acc,4)
sd_imp_acc = round(sd_imp_acc,4)
mean_sd = paste(mean_imp_acc, sd_imp_acc, sep = "±")
stat_df = data.frame(Bin = bin, Imputation_accuracy = mean_sd)
write.csv(stat_df, file = "../../output_paper/tables/Imputation_acc_real_simulated_pmra_gsa.csv", quote = F, row.names = F)






