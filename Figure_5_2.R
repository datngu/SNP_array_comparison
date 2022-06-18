
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
write.csv(stat_df, file = "../../output_paper/tables/Imputation_acc_real_simulated_PMRA.csv", quote = F, row.names = F)




# # p2



cutoffs = list("(0-0.01]" = c(0, 0.01), "(0.01-0.05]" = c(0.01, 0.05), "(0.01-0.5]" = c(0.01, 0.5),  "(0.05-0.5]" = c(0.05,0.5))


get_mean_correlation <- function(df, cutoffs, flag = "flag"){

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




# get_cor_array <- function(wgs, array, pop , trait){
#   #wgs_path = paste0("WGS_", pop, "_", trait, ".all_score" )
#   #wgs = fread()
#   cols = names(wgs)[-c(1,2)]
#   array_path = paste0("PRS_result/", array, "_", pop, "_", trait, ".all_score")
#   arr = fread(array_path)
#   res = c()
#   for(col in cols){
#     tem = cor(wgs[,..col], arr[,..col])
#     res = c(res, tem)
#   }
#  # res2 = c(array, trait, res, mean(res), sd(res))
#  # names(res2) = c("array", "trait", cols, "mean", "sd")
#   res2 = c(res)
#   names(res2) = c(cols)
#   return(res2)
# }



# get_cor_pop <- function(pop, array_list, trait_list = c("HEIGHT", "BMI", "Type_2_diabetes")){

#   res = list()

#   for( trait in trait_list){
#     wgs_path = paste0("PRS_result/95_shared_sample_wgs_", pop, "_", trait, ".all_score" )
#     wgs = fread(wgs_path)
#     all_cor = list()
#     for( array in array_list){
#     tem = get_cor_array(wgs, array, pop, trait)
#     all_cor[[array]] = tem
#     }
#     res[[trait]] = as.data.frame(do.call(cbind, all_cor))
#     res[[trait]]$trait = trait
#   }
#   df = as.data.frame(do.call(rbind, res))
#   df$pop = pop
#   return(df)
# }


# df = get_cor_pop("VNP", c("95_shared_sample_simulated", "95_shared_sample_pmra"))
# names(df) = c("sim" , "real",      "trait"  ,                    "pop")

# df$trait[df$trait == "HEIGHT"] = "Height"

# p2 =  ggplot(df,aes(x = real, y = sim, color = trait)) +  geom_point() + theme_light() + ylab("Simulated PMRA") + xlab("Real PMRA")  + geom_smooth( aes(color = NULL), method=lm, se = F) + scale_y_continuous(breaks=seq(0,1,0.01), limits = c(0.93,0.985) ) + scale_x_continuous(breaks=seq(0,1,0.01), limits = c(0.93,0.985) ) 


# # regEq <- function(lmObj, dig) {
# #     paste0("y = ",
# #         paste0(
# #             c(round(lmObj$coef[1], dig), round(sign(lmObj$coef[-1])*lmObj$coef[-1], dig)),
# #             c("", rep("*", length(lmObj$coef)-1)),
# #             paste0(c("", names(lmObj$coef)[-1]), c(ifelse(sign(lmObj$coef)[-1]==1," + "," - "), "")),
# #             collapse=""
# #         )
# #     )
# # }

# # m = lm(sim ~ real, data = df)
# # eq = regEq(m,3)

# t = paste0("R2: ", round( (cor(df[,1], df[,2]))^2,4))
# p2 <- p2 + annotate(geom = "text", x=0.97, y=0.94, label= t)



# pdf(file= "../Figure_6.pdf",  width=8, height=6)
# p2
# dev.off()


