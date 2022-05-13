#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
syntax='\nUsage:\t./this_sript.R true_genotype test_genotype out_file n_core'


if(length(args) == 0 ){
  cat("\nNo argument, Program stop! \n")
  cat(syntax)
  quit()
}

options(stringsAsFactors = FALSE)
require(data.table)
#require(foreach)
#require(doParallel)
#require(callr)
#registerDoParallel(cores=CPUNUM)

ref_path = args[1]
test_path = args[2]
out_file = args[3]
CPUNUM = args[4]

CPUNUM = as.numeric(CPUNUM)

#ref_path = "/dragennfs/area7/datnguyen/PRS_arrays_project/1000_KHV_process_GT/processed_GT_chr20_1013_KHV.txt.gz"
#test_path = "/dragennfs/area7/datnguyen/PRS_arrays_project/KHV_imputed/Axiom_GW_ASI/chr20.dose.txt.gz"
#ref_path = "/home/datn/Downloads/processed_GT_chr20_1013_KHV.txt.gz"
#test_path =  "/home/datn/Downloads/chr20.dose.txt.gz"

cat("\nreading inputs\n")
cat(ref_path)
ref = fread(ref_path, header = T)
cat("\nreading inputs\n")
cat(test_path)
test = fread(test_path, header = T)

id = intersect(ref$ID, test$ID)

ref = ref[ref$ID %in% id,]
test = test[test$ID %in% id, ]
test = test[order(match(test$ID, ref$ID)),]
test = test[,-"ID"]
ref_info = ref[, c(1:4)]
ref = ref[, -c(1:4)]

test = as.matrix(test)
ref = as.matrix(ref)


step = round(nrow(ref)/CPUNUM +1)
slipt_vec = c()
for( i in 1:CPUNUM){
  tem = rep(i, step)
  slipt_vec = c(slipt_vec, tem)
}

slipt_vec = slipt_vec[1:nrow(ref)]

test_list = list()
ref_list = list()
ref_info_list = list()
for(i in 1:CPUNUM){
  pick = slipt_vec == i
  ref_list[[i]] = ref[pick,]
  test_list[[i]] = test[pick,]
  ref_info_list[[i]] = ref_info[pick,]
}

for(i in 1:CPUNUM){

  write_out = list(ref = ref_list[[i]], test = test_list[[i]], ref_info = ref_info_list[[i]])
  save(write_out, file = paste0(out_file, "_part", i, ".Rdata"))
  cat( paste0("DONE preprocessing part ", i , "\n"))
}

cat("DONE preprocessing ALL \n")