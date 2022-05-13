#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
  cat("Syntax: Rscript path_to_bed_file")
  quit()
}

path_bed = args[1]

base_name = basename(path_bed)
out_dir = gsub("_hg38.bed", "", base_name)
dir.create(out_dir)
#setwd("/media/datn/data/1st_DC_PRS_array_project/array_annotation")
#if(!require(data.table)) install.packages("data.table")
#library(data.table)
bed = read.table(path_bed)
chr = c(1:22, "X", "Y")
bed = bed[,c(1,2)]
for(i in chr){
    tem = bed[bed[,1] == i,]
    dup = duplicated(tem[,2])
    #table(duqqp)
    tem = tem[!dup,]
    out = paste0(out_dir,"/chr", i, ".txt")
    write.table(tem, file = out, sep = "\t", quote = F, col.names = F, row.names = F)
}