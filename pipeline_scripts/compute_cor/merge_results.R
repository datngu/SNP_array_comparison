#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
syntax='\nUsage:\t./this_sript.R out_file'


if(length(args) == 0 ){
  cat("\nNo argument, Program stop! \n")
  cat(syntax)
  quit()
}

options(stringsAsFactors = FALSE)
require(data.table)

out_file = args[1]


l = list()
for(i in 1:10){
  file_in = paste0(out_file, "_part", i , ".cor")
  tem = fread(file_in)
  l[[i]] = tem
}
out = do.call(rbind, l)

fwrite(out, file = out_file, sep = "\t")

cat("DONE merging files \n")