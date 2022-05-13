#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
syntax='\nUsage:\t./this_sript.R Rdata'


if(length(args) == 0 ){
  cat("\nNo argument, Program stop! \n")
  cat(syntax)
  quit()
}

options(stringsAsFactors = FALSE)
require(data.table)


Rdata = args[1]

load(Rdata)

get_correlation <- function(ref, test, ref_info)
{
  ref_info$cor = 0
  for( i in 1:nrow(ref_info)){
  ref_info$cor[i] = cor(ref[i,], test[i,], method = "pearson")
  }
  ref_info$r_2 = ref_info$cor^2
  return(ref_info)
}

cat("\n computing correaltion\n")
cat(Rdata)

res = get_correlation(write_out$ref, write_out$test, write_out$ref_info)
outname = gsub(".Rdata", ".cor", Rdata)
fwrite(res, outname, sep = "\t")

cat("DONE\n")