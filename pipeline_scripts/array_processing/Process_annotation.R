#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
  cat("Syntax: Rscript axiom/illumina path_annotation out_bed\n")
  quit()
}
array_type = args[1]
anno_path = args[2]
out_name = args[3]

#setwd("/media/datn/data/1st_DC_PRS_array_project/array_annotation")
if(!require(data.table)) install.packages("data.table")
library(data.table)

#################
#anno_path = "/media/datn/data/1st_DC_PRS_array_project/array_annotation/Axiom_GW_ASI/Axiom_GW_ASI_SNP.na35.annot.csv"
#out_name = "Axiom_GW_ASI_GRCh37.bed"
get_bed_format_axiom <- function(anno_path, out_name, gene_build = "hg38"){
  cat("reading your input annotation...\n")
  cat(anno_path)
  cat("\n")
  anno = fread(anno_path, skip = "#")
  pick_col = c("Chromosome", "Physical Position", "Physical Position", "Ref Allele", "Alt Allele")
  anno = anno[, ..pick_col]
  fwrite(anno, file = out_name, sep = "\t", col.names = F, row.names = F, quote = F)
  cat("Your output bed file is: \n")
  cat(out_name)
  cat("\n")
}

get_bed_format_illumina <- function(anno_path, out_name, gene_build = "hg38"){
  cat("reading your input annotation...\n")
  cat(anno_path)
  cat("\n")
  anno = fread(anno_path, skip = 7)
  pick_col = c("Chr", "MapInfo", "MapInfo")
  anno = anno[, ..pick_col]
  fwrite(anno, file = out_name, sep = "\t", col.names = F, row.names = F, quote = F)
  cat("Your output bed file is: \n")
  cat(out_name)
  cat("\n")
}

if(array_type == "axiom"){
  get_bed_format_axiom(anno_path, out_name)
}else if(array_type == "illumina"){
  get_bed_format_illumina(anno_path, out_name)
}else{
  cat("Syntax: Rscript axiom/illumina path_annotation out_bed\n")
  quit()
}










# cd /media/datn/data/1st_DC_PRS_array_project/array_annotation
# GRCh37="Axiom_GW_ASI_GRCh37.bed"
# GRCh38="Axiom_GW_ASI_GRCh38.bed"
# chain_file="/media/datn/data/cross_map_lift_over/GRCh37_to_GRCh38.chain.gz"

# CrossMap.py region $chain_file $GRCh37 $GRCh38
