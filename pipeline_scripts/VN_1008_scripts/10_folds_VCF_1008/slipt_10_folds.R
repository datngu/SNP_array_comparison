setwd("/media/datn/data/1st_DC_PRS_array_project/VN_1008_scripts/10_folds_VCF_1008")

# bcftools view -h /home/datn/NAS_private/Ref_impute/VN1008/VN_1008.chr8.all.vcf.gz | grep "^#CHROM" | cut -f10- | sed 's/\t/\n/g' > sampleID_by_line.txt

data = read.table("sampleID_by_line.txt")

fold_size = 10
data$fold_id = c(1:nrow(data))
set.seed(2021)
id  = sample(1:nrow(data))
data = data[id,]
step = nrow(data) %/% fold_size + 1 
id = rep(1:fold_size, step)
data$fold_id = id[1:nrow(data)]
for(i in 1:fold_size){
  out = paste0("batch_", i, ".txt")
  tem = data[data$fold_id == i,]
  tem = as.data.frame(tem[,1])
  write.table(tem, file = out, sep = "\t", quote = F, row.names = F, col.names = F)
}
