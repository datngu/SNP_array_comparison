setwd("/media/datn/data/1st_DC_PRS_array_project/array_annotation/annotation_all_hg38")
require(data.table)

array_info = fread("array_size.txt")

all_bed = list()

chr_list = as.character(c(1:22))

for ( ar in array_info$array){
  tem = fread(paste0(ar, "_hg38.bed"))
  tem$V4 = ar
  tem$V2 = as.integer(tem$V2)
  tem$V3 = tem$V2
  tem$id = paste(tem$V1, tem$V2, sep = ":")
  pick = tem$V1 %in% chr_list & tem$V2 > 0
  tem = tem[pick,]
  all_bed[[ar]] = tem
  # all_pos = c(all_pos, all_bed[[ar]]$id)
}


df = do.call(rbind, all_bed)
df = df[,-c("V4")]
names(df) = c("chr", "start", "end", "id")
d = duplicated(df$id)
df = df[!d,]
df = as.data.frame(df)

for( ar in array_info$array){
  df[,ar] = 0
  df[,ar][which(df$id %in% all_bed[[ar]]$id)] = 1
}

class(df$start)

## sorting and indexing
for(i in chr_list){
  tem = df[df$chr == i,]
  od = order(tem$start)
  tem = tem[od,]
  if(i == "1"){
    res = tem
  }else{
    res = rbind(res, tem)
  } 
}

res = res[,-3]

fwrite(res, file = "db_array.txt", sep = "\t", row.names = F, col.names = F)
system("bgzip db_array.txt")
system("tabix -b 2 -e 2 db_array.txt.gz")
