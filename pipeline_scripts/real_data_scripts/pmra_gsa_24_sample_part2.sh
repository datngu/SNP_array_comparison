#!/bin/bash
#SBATCH --job-name=imputation_array
#SBATCH --output=imputation_array.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=512G


# IMPUTATION


impute_fun ()
{
  array_dir=$1
  ref_dir=$2
  m3vcf_dir=$3
  i=$4
  ##
  vcf=${array_dir}/chr${i}.vcf.gz
  ref_vcf=${ref_dir}/chr${i}.vcf.gz
  m3vcf=${m3vcf_dir}/chr${i}.m3vcf.gz
  map=/dragennfs/area7/datnguyen/bioinformatics_tools/shapeit4-4.1.3/maps_66/genetic_maps.b38/chr${i}.b38.gmap.gz
  ## rsID
  dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
  ## phasing
  shapeit4 --input $vcf \
  --map $map \
  --reference $ref_vcf \
  --region ${i} \
  --thread 1 \
  --seed 2022 \
  --output ${array_dir}/chr${i}_phased.vcf.gz  

  # imputing
  minimac4 --refHaps $m3vcf \
    --ChunkLengthMb 50 \
    --ChunkOverlapMb 5 \
    --haps ${array_dir}/chr${i}_phased.vcf.gz \
    --prefix ${array_dir}/chr${i}_imputed_tem \
    --ignoreDuplicates \
    --cpus 4 \
    --vcfBuffer 1100
  
  # annotate db151
  bcftools index -t ${array_dir}/chr${i}_imputed_tem*.vcf.gz
  bcftools annotate -a $dbsnp151 -c ID ${array_dir}/chr${i}_imputed_tem*.vcf.gz -o ${array_dir}/chr${i}_imputed.vcf.gz -Oz
  bcftools index -t ${array_dir}/chr${i}_imputed.vcf.gz
  rm ${array_dir}/chr${i}_imputed_tem*
}


### PROCESS IMPUTED DATA


process_imputed_file ()
{
  array_dir=$1
  i=$2
  fi=${array_dir}/chr${i}_imputed.vcf.gz

  out_dose=${array_dir}/chr${i}_imputed.dose.txt
  echo $'ID' $(bcftools view -h $fi | grep "^#CHROM" | cut -f10-) | sed -e 's/ /\t/g' > $out_dose
  bcftools  query $fi -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n' >> $out_dose
  gzip -f $out_dose
  echo "DONE DS chr${i}, array ${array_dir}"

  out_GT=${array_dir}/chr${i}_imputed.GT.txt
  echo $'ID' $(bcftools view -h $fi | grep "^#CHROM" | cut -f10-) | sed -e 's/ /\t/g' > $out_GT
  bcftools query $fi -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n' | sed 's/0\/0/0/g' | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' >> $out_GT
  gzip -f $out_GT
  echo "DONE GT chr${i}, array ${array_dir}"
}


### PLINK FILES


merge_vcf ()
{
  array_dir=$1
  #process_fam=/dragennfs/area8/datnguyen/PRS_arrays_project/KHV_10_folds_imputed/merged_all_autosome_1013_GATK.fam
  base_name=merged_all_autosome
  #dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
  ls ${array_dir}/*imputed.vcf.gz > ${array_dir}/file.txt
  mkdir ${array_dir}/plink
  bcftools concat -n -f ${array_dir}/file.txt -Oz -o ${array_dir}/plink/${base_name}_db151.vcf.gz
}

# # process fam files

# options(stringsAsFactors = FALSE)
# p = read.table("/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS/plink_WGS/process_vn1008.fam")
# s = read.table("shared_sample.txt")
# p = p[p$V1 %in% s$V1,]
# od = match(s$V1, p$V1)
# p2 = p[od,]
# write.table(p2, file = "shared_sample.fam", sep = "\t", col.names = F, quote = F, row.names = F)

process_array ()
{
  array_dir=$1
  process_fam=/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/gsa_real_data/shared_sample.fam
  base_name=merged_all_autosome
  #dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
  #ls ${array_dir}/*.vcf.gz > ${array_dir}/file.txt
  #mkdir ${array_dir}/plink
  #bcftools concat -n -f ${array_dir}/file.txt -Oz -o ${array_dir}/plink/${base_name}.vcf.gz
  #bcftools index -t ${array_dir}/plink/${base_name}.vcf.gz
  #bcftools annotate -a $dbsnp151 -c ID ${array_dir}/plink/${base_name}.vcf.gz | bgzip > ${array_dir}/plink/${base_name}_db151.vcf.gz

  plink --vcf ${array_dir}/plink/${base_name}_db151.vcf.gz \
      --make-bed  --const-fid --out ${array_dir}/plink/${base_name}_raw \
      --threads 2 \
      --memory 640000

    
    # copy processed fam with FID and Sex
    cat ${process_fam} > ${array_dir}/plink/${base_name}_raw.fam

    plink \
    --bfile ${array_dir}/plink/${base_name}_raw \
    --maf 0.0001 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --memory 640000 \
    --out ${array_dir}/plink/${base_name}.QC

    r_code='
  #!/usr/bin/env Rscript\n
  args = commandArgs(trailingOnly=TRUE)\n
  options(stringsAsFactors = FALSE)\n
  require(data.table)\n
  in_file = args[1]\n
  out_file = args[2]\n
  snp = fread(in_file, header = F)\n
  pick = snp$V1 == "."\n
  snp = snp[!pick,]\n
  pick = duplicated(snp$V1)\n
  snp_filter = snp[!pick,]\n
  fwrite(snp_filter, file = out_file, sep = "\t")\n
  cat("DONE merging files")\n
  '
  echo -e $r_code > ${array_dir}/filter_duplicated.R
  chmod 777 ./${array_dir}/filter_duplicated.R
  ./${array_dir}/filter_duplicated.R ${array_dir}/plink/${base_name}.QC.snplist ${array_dir}/plink/${base_name}.QC.snplist.nodup  
  

  plink \
      --bfile ${array_dir}/plink/${base_name}_raw \
      --threads 2 \
      --make-bed \
      --keep ${array_dir}/plink/${base_name}.QC.fam \
      --out ${array_dir}/plink/${base_name}.QC \
      --extract ${array_dir}/plink/${base_name}.QC.snplist.nodup \
      --memory 640000

  rm ${array_dir}/plink/${base_name}_raw*
  #rm ${array_dir}/plink/${base_name}.vcf.gz
  rm ${array_dir}/plink/${base_name}.QC.snplis*
  rm ${array_dir}/plink/${base_name}_db151.vcf.gz
}



# WGS
array_dir=wgs_24_sample

merge_vcf_wgs ()
{
  array_dir=$1
  #process_fam=/dragennfs/area8/datnguyen/PRS_arrays_project/KHV_10_folds_imputed/merged_all_autosome_1013_GATK.fam
  base_name=merged_all_autosome
  #dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
  ls ${array_dir}/*.vcf.gz > ${array_dir}/file.txt
  mkdir ${array_dir}/plink
  bcftools concat -n -f ${array_dir}/file.txt -Oz -o ${array_dir}/plink/${base_name}_db151.vcf.gz
}

merge_vcf_wgs $array_dir
process_array $array_dir

#######################################################
### ARRAY GSA

# real gsa
array_dir=gsa_real_array
ref_dir=reference_set
m3vcf_dir=reference_set/m3vcf

## testing
#i=22
#impute_fun $array_dir $ref_dir $m3vcf_dir $i
#process_imputed_file $array_dir $i

###
## Imputation
for i in {1..22}
do
  impute_fun $array_dir $ref_dir $m3vcf_dir $i &
done
wait

## extract GT and dose
for i in {1..22}
do
  process_imputed_file $array_dir $i &
done
wait

## Plinks
merge_vcf $array_dir
process_array $array_dir


# simulated gsa

array_dir=gsa_simulated_array
ref_dir=reference_set
m3vcf_dir=reference_set/m3vcf
#i=22

###
## Imputation
for i in {1..22}
do
  impute_fun $array_dir $ref_dir $m3vcf_dir $i &
done
wait

## extract GT and dose
for i in {1..22}
do
  process_imputed_file $array_dir $i &
done
wait

merge_vcf $array_dir
process_array $array_dir


#######################################################
### ARRAY PMRA

# real pmra
array_dir=pmra_real_array
ref_dir=reference_set
m3vcf_dir=reference_set/m3vcf

## testing
#i=22
#impute_fun $array_dir $ref_dir $m3vcf_dir $i
#process_imputed_file $array_dir $i

###
## Imputation
for i in {1..22}
do
  impute_fun $array_dir $ref_dir $m3vcf_dir $i &
done
wait

## extract GT and dose
for i in {1..22}
do
  process_imputed_file $array_dir $i &
done
wait

## Plinks
merge_vcf $array_dir
process_array $array_dir


# simulated pmra

array_dir=pmra_simulated_array
ref_dir=reference_set
m3vcf_dir=reference_set/m3vcf
#i=22

###
## Imputation
for i in {1..22}
do
  impute_fun $array_dir $ref_dir $m3vcf_dir $i &
done
wait

## extract GT and dose
for i in {1..22}
do
  process_imputed_file $array_dir $i &
done
wait

merge_vcf $array_dir
process_array $array_dir


echo "DONE ALLLL!!!!!!!!!!!!!"






