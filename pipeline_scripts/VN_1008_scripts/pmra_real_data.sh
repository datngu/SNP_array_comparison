#!/bin/bash
#SBATCH --job-name=processing_data
#SBATCH --output=processing_data.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=512G



cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/pmra_real_data

# extract samples
all_vcf=/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/VN_1008.chr10.all.vcf.gz

bcftools view -h $all_vcf | grep "^#CHROM" | cut -f10- | sed 's/\t/\n/g' > all_sample.txt


## annotation db151 pmra real data
pmra_raw_vcf=renamed_QC_95_sample_Hg38_PMRA.vcf.gz
pmra_vcf=95_sample_Hg38_PMRA_db151.vcf.gz

#dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
#bcftools annotate -a $dbsnp151 -c ID $pmra_raw_vcf -o $pmra_vcf -Oz
#bcftools index -t $pmra_vcf



bcftools view -h $pmra_vcf | grep "^#CHROM" | cut -f10- | sed 's/\t/\n/g' > pmra_sample.txt

# get shared samples

grep -Fxf all_sample.txt pmra_sample.txt > shared_sample.txt



### extract sample vcf

mkdir reference_set
mkdir 95_shared_sample_wgs
mkdir 95_shared_sample_pmra


#####################################################################################
# process reference

sample=shared_sample.txt
for i in {1..22}
do
    bcftools view /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/VN_1008.chr${i}.all.vcf.gz -S ^${sample} -o reference_set/chr${i}.vcf.gz -Oz &
    bcftools view /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/VN_1008.chr${i}.all.vcf.gz -S ${sample} -o 95_shared_sample_wgs/chr${i}.vcf.gz -Oz &
done

wait


for i in {1..22}
do
    bcftools index -t reference_set/chr${i}.vcf.gz &
done
wait


mkdir reference_set/m3vcf
# index m3vcf
for i in {1..22}
do
  Minimac3 --refHaps reference_set/chr${i}.vcf.gz --processReference --prefix reference_set/m3vcf/chr${i} --log reference_set/m3vcf/chr${i}_log.txt &
done

wait

# processing reference for imputation accuracy estimation
# note: it is 95 sample WGS,


process_ref ()
{
  vcf=$1
  out_GT=$2
  echo $'ID\tAF\tMAF\tAC' $(bcftools view -h $vcf | grep "^#CHROM" | cut -f10-) | sed -e 's/ /\t/g' > $out_GT
  bcftools query $vcf -f '%CHROM:%POS:%REF:%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/AC[\t%GT]\n' | sed 's/0\/0/0/g' | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | awk '($4 > 1) {print}' >> $out_GT
}

mkdir 95_shared_sample_wgs/ref_GT
for i in {1..22}
do
  process_ref 95_shared_sample_wgs/chr${i}.vcf.gz 95_shared_sample_wgs/ref_GT/chr${i}.txt &
done
wait

for i in {1..22}
do
  bgzip 95_shared_sample_wgs/ref_GT/chr${i}.txt &
done
wait





#####################################################################################
## process real pmra

sample=shared_sample.txt
for i in {1..22}
do
    bcftools view $pmra_vcf -r $i -S ${sample} -o 95_shared_sample_pmra/chr${i}.vcf.gz -Oz &
done

wait

for i in {1..22}
do
    bcftools index -t 95_shared_sample_pmra/chr${i}.vcf.gz &
done

wait



#####################################################################################
# processing simulated pmra


mkdir 95_shared_sample_simulated

put /media/datn/data/1st_DC_PRS_array_project/array_annotation/annotation_all_hg38/Axiom_PMRA/* /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/pmra_real_data/95_shared_sample_simulated

## simulate non-phasing data



#sample=shared_sample.txt
for i in {1..22}
do
  vcf_wgs=95_shared_sample_wgs/chr${i}.vcf.gz
  position=95_shared_sample_simulated/chr${i}.txt
  vcftools --gzvcf  $vcf_wgs --positions $position --recode --stdout | sed 's/0|0/0\/0/g' | sed 's/0|1/0\/1/g' | sed 's/1|0/1\/0/g' | sed 's/1|1/1\/1/g' | bgzip > 95_shared_sample_simulated/chr${i}.vcf.gz
  bcftools index -t 95_shared_sample_simulated/chr${i}.vcf.gz
done

#####################################################################################
#####################################################################################
#####################################################################################






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
  bcftools query $fi -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n' >> $out_dose
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

# process fam files

# options(stringsAsFactors = FALSE)
# p = read.table("/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS/plinkWGS/process_vn1008.fam")
# s = read.table("shared_sample.txt")
# p = p[p$V1 %in% s$V1,]
# od = match(s$V1, p$V1)
# p2 = p[od,]
# write.table(p2, file = "shared_sample.fam", sep = "\t", col.names = F, quote = F, row.names = F)

process_array ()
{
  array_dir=$1
  process_fam=/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/pmra_real_data/shared_sample.fam
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



# plink wgs
array_dir=95_shared_sample_wgs

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


### real pmra

array_dir=95_shared_sample_pmra
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

array_dir=95_shared_sample_simulated
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



####################################3
# compute corelations and PRS

#!/bin/bash
#SBATCH --job-name=PRS
#SBATCH --output=PRS.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=512G




## PRS
cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/pmra_real_data

array_list="95_shared_sample_simulated 95_shared_sample_pmra"

mkdir PRS_result

get_PRS ()
{
  array_dir=$1
  trait=$2
  pop=VNP
  base_name=merged_all_autosome
  sumstat=/dragennfs/area7/datnguyen/bioinformatics_tools/sumstat/GIANT_${trait}.QCed.gz
  bedfiles=${array_dir}/plink/${base_name}.QC
  #bedfiles=${array_dir}/plink_${array_dir}/${base_name}.QC
  outfn=PRS_result/${array_dir}_${pop}_${trait}
  #PRSice_linux
  /dragennfs/area7/datnguyen/bioinformatics_tools/PRSice_linux/PRSice_linux \
    --base ${sumstat} \
    --target ${bedfiles} \
    --out ${outfn} \
    --binary-target F \
    --bar-levels 0.00000005,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.2,0.3,0.5,1 --fastscore \
    --a1 A1 \
    --a2 A2 \
    --beta  \
    --bp BP \
    --chr CHR \
    --pvalue P \
    --snp SNP \
    --stat BETA \
    --clump-kb 250kb \
    --clump-p 1 \
    --clump-r2 0.1 \
    --ultra \
    --no-regress \
    --score sum \
    --thread 1
}



## imputation correlation


array_list="95_shared_sample_simulated 95_shared_sample_pmra 95_shared_sample_wgs"



for array_dir in $array_list
do
  for i in {22..1}
  do
    bash /dragennfs/area7/datnguyen/bioinformatics_tools/compute_cor/pipeline_compute_correlation_KHV.sh /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/pmra_real_data/95_shared_sample_wgs/ref_GT/chr${i}.txt.gz /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/pmra_real_data/${array_dir}/chr${i}_imputed.dose.txt.gz /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/pmra_real_data/${array_dir}/chr${i}_imputed.correlation.txt.gz &
  done
  wait
  echo "done $array_dir"
done

wait


# PRS


for array_dir in $array_list
do
  get_PRS $array_dir HEIGHT
  get_PRS $array_dir BMI
  get_PRS $array_dir Type_2_diabetes
done




### copy results
get -r PRS_result /media/datn/data2gb/GitHub/SNP_array_comparsion/pmra_real_data

# local
array_list="95_shared_sample_simulated 95_shared_sample_pmra 95_shared_sample_wgs"
for array_dir in $array_list
do
  mkdir /media/datn/data2gb/GitHub/SNP_array_comparsion/pmra_real_data/$array_dir
done

# remote

get 95_shared_sample_simulated/*cor*gz /media/datn/data2gb/GitHub/SNP_array_comparsion/pmra_real_data/95_shared_sample_simulated

get 95_shared_sample_pmra/*cor*gz /media/datn/data2gb/GitHub/SNP_array_comparsion/pmra_real_data/95_shared_sample_pmra

get 95_shared_sample_wgs/*cor*gz /media/datn/data2gb/GitHub/SNP_array_comparsion/pmra_real_data/95_shared_sample_wgs



 