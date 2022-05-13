#!/bin/bash
#SBATCH --job-name=IGSR_all
#SBATCH --output=IGSR_all.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=500000

export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/Minimac3Executable/bin:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/Minimac4/minimac4-1.0.2-Linux/bin:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/impute5_v1.1.4/impute5_v1.1.4:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/plink1.9:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/LmTag-0.0.0/bin:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/vcftools/src/cpp:$PATH
export PERL5LIB=/dragennfs/area7/datnguyen/bioinformatics_tools/vcftools/src/perl/
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/shapeit4-4.1.3/bin:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/bcftools-1.10.2:$PATH
export BCFTOOLS_PLUGINS=/dragennfs/area7/datnguyen/bioinformatics_tools/bcftools-1.10.2/plugins
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/TagIt/src:$PATH


cd /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed


impute_fun ()
{
  i=$1
  batch=$2
  array_dir=$3
  
  echo "doning chr${i}, array ${array_dir}, batch ${batch}"
  vcf_test=/dragennfs/area8/datnguyen/PRS_arrays_project/10_folds_VCF/batch_${batch}/chr${i}_test_set.vcf.gz
  vcf_ref=/dragennfs/area8/datnguyen/PRS_arrays_project/10_folds_VCF/batch_${batch}/chr${i}_reference.vcf.gz
  m3vcf=/dragennfs/area8/datnguyen/PRS_arrays_project/10_folds_VCF/m3vcf_batch_${batch}/chr${i}_reference.m3vcf.gz
  #ref_ID=/media/datn/data/1st_DC_PRS_array_project/1000_KHV_process_GT/processed_GT_chr${i}_1013_KHV_ID_only.txt

  position=${array_dir}/chr${i}.txt
  #ref_id_file=/media/datn/data/1st_DC_PRS_array_project/1000_KHV_process_GT/processed_GT_chr${i}_1013_KHV_ID_only.txt
  #/media/datn/data/1st_DC_PRS_array_project/Imputation_Arrays/phasing_imputing_array.sh $i $array_dir $vcf_test $vcf_ref $m3vcf $ref_id_file
  # extract pseudo array
  mkdir -p ${array_dir}/pseudo_array
  mkdir -p ${array_dir}/pseudo_array_imputed
  vcftools --gzvcf  $vcf_test --positions $position --recode --stdout | sed 's/0|0/0\/0/g' | sed 's/0|1/0\/1/g' | sed 's/1|0/1\/0/g' | sed 's/1|1/1\/1/g' | bgzip > ${array_dir}/pseudo_array/chr${i}_batch_${batch}.vcf.gz

  bcftools index -t ${array_dir}/pseudo_array/chr${i}_batch_${batch}.vcf.gz
  # phasing
  map=/dragennfs/area7/datnguyen/bioinformatics_tools/shapeit4-4.1.3/maps/genetic_maps.b38/chr${i}.b38.gmap.gz
  shapeit4 --input ${array_dir}/pseudo_array/chr${i}_batch_${batch}.vcf.gz \
    --map $map \
    --reference $vcf_ref \
    --region ${i} \
    --thread 4 \
    --output ${array_dir}/pseudo_array/chr${i}_batch_${batch}_phased.vcf.gz  

  # imputing
  minimac4 --refHaps $m3vcf \
    --ChunkLengthMb 200 \
    --ChunkOverlapMb 20 \
    --haps ${array_dir}/pseudo_array/chr${i}_batch_${batch}_phased.vcf.gz \
    --prefix ${array_dir}/pseudo_array_imputed/chr${i}_batch_${batch} \
    --ignoreDuplicates \
    --cpus 4 \
    --vcfBuffer 1100

  echo "DONE chr${i}; array ${array_dir}, batch ${batch}"
}


process_chr_all_array ()
{
    i=$1
    batch=$2
    for array_dir in Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_UKB_WCSG GenomeWideSNP_6.0 Axiom_PMRA Axiom_JAPONICA Axiom_PMDA infinium-omni2.5.v1.5 infinium-omni5-v1.2
    do
        impute_fun $i $batch $array_dir &
    done
    wait
}


# process_chr_all_array ()
# {
#     i=$1
#     batch=$2
#     for array_dir in chinese-genotyping-array-v1.0 infinium-core-v1.2 infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 multi-ethnic-eur-eas-sas-v1.0 oncoarray-500k cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-exome-v1.0 japanese-screening-array-v1.0 multi-ethnic-global-v1.0
#     do
#         impute_fun $i $batch $array_dir &
#     done
#     wait
# }

batch=1
for i in {1..22}  
do
  process_chr_all_array $i $batch
done

wait

batch=2
for i in {1..22}  
do
  process_chr_all_array $i $batch
done

wait

batch=3
for i in {1..22}  
do
  process_chr_all_array $i $batch
done

wait

batch=4
for i in {1..22}  
do
  process_chr_all_array $i $batch
done

wait

batch=5
for i in {1..22}  
do
  process_chr_all_array $i $batch
done

wait


batch=6
for i in {1..22}  
do
  process_chr_all_array $i $batch
done

wait

batch=7
for i in {1..22}  
do
  process_chr_all_array $i $batch
done

wait


batch=8
for i in {1..22}  
do
  process_chr_all_array $i $batch
done

wait

batch=9
for i in {1..22}  
do
  process_chr_all_array $i $batch
done

wait


batch=10
for i in {1..22}  
do
  process_chr_all_array $i $batch
done

wait
