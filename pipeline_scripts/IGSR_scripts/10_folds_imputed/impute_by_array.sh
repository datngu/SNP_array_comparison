#!/bin/bash
#SBATCH --job-name=khv_imp-10f
#SBATCH --output=khv_imp-10f.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=500000

array_dir=$1

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


cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008

# nano impute_by_array.sh

# put -r /media/datn/data/1st_DC_PRS_array_project/array_annotation/array-used_by_chr/* /.

impute_fun ()
{
  i=$1
  batch=$2
  array_dir=$3
  
  echo "doning chr${i}, array ${array_dir}, batch ${batch}"
  vcf_test=/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/10_fold_VCF_1008/batch_${batch}/chr${i}_test_set.vcf.gz
  vcf_ref=/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/10_fold_VCF_1008/batch_${batch}/chr${i}_reference.vcf.gz
  m3vcf=/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/10_fold_VCF_1008/m3vcf_batch_${batch}/chr${i}_reference.m3vcf.gz
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
  map=/dragennfs/area7/datnguyen/bioinformatics_tools/shapeit4-4.1.3/maps_66/genetic_maps.b38/chr${i}.b38.gmap.gz
  shapeit4 --input ${array_dir}/pseudo_array/chr${i}_batch_${batch}.vcf.gz \
    --map $map \
    --reference $vcf_ref \
    --region ${i} \
    --thread 1 \
    --output ${array_dir}/pseudo_array/chr${i}_batch_${batch}_phased.vcf.gz  

  # imputing
  minimac4 --refHaps $m3vcf \
    --ChunkLengthMb 50 \
    --ChunkOverlapMb 5 \
    --haps ${array_dir}/pseudo_array/chr${i}_batch_${batch}_phased.vcf.gz \
    --prefix ${array_dir}/pseudo_array_imputed/chr${i}_batch_${batch} \
    --ignoreDuplicates \
    --cpus 4 \
    --vcfBuffer 1100

  echo "DONE chr${i}; array ${array_dir}, batch ${batch}"
}



process_by_chr_array ()
{
    i=$1
    array_dir=$2
    for batch in {1..10}
    do
        impute_fun $i $batch $array_dir &
    done
    wait
}


#array_list1='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k'
#array_list2='Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2'


for i in {1..22}
do
  process_by_chr_array $i $array_dir
done


echo "DONE $array_dir"

#36866
