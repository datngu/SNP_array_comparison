#!/bin/bash
#SBATCH --job-name=IGSR_GDA
#SBATCH --output=IGSR_GDA.out
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=128000




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



#######################
cd /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed

array_dir=infinium-global-diversity-array-v1.0


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

process_10_folds ()
{
    array_dir=$1
    i=$2
    for batch in {1..10}
    do
    	impute_fun $i $batch $array_dir &
    done
    wait
}

process_10_folds $array_dir 9


wait
# ##### merging


merge_fun ()
{
	array_dir=$1
	i=$2
	out=$array_dir/chr${i}_merged_10_batchs.vcf.gz

	for n in {1..10}
	do
	  bcftools index -f --tbi $array_dir/pseudo_array_imputed/chr${i}_batch_${n}.dose.vcf.gz &
	done
	wait
	echo "merging chr $i, $array_dir"
	bcftools merge $array_dir/pseudo_array_imputed/chr${i}_batch_1.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_2.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_3.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_4.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_5.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_6.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_7.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_8.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_9.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_10.dose.vcf.gz | bgzip > $out
	bcftools index -f --tbi $out
	echo "DONE merging chr $i, $array_dir"
}

merge_fun $array_dir 9

wait




############# get DS, GT

process_population ()
{
  i=$1
  array_dir=$2
  pop=$3
  in_vcf=${array_dir}/chr${i}_merged_10_batchs.vcf.gz
  sample=/dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/VCF_sample_information/supper_population_${pop}.txt
  bcftools view $in_vcf -S $sample | bgzip > ${array_dir}/chr${i}_${pop}_imputed.vcf.gz
  echo "getting DS chr${i}, array ${array_dir} pop $pop"
  
  fi=${array_dir}/chr${i}_${pop}_imputed.vcf.gz
  out_dose=${array_dir}/chr${i}_${pop}_imputed.dose.txt
  echo $'ID' $(bcftools view -h $fi | grep "^#CHROM" | cut -f10-) | sed -e 's/ /\t/g' > $out_dose
  bcftools query $fi -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n' >> $out_dose
  gzip -f $out_dose
  echo "DONE DS chr${i}, array ${array_dir}"

  out_GT=${array_dir}/chr${i}_${pop}_imputed.GT.txt
  echo $'ID' $(bcftools view -h $fi | grep "^#CHROM" | cut -f10-) | sed -e 's/ /\t/g' > $out_GT
  bcftools query $fi -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n' | sed 's/0\/0/0/g' | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' >> $out_GT
  gzip -f $out_GT
  echo "DONE GT chr${i}, array ${array_dir}"
}

get_DS_GT_array ()
{
  array_dir=$1
  for i in 9
  do
    process_population $i $array_dir EUR &
    process_population $i $array_dir EAS &
    process_population $i $array_dir SAS &
    process_population $i $array_dir AMR &
    process_population $i $array_dir AFR &
  done
  wait
}

get_DS_GT_array $array_dir 




array_list='infinium-global-diversity-array-v1.0'



compute_cor_script=/dragennfs/area7/datnguyen/bioinformatics_tools/compute_cor_IGSR/pipeline_compute_correlation_IGSR.sh

for array_dir in $array_list
do
  for i in 9
  do
    pop=EUR
    bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt.gz &
  done
  wait
  echo "done $array_dir"
done

wait

for array_dir in $array_list
do
  for i in 9
  do
    pop=SAS
    bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt.gz &
  done
  wait
  echo "done $array_dir"
done

wait


for array_dir in $array_list
do
  for i in 9
  do
    pop=EAS
    bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt.gz &
  done
  wait
  echo "done $array_dir"
done

wait


for array_dir in $array_list
do
  for i in 9
  do
    pop=AMR
    bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt.gz &
  done
  wait
  echo "done $array_dir"
done

wait


## save RAM AFR
for array_dir in $array_list
do
  for i in 9
  do
    pop=AFR
    bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt.gz &
  done
  wait
  echo "done $array_dir"
done

wait
