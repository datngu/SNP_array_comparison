#!/bin/bash
#SBATCH --job-name=IGSR_GDA
#SBATCH --output=IGSR_GDA.out
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=60
#SBATCH --mem=512000




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

for i in {22..1}
do
	process_10_folds $array_dir $i
done

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

for i in {22..1}
do
	merge_fun $array_dir $i &
done
wait



############# get DS, GT

get_DS ()
{
	echo 
	i=$1
	array_dir=$2
	echo "getting DS chr${i}, array ${array_dir}"
	#ref_id_file=/dragennfs/area7/datnguyen/PRS_arrays_project/1000_KHV_process_GT/processed_GT_chr${i}_1013_KHV_ID_only.txt
	fi=${array_dir}/chr${i}_merged_10_batchs.vcf.gz
	out=${array_dir}/chr${i}_merged_10_batchs.dose.txt
	echo $'ID' $(bcftools view -h $fi | grep "^#CHROM" | cut -f10-) | sed -e 's/ /\t/g' > $out
	bcftools query $fi -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n' >> $out
	gzip -f $out
	echo "DONE DS chr${i}, array ${array_dir}"
}


get_GT ()
{
	i=$1
	array_dir=$2
	echo "getting GT chr${i}, array ${array_dir}"
	#ref_id_file=/dragennfs/area7/datnguyen/PRS_arrays_project/1000_KHV_process_GT/processed_GT_chr${i}_1013_KHV_ID_only.txt
	fi=${array_dir}/chr${i}_merged_10_batchs.vcf.gz
	out=${array_dir}/chr${i}_merged_10_batchs.GT_encoded.txt
	echo $'ID' $(bcftools view -h $fi | grep "^#CHROM" | cut -f10-) | sed -e 's/ /\t/g' > $out
	bcftools query $fi -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n' | sed 's/0\/0/0/g' | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' >> $out
	gzip -f $out
	echo "DONE GT chr${i}, array ${array_dir}"
}


get_DS_GT_array ()
{
  array_dir=$1
  for i in {1..22}
  do
    get_GT $i $array_dir &
    get_DS $i $array_dir &
  done
  wait
}


get_DS_GT_array $array_dir

wait

/dragennfs/area7/datnguyen/bioinformatics_tools/compute_cor_IGSR/filter_pop_IGSR.R /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed/${array_dir}

echo "done $array_dir"