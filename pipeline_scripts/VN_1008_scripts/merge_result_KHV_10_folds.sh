#!/bin/bash
#SBATCH --job-name=merge_GT_KHV_all
#SBATCH --output=merge_GT_KHV_all.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=128000


cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008


array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

# for array_dir in $array_list
# do
# 	rm $array_dir/*.tbi
# 	rm $array_dir/*.gz
# done

#array_dir=Axiom_GW_ASI

sample_id=KHV_ref_sample.txt

fi=/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/VN_1008.chr7.all.vcf.gz
bcftools view -h $fi | grep "^#CHROM" | cut -f10- | sed 's/\t/\n/g' > $sample_id



#     out_dir=/dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR
#     dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
#     file=chr${i}_2504_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz
#     bcftools annotate -a $dbsnp151 -c ID $file | bgzip > ${out_dir}/chr${i}.2504_IGSR_annotated_dbSNP151.vcf.gz

merge_fun ()
{
	array_dir=$1
	i=$2
	sammple=$3
	out_tem=$array_dir/chr${i}_merged_10_batchs_tem.vcf.gz
	out=$array_dir/chr${i}_merged_10_batchs.vcf.gz
	dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
	for n in {1..10}
	do
	  bcftools index -f --tbi $array_dir/pseudo_array_imputed/chr${i}_batch_${n}.dose.vcf.gz &
	done
	wait
	echo "merging chr $i, $array_dir"
	bcftools merge $array_dir/pseudo_array_imputed/chr${i}_batch_1.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_2.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_3.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_4.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_5.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_6.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_7.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_8.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_9.dose.vcf.gz $array_dir/pseudo_array_imputed/chr${i}_batch_10.dose.vcf.gz | bcftools view -S $sammple | bgzip > $out_tem
	bcftools index -f --tbi $out_tem
	bcftools annotate -a $dbsnp151 -c ID $out_tem | bgzip > $out
	bcftools index -f --tbi $out
	rm $out_tem
	rm ${out_tem}.tbi
	echo "DONE merging chr $i, $array_dir"
}

for i in {22..1}
do
	for array_dir in $array_list
	do
		merge_fun $array_dir $i $sample_id &
	done
	wait
done


wait

sleep 5m


process_population ()
{
	i=$1
	array_dir=$2
	pop=KHV
	#in_vcf=${array_dir}/chr${i}_merged_10_batchs.vcf.gz
	#sample=/dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/VCF_sample_information/supper_population_${pop}.txt
	#bcftools view $in_vcf -S $sample | bgzip > ${array_dir}/chr${i}_${pop}_imputed.vcf.gz
	echo "getting DS chr${i}, array ${array_dir}"
	
	fi=${array_dir}/chr${i}_merged_10_batchs.vcf.gz
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
  for i in {22..1}
  do
    process_population $i $array_dir &
  done
  wait
}




for array_dir in $array_list
do
	get_DS_GT_array $array_dir 
done

# for array_dir in $array_list
# do
# 	rm $array_dir/*.gz 
# 	rm $array_dir/*tbi
# done