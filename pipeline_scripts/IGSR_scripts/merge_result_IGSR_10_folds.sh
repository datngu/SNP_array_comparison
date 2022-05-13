#!/bin/bash
#SBATCH --job-name=merge_IGSR
#SBATCH --output=merge_IGSR.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=128000


cd /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed


array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

#array_dir=Axiom_GW_ASI


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
	for array_dir in $array_list
	do
		merge_fun $array_dir $i &
	done
	wait
done


