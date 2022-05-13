#!/bin/bash
#SBATCH --job-name=DS-GT_IGSR
#SBATCH --output=DS-GT_IGSR.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=120
#SBATCH --mem=256000


cd /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed


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
  for i in {1..22}
  do
    process_population $i $array_dir EUR &
    process_population $i $array_dir EAS &
    process_population $i $array_dir SAS &
    process_population $i $array_dir AMR &
    process_population $i $array_dir AFR &
  done
  wait
}

array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'



for array_dir in $array_list
do
	get_DS_GT_array $array_dir 
done

