#!/bin/bash
#SBATCH --job-name=cor_affy6_IGSR
#SBATCH --output=cor_igsr-affy6.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=160
#SBATCH --mem=700000

# array_dir=global-screening-array-v.3
# i=22
# bash /dragennfs/area7/datnguyen/bioinformatics_tools/compute_cor/pipeline_compute_correlation_KHV.sh /dragennfs/area7/datnguyen/PRS_arrays_project/1000_KHV_process_GT/processed_GT_chr${i}_1013_KHV.txt.gz /dragennfs/area7/datnguyen/PRS_arrays_project/KHV_imputed/${array_dir}/chr${i}.dose.txt.gz /dragennfs/area7/datnguyen/PRS_arrays_project/KHV_imputed/${array_dir}/chr${i} &

cd /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed

array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

array_list='GenomeWideSNP_6.0'



compute_cor_script=/dragennfs/area7/datnguyen/bioinformatics_tools/compute_cor_IGSR/pipeline_compute_correlation_IGSR.sh

for array_dir in $array_list
do
	for i in {22..1}
	do
		pop=EUR
		bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt &
	done
	wait
	echo "done $array_dir"
done

wait

for array_dir in $array_list
do
	for i in {22..1}
	do
		pop=SAS
		bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt &
	done
	wait
	echo "done $array_dir"
done

wait


for array_dir in $array_list
do
	for i in {22..1}
	do
		pop=EAS
		bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt &
	done
	wait
	echo "done $array_dir"
done

wait


for array_dir in $array_list
do
	for i in {22..1}
	do
		pop=AMR
		bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt &
	done
	wait
	echo "done $array_dir"
done

wait


## save RAM AFR
for array_dir in $array_list
do
	for i in {22..10}
	do
		pop=AFR
		bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt &
	done
	wait
	echo "done $array_dir"
done

wait

for array_dir in $array_list
do
	for i in {9..1}
	do
		pop=AFR
		bash $compute_cor_script /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/${pop}_process_GT/chr${i}_GT.txt.gz ${array_dir}/chr${i}_${pop}_imputed.dose.txt.gz ${array_dir}/chr${i}_${pop}_correlation.txt &
	done
	wait
	echo "done $array_dir"
done

wait

