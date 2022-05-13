#!/bin/bash
#SBATCH --job-name=cor_VGC
#SBATCH --output=cor_VGC.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=120
#SBATCH --mem=700000

# array_dir=global-screening-array-v.3
# i=22
# bash /dragennfs/area7/datnguyen/bioinformatics_tools/compute_cor/pipeline_compute_correlation_KHV.sh /dragennfs/area7/datnguyen/PRS_arrays_project/1000_KHV_process_GT/processed_GT_chr${i}_1013_KHV.txt.gz /dragennfs/area7/datnguyen/PRS_arrays_project/KHV_imputed/${array_dir}/chr${i}.dose.txt.gz /dragennfs/area7/datnguyen/PRS_arrays_project/KHV_imputed/${array_dir}/chr${i} &


#array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'
array_list='vingenchip'
for array_dir in $array_list
do
	for i in {22..1}
	do
		bash /dragennfs/area7/datnguyen/bioinformatics_tools/compute_cor/pipeline_compute_correlation_KHV.sh /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/chr${i}_GT.txt.gz /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008/${array_dir}/chr${i}_KHV_imputed.dose.txt.gz /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008/${array_dir}/chr${i}_KHV_correlation.txt.gz &
	done
	wait
	echo "done $array_dir"
done


# res_dir=/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/cor_res_vn1008
# mkdir $res_dir
# cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008
# for array_dir in $array_list
# do
#   mkdir ${res_dir}/${array_dir}
#   for i in {1..22}
#   do
#     cp ${array_dir}/chr${i}_KHV_correlation.txt.gz ${res_dir}/${array_dir}/chr${i}_VNP_correlation.txt.gz
#   done
# done

