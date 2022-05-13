#!/bin/bash
#SBATCH --job-name=filter_pop
#SBATCH --output=filter_pop.out
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=60
#SBATCH --mem=512000



cd /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed

for array_dir in Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_UKB_WCSG GenomeWideSNP_6.0 Axiom_PMRA Axiom_JAPONICA Axiom_PMDA chinese-genotyping-array-v1.0 infinium-core-v1.2 infinium-omni2.5.v1.5 infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 multi-ethnic-eur-eas-sas-v1.0 oncoarray-500k cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-exome-v1.0 infinium-omni5-v1.2 japanese-screening-array-v1.0 multi-ethnic-global-v1.0
do
	/dragennfs/area7/datnguyen/bioinformatics_tools/compute_cor_IGSR/filter_pop_IGSR.R /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed/${array_dir}
	echo "done $array_dir"
done