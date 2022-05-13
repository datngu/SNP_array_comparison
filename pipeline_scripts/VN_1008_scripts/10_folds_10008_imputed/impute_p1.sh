#!/bin/bash
#SBATCH --job-name=imp_l2-10f
#SBATCH --output=imp_l2-10f.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=80
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



array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'


#array_list1='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k'

array_list2='Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008


for array_dir in $array_list2
do
	nohup bash impute_by_array.sh $array_dir > $array_dir/impute_log.txt &
done

wait
