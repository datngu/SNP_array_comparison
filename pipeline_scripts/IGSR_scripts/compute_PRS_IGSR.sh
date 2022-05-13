#!/bin/bash
#SBATCH --job-name=PRS_scoring
#SBATCH --output=PRS_scoring.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=512000


#put -r /home/datn/bioinformatics_tools/PRSice_linux /dragennfs/area7/datnguyen/bioinformatics_tools
#put /media/datn/data/1st_DC_PRS_array_project/PRS_codes/sumstat/GIANT_BMI.QCed.gz /dragennfs/area7/datnguyen/bioinformatics_tools/sumstat
#put /media/datn/data/1st_DC_PRS_array_project/PRS_codes/sumstat/GIANT_HEIGHT.QCed.gz /dragennfs/area7/datnguyen/bioinformatics_tools/sumstat

#sumstat=/dragennfs/area7/datnguyen/bioinformatics_tools/sumstat/GIANT_${trait}.QCed.gz

cd /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed


array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

trait=BMI
trait=HEIGHT

mkdir PRS_result_all_pop

get_PRS ()
{
	array_dir=$1
	trait=$2
	pop=$3
	base_name=merged_all_autosome_${pop}
	sumstat=/dragennfs/area7/datnguyen/bioinformatics_tools/sumstat/GIANT_${trait}.QCed.gz
	bedfiles=${array_dir}/plink_${array_dir}/${base_name}.QC
	outfn=PRS_result_all_pop/${array_dir}_${pop}_${trait}
	#PRSice_linux
	/dragennfs/area7/datnguyen/bioinformatics_tools/PRSice_linux/PRSice_linux \
    --base ${sumstat} \
    --target ${bedfiles} \
    --out ${outfn} \
    --binary-target F \
    --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 --fastscore \
    --a1 A1 \
    --a2 A2 \
    --beta  \
    --bp BP \
    --chr CHR \
    --pvalue P \
    --snp SNP \
    --stat BETA \
    --clump-kb 250kb \
    --clump-p 1 \
    --clump-r2 0.1 \
    --ultra \
    --no-regress \
    --score sum \
    --thread 1
}

get_PRS_5_pop ()
{
	array_dir=$1
	get_PRS $array_dir HEIGHT EUR &
	get_PRS $array_dir HEIGHT EAS &
	get_PRS $array_dir HEIGHT SAS &
	get_PRS $array_dir HEIGHT AMR &
	get_PRS $array_dir HEIGHT AFR &
	get_PRS $array_dir BMI EUR &
	get_PRS $array_dir BMI EAS &
	get_PRS $array_dir BMI SAS &
	get_PRS $array_dir BMI AMR &
	get_PRS $array_dir BMI AFR &
	wait
}



#array_dir=Axiom_GW_ASI
#get_PRS_5_pop $array_dir HEIGHT



for array_dir in $array_list
do
	get_PRS_5_pop $array_dir
done


wait

