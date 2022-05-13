#!/bin/bash
#SBATCH --job-name=t2d_all
#SBATCH --output=t2d_all.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=90
#SBATCH --mem=512000


# WGS VN1008
#mkdir PRS_all_results


get_PRS ()
{
    trait=$1
    pop=$2
    array=$3
    bedfiles=$4
    #sumstat=$5
    #base_name=merged_all_autosome
    sumstat=/dragennfs/area7/datnguyen/bioinformatics_tools/sumstat/GIANT_${trait}.QCed.gz
    #bedfiles=plink_WGS/merged_all_autosome.QC
    outfn=${array}_${pop}_${trait}
    #PRSice_linux
    /dragennfs/area7/datnguyen/bioinformatics_tools/PRSice_linux/PRSice_linux \
    --base ${sumstat} \
    --target ${bedfiles} \
    --out ${outfn} \
    --binary-target F \
    --bar-levels 0.001,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1 --fastscore \
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
    --print-snp \
    --thread 6
}


cd /dragennfs/area8/datnguyen/PRS_arrays_project/PRS_16_cutoffs

# array=Axiom_GW_ASI
# prs vn1008 WGS
pop=VNP
array=WGS
trait=HEIGHT
get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS/plink_WGS/merged_all_autosome.QC &
trait=BMI
get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS/plink_WGS/merged_all_autosome.QC &
trait=Type_2_diabetes
get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS/plink_WGS/merged_all_autosome.QC &
wait

array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

for array in $array_list
do
    pop=VNP
    trait=HEIGHT
    get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008/${array}/plink_${array}/merged_all_autosome.QC &
    trait=BMI
    get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008/${array}/plink_${array}/merged_all_autosome.QC &
    trait=Type_2_diabetes
    get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008/${array}/plink_${array}/merged_all_autosome.QC &
    wait
done

# prs IGSR

# WGS

get_PRS_IGSR_by_pop_WGS ()
{
  pop=$1
  array=WGS
  trait=HEIGHT
  get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR/pop_${pop}/merged_all_autosome_${pop}.QC &
  trait=BMI
  get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR/pop_${pop}/merged_all_autosome_${pop}.QC &
  trait=Type_2_diabetes
  get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR/pop_${pop}/merged_all_autosome_${pop}.QC &
  wait
}

for pop in EAS EUR SAS AMR AFR
do
  get_PRS_IGSR_by_pop_WGS $pop &
done

wait

## arrays

get_PRS_IGSR_by_pop_ARRAY ()
{
  pop=$1
  array=$2
  trait=HEIGHT
  get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed/${array}/plink_${array}/merged_all_autosome_${pop}.QC &
  trait=BMI
  get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed/${array}/plink_${array}/merged_all_autosome_${pop}.QC &
  trait=Type_2_diabetes
  get_PRS ${trait} ${pop} ${array} /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed/${array}/plink_${array}/merged_all_autosome_${pop}.QC &
  wait
}


array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

for array in $array_list
do
  for pop in EAS EUR SAS AMR AFR
  do
    get_PRS_IGSR_by_pop_ARRAY $pop $array &
  done
  wait
done

wait