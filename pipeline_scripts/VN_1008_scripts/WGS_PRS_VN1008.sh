#!/bin/bash
#SBATCH --job-name=VN1008_PRS
#SBATCH --output=VN1008_PRS.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=256000




mkdir PRS_WGS

get_PRS ()
{
    trait=$1
    pop=$2
    base_name=merged_all_autosome
    sumstat=/dragennfs/area7/datnguyen/bioinformatics_tools/sumstat/GIANT_${trait}.QCed.gz
    bedfiles=plink_WGS/merged_all_autosome.QC
    outfn=PRS_WGS/WGS_${pop}_${trait}
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
    --thread 8
}


cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS

get_PRS HEIGHT VNP &
get_PRS BMI VNP &

wait
# get /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR/PRS_WGS/*all_score /media/datn/data/1st_DC_PRS_array_project/PRS_codes/PRS_6_pop


#######################


#!/bin/bash
#SBATCH --job-name=t2d_BC_PRS
#SBATCH --output=t2d_BC_PRS.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=256000

get_PRS_binary ()
{
    trait=$1
    pop=$2
    base_name=merged_all_autosome
    sumstat=/dragennfs/area7/datnguyen/bioinformatics_tools/sumstat/sumstat_GWAS_catalog_${trait}.QCed.gz
    bedfiles=plink_WGS/merged_all_autosome.QC
    outfn=PRS_WGS/WGS_${pop}_${trait}
    #PRSice_linux
    /dragennfs/area7/datnguyen/bioinformatics_tools/PRSice_linux/PRSice_linux \
    --base ${sumstat} \
    --target ${bedfiles} \
    --out ${outfn} \
    --binary-target T \
    --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 --fastscore \
    --a1 hm_effect_allele \
    --a2 hm_other_allele \
    --beta  \
    --bp hm_pos \
    --chr hm_chrom \
    --pvalue p_value \
    --snp hm_rsid \
    --stat hm_beta \
    --clump-kb 250kb \
    --clump-p 1 \
    --clump-r2 0.1 \
    --ultra \
    --no-regress \
    --score sum \
    --thread 8
}

#trait=Type_2_diabetes
#trait=Breast_cancer

cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS

get_PRS_binary Type_2_diabetes VNP &
get_PRS_binary Breast_cancer VNP &
wait