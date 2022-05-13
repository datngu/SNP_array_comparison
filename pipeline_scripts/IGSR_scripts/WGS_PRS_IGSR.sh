#!/bin/bash
#SBATCH --job-name=plink_files
#SBATCH --output=plink_files.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=256000


#cd /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf




# annotate_IGSR ()
# {
#     i=$1
#     out_dir=/dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR
#     dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
#     file=chr${i}_2504_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz
#     bcftools annotate -a $dbsnp151 -c ID $file | bgzip > ${out_dir}/chr${i}.2504_IGSR_annotated_dbSNP151.vcf.gz
# }


# for i in {1..22}
# do
#    annotate_IGSR $i &
# done
# wait




# cd /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR

# for i in {1..22}
# do
#     bcftools index -t chr${i}.2504_IGSR_annotated_dbSNP151.vcf.gz &
# done
# wait

# extract_sample ()
# {
#     i=$1
#     pop=$2
#     in_vcf=chr${i}.2504_IGSR_annotated_dbSNP151.vcf.gz
#     mkdir pop_${pop}
#     out_vcf=pop_${pop}/chr${i}_annotated_dbSNP151.vcf.gz
#     sample=/dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf/VCF_sample_information/supper_population_${pop}.txt
#     bcftools view $in_vcf -S $sample | bgzip > $out_vcf
#     echo "Done CHR $i, POP $pop"
# }





# for i in {1..10}
# do
#     extract_sample $i EUR &
#     extract_sample $i EAS &
#     extract_sample $i SAS &
#     extract_sample $i AMR &
#     extract_sample $i AFR &
# done
# wait

# for i in {11..22}
# do
#     extract_sample $i EUR &
#     extract_sample $i EAS &
#     extract_sample $i SAS &
#     extract_sample $i AMR &
#     extract_sample $i AFR &
# done
# wait



cd /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR

process_array ()
{
    pop=$1
    process_fam=/dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed/${pop}.fam
    base_name=merged_all_autosome_${pop}
    #dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
    ls pop_${pop}/chr*_annotated_dbSNP151.vcf.gz > pop_${pop}/file.txt
    bcftools concat -n -f pop_${pop}/file.txt -Oz -o pop_${pop}/${base_name}_annotated_dbSNP151.vcf.gz.vcf.gz
    #bcftools index -t ${array_dir}/plink_${array_dir}/${base_name}.vcf.gz
    #bcftools annotate -a $dbsnp151 -c ID ${array_dir}/plink_${array_dir}/${base_name}.vcf.gz | bgzip > ${array_dir}/plink_${array_dir}/${base_name}_db151.vcf.gz

    plink --vcf pop_${pop}/${base_name}_annotated_dbSNP151.vcf.gz.vcf.gz \
      --make-bed  --const-fid --out pop_${pop}/${base_name}_raw \
      --threads 2 \
      --memory 128000

    
    # copy processed fam with FID and Sex
    cat ${process_fam} > pop_${pop}/${base_name}_raw.fam

    plink \
    --bfile pop_${pop}/${base_name}_raw \
    --maf 0.0001 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --memory 128000 \
    --out pop_${pop}/${base_name}.QC

    r_code='
    #!/usr/bin/env Rscript\n
    args = commandArgs(trailingOnly=TRUE)\n
    options(stringsAsFactors = FALSE)\n
    require(data.table)\n
    in_file = args[1]\n
    out_file = args[2]\n
    snp = fread(in_file, header = F)\n
    pick = snp$V1 == "."\n
    snp = snp[!pick,]\n
    pick = duplicated(snp$V1)\n
    snp_filter = snp[!pick,]\n
    fwrite(snp_filter, file = out_file, sep = "\t")\n
    cat("DONE merging files")\n
    '
    echo -e $r_code > filter_duplicated.R
    chmod 777 ./filter_duplicated.R
    ./filter_duplicated.R pop_${pop}/${base_name}.QC.snplist pop_${pop}/${base_name}.QC.snplist.nodup 
    

    plink \
        --bfile pop_${pop}/${base_name}_raw \
        --threads 2 \
        --make-bed \
        --keep pop_${pop}/${base_name}.QC.fam \
        --out pop_${pop}/${base_name}.QC \
        --extract pop_${pop}/${base_name}.QC.snplist.nodup \
        --memory 128000

    rm pop_${pop}/${base_name}_raw*
    #rm pop_${pop}/${base_name}.vcf.gz
    rm pop_${pop}/${base_name}.QC.snplis*
    #rm ${array_dir}/plink_${array_dir}/${base_name}_db151.vcf.gz
}


process_array EAS &
process_array EUR &
process_array SAS &

wait
process_array AMR &
process_array AFR &
wait


# get PRS

mkdir PRS_WGS

get_PRS ()
{
    trait=$1
    pop=$2
    base_name=merged_all_autosome_${pop}
    sumstat=/dragennfs/area7/datnguyen/bioinformatics_tools/sumstat/GIANT_${trait}.QCed.gz
    bedfiles=pop_${pop}/${base_name}.QC
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
    --thread 1
}

cd /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR

get_PRS HEIGHT EUR &
get_PRS HEIGHT EAS &
get_PRS HEIGHT SAS &
get_PRS HEIGHT AMR &
get_PRS HEIGHT AFR &
get_PRS BMI EUR &
get_PRS BMI EAS &
get_PRS BMI SAS &
get_PRS BMI AMR &
get_PRS BMI AFR &

wait





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
    base_name=merged_all_autosome_${pop}
    sumstat=/dragennfs/area7/datnguyen/bioinformatics_tools/sumstat/sumstat_GWAS_catalog_${trait}.QCed.gz
    bedfiles=pop_${pop}/${base_name}.QC
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

trait=Type_2_diabetes
trait=Breast_cancer

cd /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR

get_PRS_binary Type_2_diabetes EUR &
get_PRS_binary Type_2_diabetes EAS &
get_PRS_binary Type_2_diabetes SAS &
get_PRS_binary Type_2_diabetes AMR &
get_PRS_binary Type_2_diabetes AFR &
get_PRS_binary Breast_cancer EUR &
get_PRS_binary Breast_cancer EAS &
get_PRS_binary Breast_cancer SAS &
get_PRS_binary Breast_cancer AMR &
get_PRS_binary Breast_cancer AFR &

wait


########### download all PRS

get /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/WGS/PRS_WGS/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores
get /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/imputed_1008/PRS_result/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores

get /dragennfs/area8/datnguyen/PRS_arrays_project/WGS_PRS_IGSR/PRS_WGS/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores
get /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed/PRS_result_all_pop/*all_score /media/datn/data2gb/GitHub/SNP_array_comparsion/data/PRS_scores