#!/bin/bash
#SBATCH --job-name=plink_files
#SBATCH --output=plink_files.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=756000



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
process_array AMR &
process_array AFR &
wait