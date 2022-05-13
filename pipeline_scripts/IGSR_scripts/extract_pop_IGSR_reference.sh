#!/bin/bash
#SBATCH --job-name=GT_true_IGSR
#SBATCH --output=GT_true_IGSR.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=128000


cd /dragennfs/area8/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf
pop=EUR


extract_sample ()
{
    i=$1
    pop=$2
    in_vcf=chr${i}_2504_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz
    mkdir ${pop}_process_GT
    out_vcf=${pop}_process_GT/chr${i}_filled_tag.vcf.gz
    out_GT=${pop}_process_GT/chr${i}_GT.txt
    sample=./VCF_sample_information/supper_population_${pop}.txt
    bcftools view $in_vcf -S $sample | bcftools +fill-tags | bgzip > $out_vcf
    echo $'ID\tAF\tMAF\tAC' $(bcftools view -h $out_vcf | grep "^#CHROM" | cut -f10-) | sed -e 's/ /\t/g' > $out_GT
    bcftools query $out_vcf -f '%CHROM:%POS:%REF:%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/AC[\t%GT]\n' | sed 's/0\/0/0/g' | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | awk '($4 > 1) {print}' >> $out_GT
    bgzip $out_GT
    echo "Done CHR $i, POP $pop"
}


for i in {1..22}
do
    extract_sample $i EUR &
done
wait
for i in {1..22}
do
    extract_sample $i EAS &
done
wait
for i in {1..22}
do
    extract_sample $i SAS &
done
wait
for i in {1..22}
do
    extract_sample $i AMR &
done
wait
for i in {1..22}
do
    extract_sample $i AFR &
done
wait
