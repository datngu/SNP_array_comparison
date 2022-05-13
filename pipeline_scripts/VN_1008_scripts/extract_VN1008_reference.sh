#!/bin/bash
#SBATCH --job-name=GT_true_vn1008
#SBATCH --output=GT_true_vn1008.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=128000


cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008

extract_sample ()
{
    i=$1
    in_vcf=VN_1008.chr${i}.all.vcf.gz
    out_vcf=VN_1008.chr${i}.all.filled_tag.vcf.gz
    out_GT=chr${i}_GT.txt
    bcftools +fill-tags $in_vcf | bgzip > $out_vcf
    echo $'ID\tAF\tMAF\tAC' $(bcftools view -h $out_vcf | grep "^#CHROM" | cut -f10-) | sed -e 's/ /\t/g' > $out_GT
    bcftools query $out_vcf -f '%CHROM:%POS:%REF:%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/AC[\t%GT]\n' | sed 's/0\/0/0/g' | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | awk '($4 > 1) {print}' >> $out_GT
    bgzip $out_GT
    rm $in_vcf
    rm ${in_vcf}.tbi
    mv $out_vcf $in_vcf
    bcftools index -t $in_vcf
    echo "Done CHR $i"
}


for i in {1..22}
do
    extract_sample $i &
done
wait
