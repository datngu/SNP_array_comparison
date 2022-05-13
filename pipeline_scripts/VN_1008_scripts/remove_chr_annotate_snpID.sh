#!/bin/bash
#SBATCH --job-name=remove_chr
#SBATCH --output=remove_chr.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=128000


cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all
mkdir phased_1008
cd phased_1008


remove_fun ()
{
	i=$1
	dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
	file=/dragennfs/area15/Ref_imputation/VN_1008/VN_1008.chr${i}.all.vcf.gz
	zcat $file | sed 's/chr//g' | bgzip > VN_1008.chr${i}.all.tem.vcf.gz
	bcftools index -t VN_1008.chr${i}.all.tem.vcf.gz
	bcftools annotate -a $dbsnp151 -c ID VN_1008.chr${i}.all.tem.vcf.gz | bgzip > VN_1008.chr${i}.all.vcf.gz
	bcftools index -t VN_1008.chr${i}.all.vcf.gz
	rm VN_1008.chr${i}.all.tem.vcf.gz*
	echo "DONE VN_1008.chr${i}.all.vcf.gz!"
}


for i in {1..22}
do
	remove_fun $i &
done
wait