## shell
cd /media/datn/data/1KG/IGSR_h38_high_cov/VCF_sample_information/10_folds_id


for i in {1..10}
do
	cat batch_${i}_EAS.txt batch_${i}_SAS.txt batch_${i}_EUR.txt batch_${i}_AMR.txt batch_${i}_AFR.txt > /media/datn/data/1KG/IGSR_h38_high_cov/10_folds_VCF/ALL_pop_batch_${i}.txt
done






#######################3

# 2504 samples
cd /media/datn/data/1KG/IGSR_h38_high_cov/VCF_no_chr

cat /media/datn/data/1KG/IGSR_h38_high_cov/VCF_sample_information/supper_population_*.txt > 2504_NYGC_high_cov_samples.txt

for i in {1..22}
do
  bcftools view chr${i}_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -S 504_NYGC_high_cov_samples.txt | bcftools +fill-tags -o chr${i}_2504_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -Oz & 
done

for i in {1..22}
do
  bcftools index -t chr${i}_2504_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz & 
done



### BY POPULATIONS

cd /media/datn/data/1KG/IGSR_h38_high_cov/VCF_no_chr

pop=EUR
mkdir $pop
for i in {1..22}
do
  bcftools view chr${i}_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -S /media/datn/data/1KG/IGSR_h38_high_cov/VCF_sample_information/supper_population_${pop}.txt | bcftools +fill-tags -o ${pop}/chr${i}_${pop}_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -Oz &
done
wait

for i in {1..22}
do
  bcftools index -t ${pop}/chr${i}_${pop}_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz & 
done
wait


pop=EAS
mkdir $pop
for i in {1..22}
do
  bcftools view chr${i}_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -S /media/datn/data/1KG/IGSR_h38_high_cov/VCF_sample_information/supper_population_${pop}.txt | bcftools +fill-tags -o ${pop}/chr${i}_${pop}_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -Oz &
done
wait

for i in {1..22}
do
  bcftools index -t ${pop}/chr${i}_${pop}_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz & 
done
wait



pop=SAS
mkdir $pop
for i in {1..22}
do
  bcftools view chr${i}_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -S /media/datn/data/1KG/IGSR_h38_high_cov/VCF_sample_information/supper_population_${pop}.txt | bcftools +fill-tags -o ${pop}/chr${i}_${pop}_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -Oz &
done
wait

for i in {1..22}
do
  bcftools index -t ${pop}/chr${i}_${pop}_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz & 
done
wait


pop=AMR
mkdir $pop
for i in {1..22}
do
  bcftools view chr${i}_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -S /media/datn/data/1KG/IGSR_h38_high_cov/VCF_sample_information/supper_population_${pop}.txt | bcftools +fill-tags -o ${pop}/chr${i}_${pop}_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -Oz &
done
wait

for i in {1..22}
do
  bcftools index -t ${pop}/chr${i}_${pop}_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz & 
done
wait

pop=AFR
mkdir $pop
for i in {1..22}
do
  bcftools view chr${i}_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -S /media/datn/data/1KG/IGSR_h38_high_cov/VCF_sample_information/supper_population_${pop}.txt | bcftools +fill-tags -o ${pop}/chr${i}_${pop}_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz -Oz &
done
wait

for i in {1..22}
do
  bcftools index -t ${pop}/chr${i}_${pop}_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz & 
done
wait




##/dragennfs/area7/datnguyen/PRS_arrays_project/2504_NYGC_unrelated_vcf

cd /media/datn/data/1KG/IGSR_h38_high_cov/10_folds_VCF


extract_folds()
{
  i=$1
  mkdir batch_${i}
  for chr in {1..22}
  do
  	echo "doing chr${chr}; batch ${i} :D"
  	#mkdir batch_${i}
  	vcf=/media/datn/data/1KG/IGSR_h38_high_cov/VCF_no_chr/chr${chr}_2504_unrelated_NYGC_high_cov_phased_no_chr_filled_tag.vcf.gz
  	sample=ALL_pop_batch_${i}.txt
  	bcftools view $vcf -S ^${sample} | bcftools +fill-tags -o batch_${i}/chr${chr}_reference.vcf.gz -Oz
  	bcftools index -t batch_${i}/chr${chr}_reference.vcf.gz
  	bcftools view $vcf -S ${sample} | bcftools +fill-tags -o batch_${i}/chr${chr}_test_set.vcf.gz -Oz
  	bcftools index -t batch_${i}/chr${chr}_test_set.vcf.gz
  	echo "DONE chr${chr}; batch ${i} :D"
  done
}


for i in {1..10}
do
  extract_folds $i & 
done










#bcftools view -h /media/datn/data/1KG/IGSR_h38_high_cov/10_folds_VCF/batch_1/chr22_test_set.vcf.gz | grep "^#CHROM" | cut -f10- | wc

