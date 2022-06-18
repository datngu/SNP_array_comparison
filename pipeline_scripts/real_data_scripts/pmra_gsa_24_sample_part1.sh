#!/bin/bash
#SBATCH --job-name=processing_data
#SBATCH --output=processing_data.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=512G



cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/real_data


# index
bcftools index -t pmra_real.vcf.gz
bcftools index -t gsa_real.vcf.gz



# extract samples
all_vcf=/dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/VN_1008.chr10.all.vcf.gz
gsa_vcf=gsa_real.vcf.gz
pmra_vcf=pmra_real.vcf.gz

bcftools view -h $pmra_vcf | grep "^#CHROM" | cut -f10- | sed 's/\t/\n/g' > all_sample.txt

bcftools view -h $gsa_vcf | grep "^#CHROM" | cut -f10- | sed 's/\t/\n/g' > gsa_sample.txt

grep -Fxf all_sample.txt gsa_sample.txt > shared_sample.txt




mkdir reference_set

#####################################################################################
# process reference

sample=shared_sample.txt

for i in {1..22}
do
    bcftools view /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/VN_1008.chr${i}.all.vcf.gz -S ^${sample} -o reference_set/chr${i}.vcf.gz -Oz &ß
done

wait


for i in {1..22}
do
    bcftools index -t reference_set/chr${i}.vcf.gz &
done
wait


mkdir reference_set/m3vcf
# index m3vcf
for i in {1..22}
do
  Minimac3 --refHaps reference_set/chr${i}.vcf.gz --processReference --prefix reference_set/m3vcf/chr${i} --log reference_set/m3vcf/chr${i}_log.txt &
done

wait


#####################################################################################
# process ground truth
mkdir wgs_24_sample

sample=shared_sample.txt
for i in {1..22}
do
    bcftools view /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/VN_1008.chr${i}.all.vcf.gz -S ${sample} -o wgs_24_sample/chr${i}.vcf.gz -Oz &
done

wait



process_ref ()
{
  vcf=$1
  out_GT=$2
  echo $'ID\tAF\tMAF\tAC' $(bcftools view -h $vcf | grep "^#CHROM" | cut -f10-) | sed -e 's/ /\t/g' > $out_GT
  bcftools query $vcf -f '%CHROM:%POS:%REF:%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/AC[\t%GT]\n' | sed 's/0\/0/0/g' | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | awk '($4 > 1) {print}' >> $out_GT
}

mkdir wgs_24_sample/ref_GT
for i in {1..22}
do
  process_ref wgs_24_sample/chr${i}.vcf.gz wgs_24_sample/ref_GT/chr${i}.txt &
done
wait

for i in {1..22}
do
  bgzip wgs_24_sample/ref_GT/chr${i}.txt &
done
wait


#############################
## GSA
# real genotype

mkdir gsa_real_array

sample=shared_sample.txt
for i in {1..22}
do
    bcftools view $gsa_vcf -r $i -S ${sample} -o gsa_real_array/chr${i}.vcf.gz -Oz &
done

wait


for i in {1..22}
do
    bcftools index -t gsa_real_array/chr${i}.vcf.gz &
done

wait


# simulated genotype

mkdir gsa_simulated_array

# copy annotation
# put /media/datn/data/1st_DC_PRS_array_project/array_annotation/annotation_all_hg38/global-screening-array-v.3/* /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/real_data/gsa_simulated_array

## simulate non-phasing data
#sample=shared_sample.txt
for i in {1..22}
do
  vcf_wgs=wgs_24_sample/chr${i}.vcf.gz
  position=gsa_simulated_array/chr${i}.txt
  vcftools --gzvcf  $vcf_wgs --positions $position --recode --stdout | sed 's/0|0/0\/0/g' | sed 's/0|1/0\/1/g' | sed 's/1|0/1\/0/g' | sed 's/1|1/1\/1/g' | bgzip > gsa_simulated_array/chr${i}.vcf.gz
  bcftools index -t gsa_simulated_array/chr${i}.vcf.gz
done


#############################
## PMRA
# real genotype

mkdir pmra_real_array

sample=shared_sample.txt
for i in {1..22}
do
    bcftools view $pmra_vcf -r $i -S ${sample} -o pmra_real_array/chr${i}.vcf.gz -Oz &
done

wait


for i in {1..22}
do
    bcftools index -t pmra_real_array/chr${i}.vcf.gz &
done

wait


# simulated genotype

mkdir pmra_simulated_array

# copy annotation
# put /media/datn/data/1st_DC_PRS_array_project/array_annotation/annotation_all_hg38/Axiom_PMRA/* /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/real_data/pmra_simulated_array

## simulate non-phasing data
#sample=shared_sample.txt
for i in {1..22}
do
  vcf_wgs=wgs_24_sample/chr${i}.vcf.gz
  position=pmra_simulated_array/chr${i}.txt
  vcftools --gzvcf  $vcf_wgs --positions $position --recode --stdout | sed 's/0|0/0\/0/g' | sed 's/0|1/0\/1/g' | sed 's/1|0/1\/0/g' | sed 's/1|1/1\/1/g' | bgzip > pmra_simulated_array/chr${i}.vcf.gz
  bcftools index -t pmra_simulated_array/chr${i}.vcf.gz
done

