#!/bin/bash
#SBATCH --job-name=10f
#SBATCH --output=10f.out
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=256000




export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/Minimac3Executable/bin:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/Minimac4/minimac4-1.0.2-Linux/bin:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/impute5_v1.1.4/impute5_v1.1.4:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/plink1.9:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/LmTag-0.0.0/bin:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/vcftools/src/cpp:$PATH
export PERL5LIB=/dragennfs/area7/datnguyen/bioinformatics_tools/vcftools/src/perl/
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/shapeit4-4.1.3/bin:$PATH
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/bcftools-1.10.2:$PATH
export BCFTOOLS_PLUGINS=/dragennfs/area7/datnguyen/bioinformatics_tools/bcftools-1.10.2/plugins
export PATH=/dragennfs/area7/datnguyen/bioinformatics_tools/TagIt/src:$PATH

# for i in {1..22}
# do
#   cp /home/datn/NAS/DatNT/phasing_vn1013/chr${i}_dbSNP151_no_Chr_PHASED_SHAPIT4_testPop_1026.gatk.norm.passVQSR.rm13Sample.vcf.gz /media/datn/data/1st_DC_PRS_array_project/1013_KHV_GATK_phased/chr${i}_GATK_1013_KHV.vcf.gz
# done

# for i in {1..22}
# do
#   bcftools index --tbi /media/datn/data/1st_DC_PRS_array_project/1013_KHV_GATK_phased/chr${i}_GATK_1013_KHV.vcf.gz &
# done


#put -r /media/datn/data/1st_DC_PRS_array_project/1013_KHV_GATK_phased /dragennfs/area8/datnguyen/PRS_arrays_project/
#put -r /media/datn/data/1st_DC_PRS_array_project/10_folds_VCF_KHV /dragennfs/area8/datnguyen/PRS_arrays_project/
# put /media/datn/data/1st_DC_PRS_array_project/VN_1008_scripts/10_folds_VCF_1008/*.txt /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/10_fold_VCF_1008

cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/10_fold_VCF_1008/



extract_ref_test()
{
  n=$1
  sample=batch_${n}.txt
  for i in {1..22}
  do
    bcftools view /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/VN_1008.chr${i}.all.vcf.gz -S ^${sample} -o batch_${n}/chr${i}_reference.vcf.gz -Oz &
    bcftools view /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/phased_1008/VN_1008.chr${i}.all.vcf.gz -S ${sample} -o batch_${n}/chr${i}_test_set.vcf.gz -Oz &
  done
  wait
}


for n in {1..10}
do
  mkdir -p batch_${n}
  extract_ref_test $n
done
wait

for n in {1..10}
do
  for fi in batch_${n}/*.vcf.gz
  do
    bcftools index -t -f $fi &
  done
  wait
done

# wait 15m

cat done!