#!/bin/bash
#SBATCH --job-name=chr17-22
#SBATCH --output=chr17-22.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=500000



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



# put -r /media/datn/data/1st_DC_PRS_array_project/2504_NYGC_unrelated_vcf /dragennfs/area7/datnguyen/PRS_arrays_project/
# put -r /media/datn/data/1st_DC_PRS_array_project/10_folds_VCF /dragennfs/area7/datnguyen/PRS_arrays_project
# put -r /media/datn/data/1st_DC_PRS_array_project/10_folds_VCF/batch_3 /dragennfs/area7/datnguyen/PRS_arrays_project/10_folds_VCF
# put -r /media/datn/data/1st_DC_PRS_array_project/10_folds_VCF/batch_4 /dragennfs/area7/datnguyen/PRS_arrays_project/10_folds_VCF
# put -r /media/datn/data/1st_DC_PRS_array_project/10_folds_VCF/batch_5 /dragennfs/area7/datnguyen/PRS_arrays_project/10_folds_VCF

# put -r /media/datn/data/1st_DC_PRS_array_project/10_folds_VCF/batch_6 /dragennfs/area7/datnguyen/PRS_arrays_project/10_folds_VCF
# put -r /media/datn/data/1st_DC_PRS_array_project/10_folds_VCF/batch_7 /dragennfs/area7/datnguyen/PRS_arrays_project/10_folds_VCF
# put -r /media/datn/data/1st_DC_PRS_array_project/10_folds_VCF/batch_8 /dragennfs/area7/datnguyen/PRS_arrays_project/10_folds_VCF
# put -r /media/datn/data/1st_DC_PRS_array_project/10_folds_VCF/batch_9 /dragennfs/area7/datnguyen/PRS_arrays_project/10_folds_VCF
# put -r /media/datn/data/1st_DC_PRS_array_project/10_folds_VCF/batch_10 /dragennfs/area7/datnguyen/PRS_arrays_project/10_folds_VCF



cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/10_fold_VCF_1008

for n in {1..10}
  do
    mkdir -p m3vcf_batch_${n}
done

#mkdir m3vcf
#mkdir 2504_NYGC_unrelated_vcf


gen_ref ()
{
  i=$1
  echo "processing chr ${i}"
  for n in {1..10}
  do
    Minimac3 --refHaps batch_${n}/chr${i}_reference.vcf.gz --processReference --prefix m3vcf_batch_${n}/chr${i}_reference --log m3vcf_batch_${n}/chr${i}_log &
  done
  wait
  echo "DONE processing chr ${i} all batch"
}

# process chr 1-8
for i in {17..22} {1..6}
do
  gen_ref $i &
done

wait

echo "done chr 17:22 - 1:6"



# process chr 1-8
for i in {7..16}
do
  gen_ref $i &
done

wait

echo "done chr 7:16"

