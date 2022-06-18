#!/bin/bash
#SBATCH --job-name=PRS
#SBATCH --output=PRS.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=512G




## PRS

#array_list="gsa_real_array gsa_simulated_array pmra_real_array pmra_simulated_array"

mkdir PRS_result

get_PRS ()
{
  array_dir=$1
  trait=$2
  pop=VNP
  base_name=merged_all_autosome
  sumstat=/dragennfs/area7/datnguyen/bioinformatics_tools/sumstat/GIANT_${trait}.QCed.gz
  bedfiles=${array_dir}/plink/${base_name}.QC
  #bedfiles=${array_dir}/plink_${array_dir}/${base_name}.QC
  outfn=PRS_result/${array_dir}_${pop}_${trait}
  #PRSice_linux
  /dragennfs/area7/datnguyen/bioinformatics_tools/PRSice_linux/PRSice_linux \
    --base ${sumstat} \
    --target ${bedfiles} \
    --out ${outfn} \
    --binary-target F \
    --bar-levels 0.00000005,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.2,0.3,0.5,1 --fastscore \
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



## imputation correlation


array_list="gsa_real_array gsa_simulated_array pmra_real_array pmra_simulated_array"



for array_dir in $array_list
do
  for i in {22..1}
  do
    bash /dragennfs/area7/datnguyen/bioinformatics_tools/compute_cor/pipeline_compute_correlation_KHV.sh /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/real_data/wgs_24_sample/ref_GT/chr${i}.txt.gz /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/real_data/${array_dir}/chr${i}_imputed.dose.txt.gz /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/real_data/${array_dir}/chr${i}_imputed.correlation.txt.gz &
  done
  wait
  echo "done $array_dir"
done

wait


# PRS


for array_dir in $array_list
do
  get_PRS $array_dir HEIGHT
  get_PRS $array_dir BMI
  get_PRS $array_dir Type_2_diabetes
done




# ### copy results
# get -r PRS_result /media/datn/data2gb/GitHub/SNP_array_comparison_full_data/real_data



# # remote

# get real_24_sample/*cor*gz  /media/datn/data2gb/GitHub/SNP_array_comparison_full_data/real_data/real_24_sample

# get simulated_24_sample/*cor*gz  /media/datn/data2gb/GitHub/SNP_array_comparison_full_data/real_data/simulated_24_sample


