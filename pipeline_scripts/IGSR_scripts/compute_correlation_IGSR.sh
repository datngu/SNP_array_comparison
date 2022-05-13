#!/bin/bash
ref_file=$1
test_file_GT=$2
out_cor=$3


/path/to/slipt_files.R $ref_file $test_file_GT $out_cor 10

for i in {1..10}
do
  /path/to/compute_corr_part.R ${out_cor}_part${i}.Rdata &
done
wait


/path/to/merge_results.R ${out_cor}

rm ${out_cor}_part*.*

echo "DONE"

# /media/datn/data/1st_DC_PRS_array_project/compute_cor/pipeline_compute_correlation_KHV.sh /home/datn/Downloads/processed_GT_chr20_1013_KHV.txt.gz /home/datn/Downloads/chr20.dose.txt.gz chr20_corr.txt
