
# local 

mkdir /media/datn/data2gb/GitHub/SNP_array_comparison/data/real_data/gsa_real_array
mkdir /media/datn/data2gb/GitHub/SNP_array_comparison/data/real_data/gsa_simulated_array
mkdir /media/datn/data2gb/GitHub/SNP_array_comparison/data/real_data/pmra_real_array
mkdir /media/datn/data2gb/GitHub/SNP_array_comparison/data/real_data/pmra_simulated_array



# ### copy results
cd /dragennfs/area8/datnguyen/PRS_arrays_project/VN_1008_all/real_data

###
get -r PRS_result /media/datn/data2gb/GitHub/SNP_array_comparison/data/real_data

get -r gsa_real_array/*cor*gz  /media/datn/data2gb/GitHub/SNP_array_comparison/data/real_data/gsa_real_array
get -r gsa_simulated_array/*cor*gz  /media/datn/data2gb/GitHub/SNP_array_comparison/data/real_data/gsa_simulated_array
get -r pmra_real_array/*cor*gz  /media/datn/data2gb/GitHub/SNP_array_comparison/data/real_data/pmra_real_array
get -r pmra_simulated_array/*cor*gz  /media/datn/data2gb/GitHub/SNP_array_comparison/data/real_data/pmra_simulated_array




# # remote

# get real_24_sample/*cor*gz  /media/datn/data2gb/GitHub/SNP_array_comparison_full_data/real_data/real_24_sample

# get simulated_24_sample/*cor*gz  /media/datn/data2gb/GitHub/SNP_array_comparison_full_data/real_data/simulated_24_sample

