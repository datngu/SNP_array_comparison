#!/bin/bash
#SBATCH --job-name=Plink_files_IGSR
#SBATCH --output=Plink_files_IGSR.log
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=120
#SBATCH --mem=256000


cd /dragennfs/area8/datnguyen/PRS_arrays_project/IGSR_imputed


array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

#array_list='Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

array_list1='Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG'

array_list2='chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

# array_dir=Axiom_GW_ASI

#chr2_AFR_imputed.vcf.gz
annotate_dbSNP ()
{
	array_dir=$1
	pop=$2
	base_name=merged_all_autosome_${pop}
	dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
	ls ${array_dir}/*_${pop}_imputed.vcf.gz > ${array_dir}/file_${pop}.txt
	mkdir ${array_dir}/plink_${array_dir}
	bcftools concat -n -f ${array_dir}/file_${pop}.txt -Oz -o ${array_dir}/plink_${array_dir}/${base_name}.vcf.gz
	bcftools index -t ${array_dir}/plink_${array_dir}/${base_name}.vcf.gz
	bcftools annotate -a $dbsnp151 -c ID ${array_dir}/plink_${array_dir}/${base_name}.vcf.gz | bgzip > ${array_dir}/plink_${array_dir}/${base_name}_db151.vcf.gz
}

for array_dir in $array_list1
do
	annotate_dbSNP $array_dir EUR &
	annotate_dbSNP $array_dir EAS &
	annotate_dbSNP $array_dir SAS &
	annotate_dbSNP $array_dir AMR &
	annotate_dbSNP $array_dir AFR &
done

wait


for array_dir in $array_list2
do
	annotate_dbSNP $array_dir EUR &
	annotate_dbSNP $array_dir EAS &
	annotate_dbSNP $array_dir SAS &
	annotate_dbSNP $array_dir AMR &
	annotate_dbSNP $array_dir AFR &
done

wait


	# plink --vcf ${array_dir}/plink_${array_dir}/${base_name}.vcf.gz \
 #      --make-bed  --const-fid --out ${array_dir}/plink_${array_dir}/${base_name}_raw \
 #      --threads 2 \
 #      --memory 128000


process_array ()
{
	array_dir=$1
	pop=$2
	process_fam=${pop}.fam
	base_name=merged_all_autosome_${pop}
	#dbsnp151=/dragennfs/area7/datnguyen/bioinformatics_tools/dbSNP/dbsnp151.vcf.gz
	#ls ${array_dir}/*.vcf.gz > ${array_dir}/file.txt
	#mkdir ${array_dir}/plink_${array_dir}
	#bcftools concat -n -f ${array_dir}/file.txt -Oz -o ${array_dir}/plink_${array_dir}/${base_name}.vcf.gz
	#bcftools index -t ${array_dir}/plink_${array_dir}/${base_name}.vcf.gz
	#bcftools annotate -a $dbsnp151 -c ID ${array_dir}/plink_${array_dir}/${base_name}.vcf.gz | bgzip > ${array_dir}/plink_${array_dir}/${base_name}_db151.vcf.gz

	plink --vcf ${array_dir}/plink_${array_dir}/${base_name}_db151.vcf.gz \
      --make-bed  --const-fid --out ${array_dir}/plink_${array_dir}/${base_name}_raw \
      --threads 2 \
      --memory 128000

    
    # copy processed fam with FID and Sex
    cat ${process_fam} > ${array_dir}/plink_${array_dir}/${base_name}_raw.fam

    plink \
    --bfile ${array_dir}/plink_${array_dir}/${base_name}_raw \
    --maf 0.0001 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --memory 128000 \
    --out ${array_dir}/plink_${array_dir}/${base_name}.QC

    r_code='
	#!/usr/bin/env Rscript\n
	args = commandArgs(trailingOnly=TRUE)\n
	options(stringsAsFactors = FALSE)\n
	require(data.table)\n
	in_file = args[1]\n
	out_file = args[2]\n
	snp = fread(in_file, header = F)\n
	pick = snp$V1 == "."\n
	snp = snp[!pick,]\n
	pick = duplicated(snp$V1)\n
	snp_filter = snp[!pick,]\n
	fwrite(snp_filter, file = out_file, sep = "\t")\n
	cat("DONE merging files")\n
	'
	echo -e $r_code > filter_duplicated.R
	chmod 777 ./filter_duplicated.R
	./filter_duplicated.R ${array_dir}/plink_${array_dir}/${base_name}.QC.snplist ${array_dir}/plink_${array_dir}/${base_name}.QC.snplist.nodup	
	

	plink \
	    --bfile ${array_dir}/plink_${array_dir}/${base_name}_raw \
	    --threads 2 \
	    --make-bed \
	    --keep ${array_dir}/plink_${array_dir}/${base_name}.QC.fam \
	    --out ${array_dir}/plink_${array_dir}/${base_name}.QC \
	    --extract ${array_dir}/plink_${array_dir}/${base_name}.QC.snplist.nodup \
	    --memory 128000

	rm ${array_dir}/plink_${array_dir}/${base_name}_raw*
	rm ${array_dir}/plink_${array_dir}/${base_name}.vcf.gz
	rm ${array_dir}/plink_${array_dir}/${base_name}.QC.snplis*
	#rm ${array_dir}/plink_${array_dir}/${base_name}_db151.vcf.gz
}


for array_dir in $array_list
do
	process_array $array_dir EUR
	process_array $array_dir EAS
	process_array $array_dir SAS
	process_array $array_dir AMR
	process_array $array_dir AFR
	echo "DONE $array_dir"
done
echo "DONE ALL"