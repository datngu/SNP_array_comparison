


# igsr 10f res
cd /dragennfs/area8/datnguyen/PRS_arrays_project

mkdir igsr_10f_correlation_res

array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'


cd /dragennfs/area8/datnguyen/PRS_arrays_project
for array_dir in $array_list
do
	gzip -f IGSR_imputed/${array_dir}/*correlation.txt &
	wait
done

#array_dir=infinium-global-diversity-array-v1.0
cd /dragennfs/area8/datnguyen/PRS_arrays_project

for array_dir in $array_list
do
	mkdir -p igsr_10f_correlation_res/${array_dir}
	cp IGSR_imputed/${array_dir}/*correlation.txt.gz igsr_10f_correlation_res/${array_dir}
done



get -r /dragennfs/area8/datnguyen/PRS_arrays_project/igsr_10f_correlation_res /media/datn/data/1st_DC_PRS_array_project/



array_dir=infinium-global-diversity-array-v1.0
cd /dragennfs/area8/datnguyen/PRS_arrays_project
cp IGSR_imputed/${array_dir}/chr9*correlation.txt.gz igsr_10f_correlation_res/${array_dir}

get /dragennfs/area8/datnguyen/PRS_arrays_project/igsr_10f_correlation_res/infinium-global-diversity-array-v1.0/chr9*correlation.txt.gz /media/datn/data/1st_DC_PRS_array_project/igsr_10f_correlation_res/infinium-global-diversity-array-v1.0



#KHV_10f

cd /dragennfs/area8/datnguyen/PRS_arrays_project

mkdir khv_10f_correlation_res

array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'



#array_dir=infinium-global-diversity-array-v1.0
cd /dragennfs/area8/datnguyen/PRS_arrays_project

for array_dir in $array_list
do
	mkdir -p khv_10f_correlation_res/${array_dir}
	mv KHV_10_folds_imputed/${array_dir}/*correlation.txt.gz igsr_10f_correlation_res/${array_dir}
done





get -r /dragennfs/area8/datnguyen/PRS_arrays_project/khv_10f_correlation_res /media/datn/data/1st_DC_PRS_array_project/



## local

cd /media/datn/data/1st_DC_PRS_array_project/

array_list='Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_JAPONICA infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 japanese-screening-array-v1.0 multi-ethnic-eur-eas-sas-v1.0 multi-ethnic-global-v1.0 oncoarray-500k Axiom_PMDA Axiom_PMRA Axiom_UKB_WCSG chinese-genotyping-array-v1.0 cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-core-v1.2 infinium-global-diversity-array-v1.0 infinium-omni2.5.v1.5 infinium-omni5-v1.2 GenomeWideSNP_6.0'

for array_dir in $array_list
do
	mv khv_10f_correlation_res/${array_dir}/*correlation.txt.gz igsr_10f_correlation_res/${array_dir}
done


rm -r khv_10f_correlation_res





