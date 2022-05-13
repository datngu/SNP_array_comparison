# download array annotations
cd /media/datn/data/1st_DC_PRS_array_project/array_annotation/zipped_csv

wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/infinium-psycharray/v1-3/infinium-psycharray-24-v1-3-a1-manifest-file-csv.zip -O infinium-psycharray-v1.3.zip
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/omnizhonghua-8/v1-4/infinium-omnizhonghua-8-v1-4-manifest-file-csv.zip -O infinium-omnizhonghua-v1.4.zip
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/infinium-omni5-4/v1-2/infinium-omni5-4-v1-2-a2-manifest-file-csv.zip -O infinium-omni5-v1.2.zip
wget https://webdata.illumina.com/downloads/productfiles/humanomni25/v1-5/infinium-omni2-5-8v1-5-a1-manifest-file-csv.zip -O infinium-omni2.5.v1.5.zip
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/multiethnic-global-8/v1-0/build38/multi-ethnic-global-8-d2-csv.zip -O multi-ethnic-global-v1.0.zip
wget https://sapac.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/japanese-screening-array/japanese-screening-array-24v1-0_B2-manifest-file-csv.zip -O japanese-screening-array-v1.0.zip
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/chinese-genotyping-array/chinese-genotyping-array-24v1-0_A1-manifest-file-csv.zip -O chinese-genotyping-array-v1.0.zip
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/infiniumexome24/v1-0/infinium-exome-24-v1-0-manifest-file-csv.zip -O infinium-exome-v1.0.zip
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/v3-0/GSA-24v3-0-A2-manifest-file-csv.zip -O global-screening-array-v.3.zip
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/infiniumcore-24/v1-2/infiniumcore-24-v1-2-a1-manifest-file-csv.zip -O infinium-core-v1.2.zip
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/cytosnp/v2-1/humancytosnp-12-v2-1-ns550-a1-manifest-file-csv.zip -O human-cytosnp-12-v2.1.zip
wget https://webdata.illumina.com/downloads/productfiles/cytosnp-850k/v1-2/cytosnp-850k-v1-2-ns550-d2-manifest-file-csv.zip -O cytosnp-850k-v1.2.zip
get ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/multiethnic-amr-afr-8/v1-0/multi-ethnic-amr-afr-8-v1-0-a1-manifest-file-csv.zip -O multi-ethnic-amr-afr-v1.0.zip
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/multiethnic-eur-eas-sas-8/v1-0/multi-ethnic-eur-eas-sas-8v-1-0-a1-manifest-file-csv.zip -O multi-ethnic-eur-eas-sas-v1.0.zip
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/multiethnic-global-8/v1-0/build38/multi-ethnic-global-8-d2-csv.zip -O multi-ethnic-global-v1.0.zip
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/Downloads/ProductFiles/OncoArray-500K/oncoarray-500k-manifest-file-c-csv.zip -O oncoarray-500k.zip
wget https://webdata.illumina.com/downloads/productfiles/global-diversity-array/infinium-global-diversity-array-8-v1-0-D2-manifest-file-csv.zip -O infinium-global-diversity-array-v1.0.zip
### affymextrix
cd /media/datn/data/1st_DC_PRS_array_project/array_annotation/zipped_csv

wget http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/GenomeWideSNP_6.na35.annot.csv.zip -O GenomeWideSNP_6.0.zip
wget https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/Axiom_PMDA.na36.r7.a8.annot.csv.zip -O Axiom_PMDA.zip
wget https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/Axiom_PMRA.na35.r3.a1.annot.csv.zip -O Axiom_PMRA.zip
wget https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/Axiom_UKB_WCSG.na35.annot.csv.zip -O Axiom_UKB_WCSG.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/Axiom_GW_ASI_SNP.na35.annot.csv.zip -O Axiom_GW_ASI.zip
wget https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/Axiom_GW_CHB2-na35-annot-csv.zip -O Axiom_GW_CHB.zip
wget https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/Axiom_GW_EUR-na35-annot-csv.zip -O Axiom_GW_EUR.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/Axiom_GW_PanAFR.na35.annot.csv.zip -O Axiom_GW_PanAFR.zip
wget https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/Axiom_JAPONICA.na36.r1.a1.annot.csv.zip -O Axiom_JAPONICA


# AXIOM
cd /media/datn/data/1st_DC_PRS_array_project/array_annotation/annotaion_csv/axiom
chain_file="/media/datn/data/cross_map_lift_over/GRCh37_to_GRCh38.chain.gz"
#hg19

for array in Axiom_GW_ASI Axiom_GW_CHB Axiom_GW_EUR Axiom_GW_PanAFR Axiom_UKB_WCSG GenomeWideSNP_6.0 Axiom_PMRA
do
  /media/datn/data/1st_DC_PRS_array_project/array_annotation/Process_annotation.R axiom ${array}.csv.gz ${array}_hg19.bed
  CrossMap.py region $chain_file ${array}_hg19.bed ${array}_hg38.bed
done

#hg38

for array in Axiom_JAPONICA Axiom_PMDA
do
  /media/datn/data/1st_DC_PRS_array_project/array_annotation/Process_annotation.R axiom ${array}_hg38.csv.gz ${array}_hg38.bed  
done



# ILLUMINA
cd /media/datn/data/1st_DC_PRS_array_project/array_annotation/annotaion_csv/illumina
chain_file="/media/datn/data/cross_map_lift_over/GRCh37_to_GRCh38.chain.gz"
# hg19

for array in chinese-genotyping-array-v1.0 infinium-core-v1.2 infinium-omni2.5.v1.5 infinium-omnizhonghua-v1.4 infinium-psycharray-v1.3 multi-ethnic-eur-eas-sas-v1.0 oncoarray-500k
do
  /media/datn/data/1st_DC_PRS_array_project/array_annotation/Process_annotation.R illumina ${array}.csv.gz ${array}_hg19.bed
  CrossMap.py region $chain_file ${array}_hg19.bed ${array}_hg38.bed  
done

# hg38

for array in cytosnp-850k-v1.2 global-screening-array-v.3 human-cytosnp-12-v2.1 infinium-exome-v1.0 infinium-omni5-v1.2 japanese-screening-array-v1.0 multi-ethnic-global-v1.0 infinium-global-diversity-array
do
  /media/datn/data/1st_DC_PRS_array_project/array_annotation/Process_annotation.R illumina ${array}_hg38.csv.gz ${array}_hg38.bed
done

## extract chromosome

cd /media/datn/data/1st_DC_PRS_array_project/array_annotation/annotation_all_hg38

for fi in *.bed 
do 
  /media/datn/data/1st_DC_PRS_array_project/array_annotation/Extract_chr_from_bed.R $fi
done