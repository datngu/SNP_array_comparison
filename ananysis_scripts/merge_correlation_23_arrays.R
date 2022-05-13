


#chr="22"
#out_fn = "test.txt.gz"

merge_23_array_by_chr <- function(pop, chr, array_list){
	require(data.table)
	#if(is.null(out_fn)) 
	out_fn = paste0("chr", chr, "_", pop, "_merged_23_array_correlation.txt.gz")
	array = array_list[1]
	path = paste0(array, "/chr", chr, "_", pop, "_correlation.txt.gz")
	df = fread(path)
	df = df[,c(1:4)]
	
	for(array in array_list){
		path = paste0(array, "/chr", chr, "_", pop, "_correlation.txt.gz")
		tem = fread(path)
		tem[,array] = tem$cor
		pick_col =  c("ID", array)
		tem = tem[, ..pick_col]
		df = merge(df, tem, by = "ID")
	}
	fwrite(df, file = out_fn, sep = "\t")
}


setwd("/media/datn/data/1st_DC_PRS_array_project/Imputation_correlation_VNP_1008")
chr_list = as.character(c(1:22))
array_list = c("Axiom_GW_ASI", "Axiom_GW_CHB", "Axiom_GW_EUR", "Axiom_GW_PanAFR", "Axiom_JAPONICA", "infinium-omnizhonghua-v1.4", "infinium-psycharray-v1.3", "japanese-screening-array-v1.0", "multi-ethnic-eur-eas-sas-v1.0", "multi-ethnic-global-v1.0", "oncoarray-500k", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB_WCSG", "chinese-genotyping-array-v1.0", "cytosnp-850k-v1.2", "global-screening-array-v.3", "human-cytosnp-12-v2.1", "infinium-core-v1.2", "infinium-global-diversity-array-v1.0", "infinium-omni2.5.v1.5", "infinium-omni5-v1.2", "GenomeWideSNP_6.0")


for(chr in chr_list){
	merge_23_array_by_chr(pop = pop, chr = chr, array_list = array_list)
}



setwd("/media/datn/data/1st_DC_PRS_array_project/igsr_10f_correlation_res")
chr_list = as.character(c(1:22))
array_list = c("Axiom_GW_ASI", "Axiom_GW_CHB", "Axiom_GW_EUR", "Axiom_GW_PanAFR", "Axiom_JAPONICA", "infinium-omnizhonghua-v1.4", "infinium-psycharray-v1.3", "japanese-screening-array-v1.0", "multi-ethnic-eur-eas-sas-v1.0", "multi-ethnic-global-v1.0", "oncoarray-500k", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB_WCSG", "chinese-genotyping-array-v1.0", "cytosnp-850k-v1.2", "global-screening-array-v.3", "human-cytosnp-12-v2.1", "infinium-core-v1.2", "infinium-global-diversity-array-v1.0", "infinium-omni2.5.v1.5", "infinium-omni5-v1.2", "GenomeWideSNP_6.0")

for(chr in chr_list){
	merge_23_array_by_chr(pop = "EAS", chr = chr, array_list = array_list)
	merge_23_array_by_chr(pop = "EUR", chr = chr, array_list = array_list)
	merge_23_array_by_chr(pop = "SAS", chr = chr, array_list = array_list)
	merge_23_array_by_chr(pop = "AMR", chr = chr, array_list = array_list)
	merge_23_array_by_chr(pop = "AFR", chr = chr, array_list = array_list)
}