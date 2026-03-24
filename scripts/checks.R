#checking if used the correct bulk counts

verify_oligo_bulk_vst <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/vst_norm_peak_counts/20260227_oligo_noage_vst_peak_counts_N99464.csv")
#im using the correct vst 

#load the norm count 
verify_oligo_bulk_countnorm <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/count_norm_peak_counts/20260227_oligo_noage_norm_peak_counts_N99464.csv")

new_bulk <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/vst_norm_peak_counts/20260312_oligo_noage_vst_peak_counts_N99464.csv")

dim(new_bulk)
head(colnames(new_bulk))

# extract count matrix
bulk_counts_new <- new_bulk[, grepl("^IGF", colnames(new_bulk))]
rownames(bulk_counts_new) <- paste(new_bulk$seqnames, new_bulk$start, new_bulk$end, sep = "_")

# check dimensions and column names
dim(bulk_counts_new)
colnames(bulk_counts_new)

colSums(bulk_counts_new)
colSums(bulk_counts)




#check dimensions of new pseudobuljed
dim(bulk_counts)
colnames(bulk_counts)
colSums(bulk_counts)
head(bulk_counts[, 1:3])


load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_pseudobulked_20260316.RData")

colSums(bulk_counts)
dim(bulk_counts)


head(bulk_counts)
head(oligo_bulk_data[, grepl("^IGF", colnames(oligo_bulk_data))])




