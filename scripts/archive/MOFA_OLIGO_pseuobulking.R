library(Seurat)

# PSEUDOBULKING SCRIPT

#loading in new data

# Load the RDS object
sc_data <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/data_out_LA_dup_rm/081_analysis/001_deseq2/20260227_diff_analysis.rds")
str(sc_data, max.level = 2)

# See which cell types are available
names(sc_data$per_cluster)

#load in bulk data (VST) (this has been updated, was previosly not vst normalised)
#oligo_bulk_data <-read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/vst_norm_peak_counts/20260227_oligo_noage_vst_peak_counts_N99464.csv")
oligo_bulk_data <-read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/vst_norm_peak_counts/20260312_oligo_noage_vst_peak_counts_N99464.csv")

#load in bulk meta data
bulk_meta <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/metadata/sample_metadata_bulk.csv")

#access sc meta_data
sc_data$per_cluster$Oligodendrocytes.OPALIN.$colData
names(sc_data)

#drop sample missing from BULK (IGF136865 from sc)
#extract vst counts
opalin_vst  <- sc_data$per_cluster$Oligodendrocytes.OPALIN.$vst_counts
plekhg1_vst <- sc_data$per_cluster$Oligodendrocytes.PLECHG1.$vst_counts
opc_vst <- sc_data$per_cluster$Oligodendrocyte.precursor.cells..OPCs.$vst_counts

head(colnames(sc_data$per_cluster$Oligodendrocytes.PLECHG1.$vst_counts))


# drop IGF136865
opalin_vst  <- opalin_vst[, !grepl("IGF136865", colnames(opalin_vst))]
plekhg1_vst <- plekhg1_vst[, !grepl("IGF136865", colnames(plekhg1_vst))]
opc_vst <- opc_vst[, !grepl("IGF136865", colnames(opc_vst))]
                       colnames(opc_vst)
dim(opc_vst)

dim(oligo_bulk_data)
head(oligo_bulk_data)

#make count matrices for bulk data
bulk_counts <- oligo_bulk_data[, grepl("^IGF", colnames(oligo_bulk_data))]
rownames(bulk_counts) <- paste(oligo_bulk_data$seqnames, oligo_bulk_data$start, oligo_bulk_data$end, sep = "_")

dim(bulk_counts)

sample_mapping <- data.frame(
  pid       = c("PD1222", "PD1231", "PD726", "PD833", "PD936",
                "PDC167", "PDC184", "PDC197", "PDC87"),
  sc_sample = c("IGF136866", "IGF136868", "IGF136872", "IGF136870", "IGF136873",
                "IGF136867", "IGF136869", "IGF136871", "IGF136874"),
  pd_status = c("pd", "pd", "pd", "pd", "pd",
                "control", "control", "control", "control")
)

print(sample_mapping)

#add bulk_sample column to mapping to link sc and bulk data

oligo_bulk_meta <- bulk_meta[bulk_meta$celltype == "1:oligodendrocytes", ]
sample_mapping$bulk_sample <- oligo_bulk_meta$sample[match(sample_mapping$pid, oligo_bulk_meta$pid)]

print(sample_mapping)

colnames(opalin_vst)
colnames(bulk_counts)

#rename columns to PID

# strip suffix from sc columns
colnames(opalin_vst)  <- gsub("_Oligodendrocytes.*|_Oligodendrocyte.*", "", colnames(opalin_vst))
colnames(plekhg1_vst) <- gsub("_Oligodendrocytes.*|_Oligodendrocyte.*", "", colnames(plekhg1_vst))
colnames(opc_vst)     <- gsub("_Oligodendrocytes.*|_Oligodendrocyte.*", "", colnames(opc_vst))

# rename both to PID
colnames(bulk_counts) <- sample_mapping$pid[match(colnames(bulk_counts), sample_mapping$bulk_sample)]
colnames(opalin_vst)  <- sample_mapping$pid[match(colnames(opalin_vst),  sample_mapping$sc_sample)]
colnames(plekhg1_vst) <- sample_mapping$pid[match(colnames(plekhg1_vst), sample_mapping$sc_sample)]
colnames(opc_vst)     <- sample_mapping$pid[match(colnames(opc_vst),     sample_mapping$sc_sample)]

colnames(bulk_counts)
colnames(opalin_vst)

#order them to be consistent
bulk_counts <- bulk_counts[, sort(colnames(bulk_counts))]
opalin_vst  <- opalin_vst[, sort(colnames(opalin_vst))]
plekhg1_vst <- plekhg1_vst[, sort(colnames(plekhg1_vst))]
opc_vst     <- opc_vst[, sort(colnames(opc_vst))]

colnames(bulk_counts)
colnames(opalin_vst)

save(bulk_counts, oligo_bulk_data, bulk_meta, sample_mapping,
     opalin_vst, plekhg1_vst, opc_vst,
     file = "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/mofa_OLIGO_pseudobulked_20260316.RData")


data.frame(sample = colnames(bulk_counts),
           total_counts = colSums(bulk_counts))


save(bulk_counts, oligo_bulk_data, bulk_meta, sample_mapping,
     opalin_vst, plekhg1_vst, opc_vst,
     file = "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_pseudobulked_20260316.RData")







