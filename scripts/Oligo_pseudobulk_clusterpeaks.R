.libPaths(c("/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))

library(Seurat)
library(Matrix)
# load cluster-specific SC data
opalin_sc_data  <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/cluster_1_seurat_oligo_opalin.rds")
plekhg1_sc_data <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/cluster_5_norm_oligo_plekhg1.rds")
opc_sc_data     <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/cluster_7_norm_opc.rds")

# load bulk data
oligo_bulk_data <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/vst_norm_peak_counts/20260312_oligo_noage_vst_peak_counts_N99464.csv")

# load bulk metadata
bulk_meta <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/metadata/sample_metadata_bulk.csv")

sample_mapping <- data.frame(
  pid       = c("PD1222", "PD1231", "PD726", "PD833", "PD936",
                "PDC167", "PDC184", "PDC197", "PDC87"),
  sc_sample = c("IGF136866", "IGF136868", "IGF136872", "IGF136870", "IGF136873",
                "IGF136867", "IGF136869", "IGF136871", "IGF136874"),
  pd_status = c("pd", "pd", "pd", "pd", "pd",
                "control", "control", "control", "control")
)

oligo_bulk_meta <- bulk_meta[bulk_meta$celltype == "1:oligodendrocytes", ]
sample_mapping$bulk_sample <- oligo_bulk_meta$sample[match(sample_mapping$pid, oligo_bulk_meta$pid)]

print(sample_mapping)


bulk_counts <- oligo_bulk_data[, grepl("^IGF", colnames(oligo_bulk_data))]
rownames(bulk_counts) <- paste(oligo_bulk_data$seqnames, oligo_bulk_data$start, oligo_bulk_data$end, sep = "_")
colnames(bulk_counts) <- sample_mapping$pid[match(colnames(bulk_counts), sample_mapping$bulk_sample)]
bulk_counts <- bulk_counts[, sort(colnames(bulk_counts))]

# check
dim(bulk_counts)
colSums(bulk_counts)

# extract normalised data layer from each subtype
opalin_norm  <- GetAssayData(opalin_sc_data,  assay = "peaks_cluster", layer = "data")
plekhg1_norm <- GetAssayData(plekhg1_sc_data, assay = "peaks_cluster", layer = "data")
opc_norm     <- GetAssayData(opc_sc_data,     assay = "peaks_cluster", layer = "data")

# check cell counts per donor
cat("Opalin+ cells per donor:\n")
print(table(gsub("_.*", "", colnames(opalin_norm))))
cat("Plekhg1+ cells per donor:\n")
print(table(gsub("_.*", "", colnames(plekhg1_norm))))
cat("OPC cells per donor:\n")
print(table(gsub("_.*", "", colnames(opc_norm))))

# pseudobulk by averaging across cells per donor
pseudobulk_mean <- function(norm_mat, sample_mapping) {
  cell_samples <- gsub("_.*", "", colnames(norm_mat))
  pseudo <- sapply(unique(cell_samples), function(s) {
    cells <- cell_samples == s
    rowMeans(norm_mat[, cells, drop = FALSE])
  })
  # drop IGF136865
  pseudo <- pseudo[, colnames(pseudo) != "IGF136865"]
  # rename to PID
  colnames(pseudo) <- sample_mapping$pid[match(colnames(pseudo), sample_mapping$sc_sample)]
  # sort columns
  pseudo[, sort(colnames(pseudo))]
}

opalin_pseudo_norm  <- pseudobulk_mean(opalin_norm,  sample_mapping)
plekhg1_pseudo_norm <- pseudobulk_mean(plekhg1_norm, sample_mapping)
opc_pseudo_norm     <- pseudobulk_mean(opc_norm,     sample_mapping)

# check dimensions and column sums
cat("Opalin+ dimensions:", dim(opalin_pseudo_norm), "\n")
cat("Plekhg1+ dimensions:", dim(plekhg1_pseudo_norm), "\n")
cat("OPC dimensions:", dim(opc_pseudo_norm), "\n")

colSums(opalin_pseudo_norm)
colSums(plekhg1_pseudo_norm)
colSums(opc_pseudo_norm)

# check column names match
colnames(bulk_counts)
colnames(opalin_pseudo_norm)
colnames(plekhg1_pseudo_norm)
colnames(opc_pseudo_norm)

save(bulk_counts, oligo_bulk_data, bulk_meta, sample_mapping,
     opalin_pseudo_norm, plekhg1_pseudo_norm, opc_pseudo_norm,
     file = "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_pseudobulked.RData")
