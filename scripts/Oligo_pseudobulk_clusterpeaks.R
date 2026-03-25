.libPaths(c("/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))

#PSEUDOLBULKING NEW SC DATA AND DESEQ2 NORMALISATION


library(Seurat)
library(Matrix)
library(DESeq2)

# load cluster-specific SC data
opalin_sc_data  <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/cluster_1_seurat_oligo_opalin.rds")
plekhg1_sc_data <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/cluster_5_norm_oligo_plekhg1.rds")
opc_sc_data     <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/cluster_7_norm_opc.rds")

# load bulk data
oligo_bulk_data <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/vst_norm_peak_counts/20260312_oligo_noage_vst_peak_counts_N99464.csv")

# load bulk metadata
bulk_meta <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/metadata/sample_metadata_bulk.csv")

#sample mapping 
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

#extract bulk count matrix 

bulk_counts <- oligo_bulk_data[, grepl("^IGF", colnames(oligo_bulk_data))]
rownames(bulk_counts) <- paste(oligo_bulk_data$seqnames, oligo_bulk_data$start, oligo_bulk_data$end, sep = "_")
colnames(bulk_counts) <- sample_mapping$pid[match(colnames(bulk_counts), sample_mapping$bulk_sample)]
bulk_counts <- bulk_counts[, sort(colnames(bulk_counts))]

# check
dim(bulk_counts)
colSums(bulk_counts)

# extract RAW COUNTS FROM COUNTS LAYER OF SC data
opalin_raw  <- GetAssayData(opalin_sc_data,  assay = "peaks_cluster", layer = "counts")
plekhg1_raw <- GetAssayData(plekhg1_sc_data, assay = "peaks_cluster", layer = "counts")
opc_raw     <- GetAssayData(opc_sc_data,     assay = "peaks_cluster", layer = "counts")

#pseudobulk by summing counts per sample
pseudobulk_sum <- function(count_mat, sample_mapping) {
  cell_samples <- gsub("_.*", "", colnames(count_mat))
  pseudo <- sapply(unique(cell_samples), function(s) {
    cells <- cell_samples == s
    rowSums(count_mat[, cells, drop = FALSE])
  })
  pseudo <- pseudo[, colnames(pseudo) != "IGF136865"]
  colnames(pseudo) <- sample_mapping$pid[match(colnames(pseudo), sample_mapping$sc_sample)]
  pseudo[, sort(colnames(pseudo))]
}

opalin_pseudo  <- pseudobulk_sum(opalin_raw,  sample_mapping)
plekhg1_pseudo <- pseudobulk_sum(plekhg1_raw, sample_mapping)
opc_pseudo     <- pseudobulk_sum(opc_raw,     sample_mapping)



# check cell counts per sample

cat("Opalin+ dimensions:", dim(opalin_pseudo), "\n")
cat("Plekhg1+ dimensions:", dim(plekhg1_pseudo), "\n")
cat("OPC dimensions:", dim(opc_pseudo), "\n")

colSums(opalin_pseudo)
colSums(plekhg1_pseudo)
colSums(opc_pseudo)



#DESeq2 VST normalisation

vst_normalise <- function(pseudo_mat) {
  # round to integers for DESeq2
  pseudo_mat <- round(as.matrix(pseudo_mat))
  
  # create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = pseudo_mat,
    colData   = data.frame(
      condition = rep("sample", ncol(pseudo_mat)),
      row.names = colnames(pseudo_mat)
    ),
    design = ~ 1
  )
  
  # estimate size factors first
  dds <- estimateSizeFactors(dds)
  
  # VST normalisation
  vst_mat <- assay(vst(dds, blind = TRUE))
  
  return(vst_mat)
}

#filter zero variance peaks after vst
filter_zero_var <- function(vst_mat) {
  keep <- apply(vst_mat, 1, var) > 0
  vst_mat[keep, , drop = FALSE]
}

opalin_vst  <- filter_zero_var(vst_normalise(opalin_pseudo))
plekhg1_vst <- filter_zero_var(vst_normalise(plekhg1_pseudo))
opc_vst     <- filter_zero_var(vst_normalise(opc_pseudo))


dim(opalin_vst)
summary(as.vector(opalin_vst))
head(opalin_vst[, 1:3])


#check peaks removed
tmp_opalin <- vst_normalise(opalin_pseudo)
sum(apply(tmp_opalin, 1, var) == 0)

opalin_vst <- filter_zero_var(tmp_opalin)

colnames(opalin_vst)
colnames(plekhg1_vst)
colnames(opc_vst)

colSums(opalin_vst)
colSums(plekhg1_vst)
colSums(opc_vst)

save(bulk_counts, oligo_bulk_data, bulk_meta, sample_mapping,
     opalin_vst, plekhg1_vst, opc_vst,
     opalin_pseudo, plekhg1_pseudo, opc_pseudo,
     file = "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_pseudobulked_vst.RData")



