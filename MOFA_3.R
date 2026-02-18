setwd("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc")
scratch_lib <- "~/scratch/Rlibs"
dir.create(scratch_lib, showWarnings = FALSE)

# Prepend to lib paths
.libPaths(c(scratch_lib, .libPaths()))
print(.libPaths())  # confirm your scratch folder is first

BiocManager::install("Rsamtools", lib=scratch_lib)

#setting the library path
# Set library path
.libPaths(c("/cephfs/volumes/hpc_home/k25093549/002238a3-baec-473e-946f-6759b7d9f32d/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))

#NEW PATH, PROPER PATH
.libPaths(c("/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
.libPaths()



.libPaths()
library(MOFA2)
library(data.table)
library(tidyverse)
library(Seurat)
library(Signac)
library(BiocManager)

library(Rsamtools)



setwd("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data")




#creating a pseudobulk
sc_obj <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/20260204_merged_clustered_sc_data_clean.rds")
sc_obj

cat("=== SEURAT OBJECT INFO ===\n")
print(sc_obj)

cat("\n=== ASSAYS ===\n")
print(Assays(sc_obj))

cat("\n=== METADATA COLUMNS ===\n")
print(colnames(sc_obj@meta.data))

cat("\n=== FIRST FEW METADATA ROWS ===\n")
print(head(sc_obj@meta.data, 3))

cat("\n=== CELL TYPE DISTRIBUTION ===\n")
print(table(sc_obj@meta.data$celltype))

cat("\n=== INDIVIDUAL DISTRIBUTION ===\n")
print(table(sc_obj@meta.data$orig.ident))

#check distribution of samples

# The key columns are:
# - "sample" for individual ID (e.g., IGF136865)
# - "celltype" for cell type
# - "pd_status" for disease status

# Check sample distribution
cat("=== SAMPLES ===\n")
print(table(sc_obj@meta.data$sample))

cat("\n=== PD STATUS ===\n")
print(table(sc_obj@meta.data$pd_status))

cat("\n=== SAMPLES × PD STATUS ===\n")
print(table(sc_obj@meta.data$sample, sc_obj@meta.data$pd_status))

cat("\n=== SAMPLES × CELL TYPE ===\n")
print(table(sc_obj@meta.data$sample, sc_obj@meta.data$celltype))

library(dplyr)

# Add a simplified cell type column that matches bulk data
sc_obj$celltype_bulk <- sc_obj$celltype

#convert factor to charcter
sc_obj$celltype_bulk <- as.character(sc_obj$celltype)

# Group MSNs and Interneurons as "Neurons"
sc_obj$celltype_bulk[sc_obj$celltype_bulk %in% c("MSNs", 
                                                 "MSNs (DRD1+ enriched)",                                                "MSNs (DRD2+ enriched)", 
                                                 "Interneurons")] <- "Neurons"
# Group oligodendrocytes (including OPCs)
sc_obj$celltype_bulk[sc_obj$celltype_bulk %in% c("Oligodendrocytes (OPALIN+)", 
                                                 "Oligodendrocytes (PLEKHG1+)",
                                                 "Oligodendrocyte precursor cells (OPCs)")] <- "Oligodendrocytes"

# Keep Microglia as is
sc_obj$celltype_bulk[sc_obj$celltype_bulk == "Microglia"] <- "Microglia"



# Mark Astrocytes and Vascular for exclusion
sc_obj$celltype_bulk[sc_obj$celltype_bulk %in% c("Astrocytes", "Vascular")] <- "Exclude"

# Check the new groupings
cat("NEW CELL TYPE GROUPINGS\n")
print(table(sc_obj$celltype_bulk))

# Filter to keep only the 3 main cell types
sc_obj_filtered <- subset(sc_obj, subset = celltype_bulk %in% c("Neurons", "Oligodendrocytes", "Microglia"))

cat("\nFILTERED SAMPLES × CELL TYPES\n")
print(table(sc_obj_filtered$sample, sc_obj_filtered$celltype_bulk))

cat("\nTOTAL CELLS PER SAMPLE\n")
print(table(sc_obj_filtered$sample))

#Pseudobulk 

# Create pseudobulk for each cell type
cell_types <- c("Neurons", "Oligodendrocytes", "Microglia")

pseudobulk_matrices <- list()

for (ct in cell_types) {
  cat("Processing:", ct, "\n")
  
  # Subset to this cell type
  cells_subset <- subset(sc_obj_filtered, subset = celltype_bulk == ct)
  
  # Check how many cells per individual
  print(table(cells_subset$sample))
  
  # Aggregate counts by individual (sample)
  pseudobulk <- AggregateExpression(
    cells_subset,
    group.by = "sample",        # Aggregate by individual
    assays = "peaks",            # Use peaks assay
    slot = "counts",             # Use raw counts
    return.seurat = FALSE
  )
  
  # Extract the matrix
  pseudobulk_matrix <- pseudobulk$peaks
  
  # Store it
  pseudobulk_matrices[[ct]] <- pseudobulk_matrix
  
  cat("Dimensions:", dim(pseudobulk_matrix), "\n")
  cat("Sample names:", colnames(pseudobulk_matrix), "\n\n")
}

# Check what you have
names(pseudobulk_matrices)
lapply(pseudobulk_matrices, dim)

# Check the pseudobulk matrices
cat("PSEUDOBULK SUMMARY\n")
names(pseudobulk_matrices)

cat("\nDIMENSIONS\n")
lapply(pseudobulk_matrices, dim)

cat("\nSAMPLE NAMES\n")
lapply(pseudobulk_matrices, colnames)

cat("\nFIRST FEW PEAK NAMES\n")
head(rownames(pseudobulk_matrices$Neurons))

cat("\nTOTAL COUNTS\n")
lapply(pseudobulk_matrices, sum)

# List files in the bulk data directory
bulk_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/norm_peak_counts"

cat("FILES IN BULK DATA DIRECTORY:\n")
list.files(bulk_dir)

cat("\nFULL FILE PATHS:\n")
list.files(bulk_dir, full.names = TRUE) 

bulk_neurons <- read.csv(file.path(bulk_dir, "20260210_neuro_noage_norm_peak_counts_N93431.csv"), 
                         row.names = 1)
bulk_oligos <- read.csv(file.path(bulk_dir, "20260210_oligo_noage_norm_peak_counts_N93431.csv"), 
                        row.names = 1)
bulk_microglia <- read.csv(file.path(bulk_dir, "20260210_micro_noage_norm_peak_counts_N93431.csv"), 
                           row.names = 1)

cat("BULK DATA DIMENSIONS\n")
cat("Neurons:", dim(bulk_neurons), "\n")
cat("Oligodendrocytes:", dim(bulk_oligos), "\n")
cat("Microglia:", dim(bulk_microglia), "\n")

cat("\nBULK SAMPLE NAMES\n")
cat("Neurons:", colnames(bulk_neurons), "\n")
cat("Oligodendrocytes:", colnames(bulk_oligos), "\n")
cat("Microglia:", colnames(bulk_microglia), "\n")

# Check peak names
cat("\nFIRST FEW BULK PEAK NAMES\n")
head(rownames(bulk_neurons))

# List files in metadata directory
meta_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/metadata"

list.files(meta_dir, full.names = TRUE)

# Load bulk metadata
bulk_meta <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/metadata/sample_metadata_bulk.csv")

# Check the structure
head(bulk_meta)

# Check column names
colnames(bulk_meta)

# See all samples
print(bulk_meta)

#GOING TO LINK THESE DATASETS BY THE PID NUMBER

# Check unique person IDs in single-cell data
cat("SINGLE-CELL SAMPLES AND METADATA\n")
sc_metadata <- unique(sc_obj_filtered@meta.data[, c("sample", "pd_status", "age", "pmi")])
sc_metadata <- sc_metadata[order(sc_metadata$sample), ]
print(sc_metadata)

# Check unique PIDs in bulk
cat("\nBULK PERSON IDs\n")
unique(bulk_meta$pid)

# Count samples per person in bulk
cat("\nBULK SAMPLES PER PERSON\n")
table(bulk_meta$pid)

#trying to match them up using age and pmi

# Create mapping based on age and pmi
# Single-cell samples with their metadata
sc_meta_unique <- data.frame(
  sc_sample = c("IGF136865", "IGF136866", "IGF136867", "IGF136868", "IGF136869", 
                "IGF136870", "IGF136871", "IGF136872", "IGF136873", "IGF136874"),
  age = c(82, 80, 90, 66, 80, 70, 75, 65, 67, 94),
  pmi = c(20, 24, 22, 32, 10, 10, 46, 39, 45, 24),
  pd_status_sc = c("control", "pd", "control", "pd", "control", "pd", "control", "pd", "pd", "control")
)

# Get unique bulk person info
bulk_person <- unique(bulk_meta[, c("pid", "age", "pmi", "group")])

# Merge to match samples
sample_mapping <- merge(sc_meta_unique, bulk_person, by = c("age", "pmi"))

# Sort by sc_sample
sample_mapping <- sample_mapping[order(sample_mapping$sc_sample), ]

cat("SAMPLE MAPPING\n")
print(sample_mapping)

# Check which single-cell sample is missing from bulk
missing_sc <- setdiff(sc_meta_unique$sc_sample, sample_mapping$sc_sample)
cat("\nSINGLE-CELL SAMPLE MISSING FROM BULK:\n")
print(missing_sc)

#couldnt find IGF136868 in the bulk data
# Issue 1: IGF136868 (age 66, pmi 32)
cat("IGF136868 (SC sample - age 66, pmi 32)\n")
cat("Single-cell pd_status:", unique(sc_obj_filtered@meta.data[sc_obj_filtered@meta.data$sample == "IGF136868", "pd_status"]), "\n")
cat("\nBulk matches:\n")
print(bulk_meta[bulk_meta$age == 66 & bulk_meta$pmi == 32, c("sample", "pid", "age", "pmi", "group")])


#remove IGF136869 from pseudobulk matrices (this sample was removed from bulk when cleaned)
samples_to_keep <- setdiff(colnames(pseudobulk_matrices$Neurons), "IGF136869")
cat("Samples to keep:", samples_to_keep, "\n")

pseudobulk_matrices$Neurons <- pseudobulk_matrices$Neurons[, samples_to_keep]
pseudobulk_matrices$Oligodendrocytes <- pseudobulk_matrices$Oligodendrocytes[, samples_to_keep]
pseudobulk_matrices$Microglia <- pseudobulk_matrices$Microglia[, samples_to_keep]

cat("\nNEW PSEUDOBULK DIMENSIONS\n")
lapply(pseudobulk_matrices, dim)
cat("\nNEW SAMPLE NAMES\n")
colnames(pseudobulk_matrices$Neurons)

#results in 9 pseudobulk samples matching the bulk samples

#tackle 2 issues of weird data
#need to recode according to pid (absolute truth)

#Create sample mapping based on age/pmi and correct pd_status from PID
sample_mapping <- data.frame(
  sc_sample = c("IGF136865", "IGF136866", "IGF136867", "IGF136868", 
                "IGF136870", "IGF136871", "IGF136872", "IGF136873", "IGF136874"),
  bulk_pid = c("PDC126", "PD1222", "PDC167", "PD1231", 
               "PDC833", "PDC197", "PD726", "PD936", "PDC87"),
  age = c(82, 80, 90, 66, 70, 75, 65, 67, 94),
  pmi = c(20, 24, 22, 32, 10, 46, 39, 45, 24)
)

# Derive correct pd_status from PID naming
# PD#### = Parkinson's disease
# PDC#### = Control
sample_mapping$pd_status <- ifelse(grepl("^PD[0-9]", sample_mapping$bulk_pid), "pd", "control")

cat("CORRECTED SAMPLE MAPPING\n")
print(sample_mapping)

sample_mapping$bulk_neuron_id <- NA
sample_mapping$bulk_oligo_id <- NA
sample_mapping$bulk_micro_id <- NA


for(i in 1:nrow(sample_mapping)) {
  pid <- sample_mapping$bulk_pid[i]
  
  sample_mapping$bulk_neuron_id[i] <- bulk_meta$sample[bulk_meta$pid == pid & bulk_meta$celltype == "0:neurons"]
  sample_mapping$bulk_oligo_id[i] <- bulk_meta$sample[bulk_meta$pid == pid & bulk_meta$celltype == "1:oligodendrocytes"]
  sample_mapping$bulk_micro_id[i] <- bulk_meta$sample[bulk_meta$pid == pid & bulk_meta$celltype == "2:microglia"]
}

cat("\nFINAL MAPPING WITH BULK SAMPLE IDs\n")
print(sample_mapping)

cat("\nPD STATUS COUNTS\n")
table(sample_mapping$pd_status)



#a bulk sample is missing the microglia
# Check which PID is causing the problem
for(i in 1:nrow(sample_mapping)) {
  pid <- sample_mapping$bulk_pid[i]
  
  cat("Checking PID:", pid, "\n")
  
  # Check what samples exist for this PID
  cat("  Samples in bulk_meta:", bulk_meta$sample[bulk_meta$pid == pid], "\n")
  cat("  Cell types:", as.character(bulk_meta$celltype[bulk_meta$pid == pid]), "\n\n")
}


#PD833 should be a control (pasting error)
# Check PDC833 and PD833 in bulk metadata
cat("Checking PD833 and PDC833:\n")
print(bulk_meta[bulk_meta$pid %in% c("PD833", "PDC833"), c("sample", "pid", "celltype", "group")])

#recode PD833 TO PDC833
bulk_meta$pid[bulk_meta$pid == "PD833"] <- "PDC833"

# Verify the change
cat("CHECKING PDC833 samples after recoding:\n")
print(bulk_meta[bulk_meta$pid == "PDC833", c("sample", "pid", "celltype", "group")])

#VERIFY PD833 NO LONGER EXSISTS
cat("\nPD833 should be empty now:\n")
print(bulk_meta[bulk_meta$pid == "PD833", ])

# Check unique PIDs
cat("\nUNIQUE PIDs after recoding:\n")
print(unique(bulk_meta$pid))

#recode group
bulk_meta$group[bulk_meta$pid == "PDC833"] <- "0:control"
print(bulk_meta[bulk_meta$pid == "PDC833", c("sample", "pid", "celltype", "group")])

#check currrent order of psuedobulk
print(colnames(pseudobulk_matrices$Neurons))

cat("\nBULK SAMPLE IDs FROM MAPPING:\n")
print(sample_mapping[, c("sc_sample", "bulk_neuron_id", "bulk_oligo_id", "bulk_micro_id")])

#reorder the bulk to match pseudobulk
bulk_neurons_ordered <- bulk_neurons[, sample_mapping$bulk_neuron_id]
bulk_oligos_ordered <- bulk_oligos[, sample_mapping$bulk_oligo_id]
bulk_microglia_ordered <- bulk_microglia[, sample_mapping$bulk_micro_id]
#converting to charcters
sample_mapping$bulk_neuron_id <- as.character(sample_mapping$bulk_neuron_id)
sample_mapping$bulk_oligo_id <- as.character(sample_mapping$bulk_oligo_id)
sample_mapping$bulk_micro_id <- as.character(sample_mapping$bulk_micro_id)

colnames(bulk_neurons)

sample_mapping$bulk_neuron_id
#issue with mapping
# Rebuild mapping
sample_mapping <- data.frame(
  sc_sample = c("IGF136865", "IGF136866", "IGF136867", "IGF136868", 
                "IGF136870", "IGF136871", "IGF136872", "IGF136873", "IGF136874"),
  bulk_pid = c("PDC126", "PD1222", "PDC167", "PD1231", 
               "PDC833", "PDC197", "PD726", "PD936", "PDC87"),
  age = c(82, 80, 90, 66, 70, 75, 65, 67, 94),
  pmi = c(20, 24, 22, 32, 10, 46, 39, 45, 24)
)

# Get bulk IDs
for(i in 1:nrow(sample_mapping)) {
  pid <- sample_mapping$bulk_pid[i]
  sample_mapping$bulk_neuron_id[i] <- bulk_meta$sample[bulk_meta$pid == pid & bulk_meta$celltype == "0:neurons"]
  sample_mapping$bulk_oligo_id[i] <- bulk_meta$sample[bulk_meta$pid == pid & bulk_meta$celltype == "1:oligodendrocytes"]
  sample_mapping$bulk_micro_id[i] <- bulk_meta$sample[bulk_meta$pid == pid & bulk_meta$celltype == "2:microglia"]
}

print(sample_mapping)

#this shoudl work
# Reorder bulk to match pseudobulk

#convert bulk data to matrices from dat.fram
bulk_neurons <- as.matrix(bulk_neurons)
bulk_oligos <- as.matrix(bulk_oligos)
bulk_microglia <- as.matrix(bulk_microglia)


colnames(bulk_neurons)


# Reload bulk neurons to see ALL columns
bulk_neurons_full <- read.csv(file.path(bulk_dir, "20260210_neuro_noage_norm_peak_counts_N93431.csv"), row.names = 1)

colnames(bulk_neurons_full)

# See which samples are actually in bulk files
bulk_meta[bulk_meta$pid == "PDC126", ]

# Check if PDC126 neuron sample exists
bulk_meta[bulk_meta$sample == "IGF136815", ]


#need to remove this as well from pseudobulk, excluded from bulk anlalysis
#removing 1GF 136865 (PFC126) not in bulk from scratch 
samples_to_keep <- setdiff(colnames(pseudobulk_matrices$Neurons), "IGF136865")

pseudobulk_matrices$Neurons <- pseudobulk_matrices$Neurons[, samples_to_keep]
pseudobulk_matrices$Oligodendrocytes <- pseudobulk_matrices$Oligodendrocytes[, samples_to_keep]
pseudobulk_matrices$Microglia <- pseudobulk_matrices$Microglia[, samples_to_keep]

sample_mapping <- sample_mapping[sample_mapping$sc_sample != "IGF136865", ]
lapply(pseudobulk_matrices, dim)
print(sample_mapping)

#reodering should now work

# Convert to matrix
bulk_neurons <- as.matrix(bulk_neurons)
bulk_oligos <- as.matrix(bulk_oligos)
bulk_microglia <- as.matrix(bulk_microglia)

# Reorder
bulk_neurons_ordered <- bulk_neurons[, sample_mapping$bulk_neuron_id]
bulk_oligos_ordered <- bulk_oligos[, sample_mapping$bulk_oligo_id]
bulk_microglia_ordered <- bulk_microglia[, sample_mapping$bulk_micro_id]

# Rename columns to match pseudobulk
colnames(bulk_neurons_ordered) <- sample_mapping$sc_sample
colnames(bulk_oligos_ordered) <- sample_mapping$sc_sample
colnames(bulk_microglia_ordered) <- sample_mapping$sc_sample

# Check match
colnames(pseudobulk_matrices$Neurons)
colnames(bulk_neurons_ordered)


# 2)try and match peaks
# Get all unique peaks
all_peaks <- unique(c(
  rownames(bulk_neurons_ordered),
  rownames(bulk_oligos_ordered),
  rownames(bulk_microglia_ordered),
  rownames(pseudobulk_matrices$Neurons)
))

length(all_peaks)
# Check peak name formats
head(rownames(bulk_neurons_ordered))  # Bulk peaks
head(rownames(pseudobulk_matrices$Neurons))  # Pseudobulk peaks
#peak names dont match

# Check raw peak files
raw_neurons <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/raw_peak_counts/20260123_neuron_raw_peak_counts_N95726.csv", 
                        row.names = 1, nrows = 10)

head(rownames(raw_neurons))

save(pseudobulk_matrices, bulk_neurons_ordered, bulk_oligos_ordered, 
     bulk_microglia_ordered, sample_mapping, 
     file = "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/mofa_prep_data.RData")

#load saved objects from previous session
load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/mofa_prep_data.RData")


setwd("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/norm_peak_counts")

#check updated bulk data 
bulk_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/norm_peak_counts"

bulk_neurons <- read.csv(file.path(bulk_dir, "20260218_neuro_noage_norm_peak_counts_N93969.csv"), row.names = 1)
bulk_oligos <- read.csv(file.path(bulk_dir, "20260218_oligo_noage_raw_peak_counts_N99464.csv"), row.names = 1)
bulk_microglia <- read.csv(file.path(bulk_dir, "20260218_micro_noage_norm_peak_counts_N91620.csv"), row.names = 1)

# Check structure
dim(bulk_neurons)
colnames(bulk_neurons)[1:5]
head(rownames(bulk_neurons))

#load updated metadata
bulk_meta <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/metadata/sample_metadata_bulk.csv")

# Check if typos are fixed
unique(bulk_meta$pid)

#PDC833 is now PD833 in updated metdata

# Create peak IDs from coordinates
bulk_neurons$peak_id <- paste(bulk_neurons$chr, bulk_neurons$start, bulk_neurons$end, sep = "-")
bulk_oligos$peak_id <- paste(bulk_oligos$chr, bulk_oligos$start, bulk_oligos$end, sep = "-")
bulk_microglia$peak_id <- paste(bulk_microglia$chr, bulk_microglia$start, bulk_microglia$end, sep = "-")

#oligos doesnt have coordinate columns

colnames(bulk_oligos)[1:10]





