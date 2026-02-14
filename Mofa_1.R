setwd("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data")

sc_data <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/20260204_merged_clustered_sc_data_clean.rds")
ls()

class(sc_data)




# Install dependencies first
install.packages("BiocManager")
library(BiocManager)
BiocManager::install(c("multtest", "limma", "DESeq2"))

# Then try Seurat again
install.packages("Seurat", dependencies = TRUE)

colnames(sc_data@meta.data)
table(sc_data@meta.data$celltype)

#subsetting microglia
# First install the BiocManager dependency
BiocManager::install("Rsamtools")

# Install the missing dependencies first
install.packages(c("httr", "fitdistrplus", "MASS", "sctransform", "plotly"))

# Then try Seurat again
install.packages("Seurat")

#SEURAT not installing, using other method

#Get cell indcies
microglia_idx <- which(sc_data@meta.data$celltype == "Microglia")
length(microglia_idx)
names(sc_data@assays)
assay_name <- names(sc_data@assays)[1]
print(assay_name)

all_counts <- sc_data@assays[[assay_name]]@counts

# Check dimensions
dim(all_counts)  # Should show (number of peaks x total cells)
microglia_counts <- all_counts[, microglia_idx]
dim(microglia_counts)  # Should be (peaks x 4156)

#extract metadata from microglia
microglia_meta <- sc_data@meta.data[microglia_idx, ]

#verify extraction
nrow(microglia_meta)  # Should be 4156
table(microglia_meta$celltype)  # Should show only "Microglia"

#confirming how the data looks
head(rownames(microglia_counts))  # Peak names
head(colnames(microglia_counts)) # Cell barcodes

.libPaths("/scratch/users/k25093549/R_packages")

#loading bulk data in 
bulk_microglia <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/norm_peak_counts/20260210_micro_noage_norm_peak_counts_N93431.csv")
# Check the structure
head(bulk_microglia)
dim(bulk_microglia)
colnames(bulk_microglia)

#checking if first column is peak names 
head(bulk_microglia[,1])

head(rownames(bulk_microglia))
# Also show the first few columns to understand structure
bulk_microglia[1:5, 1:5]

# And column names
colnames(bulk_microglia)

#load in meta data (will extract only for microglia)
peak_meta <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/metadata/sample_metadata_bulk.csv")

head(peak_meta)
dim(peak_meta)
colnames(peak_meta)
table(microglia_meta$sample)
















