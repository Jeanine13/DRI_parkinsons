library (MOFA2)
.libPaths("/scratch/users/k25093549/R_packages")
library(DESeq2)

sc_data <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/20260204_merged_clustered_sc_data_clean.rds")
class(sc_data)

library(seurat)
library(HDF5Array)
module load R/4.5.0
version
install.packages("Seurat")
