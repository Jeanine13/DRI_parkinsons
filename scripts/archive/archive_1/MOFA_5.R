#NEW PATH, PROPER PATH
.libPaths(c("/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
.libPaths()

library(MOFA2)
library(data.table)

library(Seurat)
library(Signac)
library(BiocManager)

library(Rsamtools)
#BULK DATA--

#load in bulk data
bulk_oligo <-read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/norm_peak_counts/20260218_oligo_noage_norm_peak_counts_N99464.csv")
bulk_meta <-read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/metadata/sample_metadata_bulk.csv")

# fix typo in bulk metadata
bulk_meta$pid[bulk_meta$pid == "PD833"]     <- "PDC833"
bulk_meta$group[bulk_meta$pid == "PDC833"]  <- "0:control"

#produce a clean bulk matrix
#define annotations columns vs count columns (removing annotations)
annot_cols <- c("X", "chr", "start", "end", "seqnames", "width", "strand", 
                "annotation", "geneChr", "geneStart", "geneEnd", "geneLength", 
                "geneStrand", "geneId", "transcriptId", "distanceToTSS", 
                "ENSEMBL", "SYMBOL", "GENENAME")
# extract just the count matrix
bulk_counts <- bulk_oligo[, !colnames(bulk_oligo) %in% annot_cols]

# extract annotation separately
bulk_annot <- bulk_oligo[, annot_cols]


#createa peak ID from corrindates to use as row names
rownames(bulk_counts) <- paste(bulk_oligo$seqnames, bulk_oligo$start, bulk_oligo$end, sep = "_")
rownames(bulk_annot) <- rownames(bulk_counts)

#SC DATA
library(Signac)
library(Seurat)
#read in sc data
sc_data <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/20260214_merged_sub_clustered_ga_nox.rds")

head(sc_data@meta.data)
table(sc_data$celltype)

# get counts matrix and metadata
counts <- GetAssayData(sc_data, assay = "peaks", layer = "counts")
meta <- sc_data@meta.data


# subset to oligo subtypes only
oligo_types <- c("Oligodendrocytes OPALIN+", "Oligodendrocytes PLECHG1+", 
                 "Oligodendrocyte precursor cells (OPCs)")
meta_oligo <- meta[meta$celltype %in% oligo_types, ]

# pseudobulk, sum counts per sample x celltype
meta_oligo$group_id <- paste(meta_oligo$sample, meta_oligo$celltype, sep = "_")
group_levels <- unique(meta_oligo$group_id)

pb_counts <- sapply(group_levels, function(g) {
  Matrix::rowSums(counts[, meta_oligo$group_id == g, drop = FALSE])
})


dim(pb_counts)

head(colnames(pb_counts))

#split into three matrices, 1 per subtype

opalin_counts  <- pb_counts[, grepl("OPALIN",  colnames(pb_counts))]
plekhg1_counts <- pb_counts[, grepl("PLECHG1", colnames(pb_counts))]
opc_counts     <- pb_counts[, grepl("OPCs",    colnames(pb_counts))]
dim(opalin_counts)
dim(plekhg1_counts)
dim(opc_counts)



#SAMPLE MATCHING BETWEEN SC AND BULK
#match sc to bulk

# get unique sample metadata from sc
sc_meta_unique <- unique(meta_oligo[, c("sample", "age", "pmi", "pd_status")])

# get unique person info from bulk
bulk_person <- unique(bulk_meta[, c("pid", "age", "pmi", "group")])

#view data, possible missing samples in sc
meta_samples <- bulk_meta$sample[bulk_meta$celltype == "1:oligodendrocytes"]
# samples in bulk counts
count_samples <- colnames(bulk_counts)
# which is in meta but not in counts
setdiff(meta_samples, count_samples)
 #MDropped sampled is IGF136816

sample_mapping <- merge(sc_meta_unique, bulk_person, by = c("age", "pmi"))

#drop it from sc
sample_mapping <- sample_mapping[sample_mapping$sample != "IGF136865", ]

samples_to_drop <- c("IGF136865", "IGF136869")


#dropping from oligo subtypes

opalin_counts  <- opalin_counts[, !grepl(paste(samples_to_drop, collapse="|"), colnames(opalin_counts))]
plekhg1_counts <- plekhg1_counts[, !grepl(paste(samples_to_drop, collapse="|"), colnames(plekhg1_counts))]
opc_counts     <- opc_counts[, !grepl(paste(samples_to_drop, collapse="|"), colnames(opc_counts))]

dim(opalin_counts)
dim(plekhg1_counts)
dim(opc_counts)


meta[meta$sample == "IGF136869", c("sample", "age", "pmi")] |> unique()
bulk_meta[bulk_meta$age == 80 & bulk_meta$pmi == 10, ]


#cant find IGF136869 in the bulk data but is in the sc data
#confim it in the single cell
unique(meta[meta$sample == "IGF136869", c("sample", "age", "pmi", "pd_status")])

#proceed without (dropped form sc data), have 8 samples
library(DESeq2)

#--------------Troubleshooting----------
#fixing issues with rownames (build col_data with opalin)
rownames(colData_sc) <- gsub("_Oligodendrocytes.*|_Oligodendrocyte.*", "", colnames(opalin_counts))

print(colData_sc)
colnames(opalin_counts)  <- gsub("_Oligodendrocytes.*|_Oligodendrocyte.*", "", colnames(opalin_counts))
colnames(plekhg1_counts) <- gsub("_Oligodendrocytes.*|_Oligodendrocyte.*", "", colnames(plekhg1_counts))
colnames(opc_counts)     <- gsub("_Oligodendrocytes.*|_Oligodendrocyte.*", "", colnames(opc_counts))
colnames(opalin_counts)

colnames(opalin_counts)
colnames(plekhg1_counts)
colnames(opc_counts)
#------------------------------
#DESeq normalisation 

#run DESeq2 on pseudobulk matrices

sc_samples <- gsub("_Oligodendrocytes.*|_Oligodendrocyte.*", "", colnames(opalin_counts))

colData_sc <- sample_mapping[match(sc_samples, sample_mapping$sample), c("sample", "pd_status")]
rownames(colData_sc) <- colnames(opalin_counts)
print(colData_sc)

#fix peaks names for figure so can compare
head(rownames(opalin_counts))
#reformat bulk peak names
rownames(bulk_counts) <- gsub("chr", "", rownames(bulk_counts))
rownames(bulk_counts) <- gsub("_", "-", rownames(bulk_counts))
head(rownames(bulk_counts))




