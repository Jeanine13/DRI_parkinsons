# CHIPSEEKER ANNOTATION SCRIPT
# RUN THIS SCRIPT WHEN PEAK OBJECTS ARE UPDATED
# CREATES ANNOTATION OBJECTS FOR ALL FOUR DATASETS WITH PEAK KEYS AND GENE SYMBOLS
# USED LATER FOR MOFA WEIGHT INTERPRETATION

.libPaths(c(
  "/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5",
  .libPaths()
))

# DEFINE OUTPUT DIRECTORY
mofa_out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/new_sc_mofa_clean_14042026"
dir.create(mofa_out_dir, recursive = TRUE, showWarnings = FALSE)

# LOAD GRANGES AND VST DATA
load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_granges_vst.RData")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
library(dplyr)
library(GenomeInfoDb)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# FIX CHROMOSOME NAMING TO MATCH TXDB
seqlevelsStyle(bulk_gr)    <- "UCSC"
seqlevelsStyle(opalin_gr)  <- "UCSC"
seqlevelsStyle(plekhg1_gr) <- "UCSC"
seqlevelsStyle(opc_gr)     <- "UCSC"

# ANNOTATE ALL PEAKS WITH GENIC LOCATIONS AND NEAREST GENE (ENTREZID)
bulk_anno_df <- as.data.frame(
  annotatePeak(
    bulk_gr,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )
)

opalin_anno_cp <- as.data.frame(
  annotatePeak(
    opalin_gr,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )
)

plekhg1_anno_cp <- as.data.frame(
  annotatePeak(
    plekhg1_gr,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )
)

opc_anno_cp <- as.data.frame(
  annotatePeak(
    opc_gr,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )
)

# ADD HUMAN-READABLE GENE SYMBOLS AND GENE NAMES USING ENSDB
add_gene_symbols <- function(anno_df) {
  annotations_edb <- AnnotationDbi::select(
    EnsDb.Hsapiens.v86,
    keys = as.character(anno_df$geneId),
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
  )
  
  annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)
  annotations_edb <- annotations_edb[!duplicated(annotations_edb$ENTREZID), ]
  
  anno_df %>%
    left_join(annotations_edb, by = c("geneId" = "ENTREZID"))
}

bulk_anno_df    <- add_gene_symbols(bulk_anno_df)
opalin_anno_cp  <- add_gene_symbols(opalin_anno_cp)
plekhg1_anno_cp <- add_gene_symbols(plekhg1_anno_cp)
opc_anno_cp     <- add_gene_symbols(opc_anno_cp)

# CREATE PEAK KEYS FOR MERGING WITH MOFA WEIGHTS
# BULK PEAKS ARE IN chr1_start_end FORMAT
# SC PEAKS ARE IN 1-start-end FORMAT
bulk_anno_df$peak_key <- paste(
  bulk_anno_df$seqnames,
  bulk_anno_df$start,
  bulk_anno_df$end,
  sep = "_"
)

opalin_anno_cp$peak_key <- paste(
  gsub("chr", "", opalin_anno_cp$seqnames),
  opalin_anno_cp$start,
  opalin_anno_cp$end,
  sep = "-"
)

plekhg1_anno_cp$peak_key <- paste(
  gsub("chr", "", plekhg1_anno_cp$seqnames),
  plekhg1_anno_cp$start,
  plekhg1_anno_cp$end,
  sep = "-"
)

opc_anno_cp$peak_key <- paste(
  gsub("chr", "", opc_anno_cp$seqnames),
  opc_anno_cp$start,
  opc_anno_cp$end,
  sep = "-"
)

# VERIFY ANNOTATION OBJECTS
cat("Bulk annotation peaks:", nrow(bulk_anno_df), "\n")
cat("Opalin annotation peaks:", nrow(opalin_anno_cp), "\n")
cat("Plekhg1 annotation peaks:", nrow(plekhg1_anno_cp), "\n")
cat("OPC annotation peaks:", nrow(opc_anno_cp), "\n")

# OPTIONAL PREVIEW OF KEY COLUMNS
cat("\nPreview of Opalin annotation:\n")
print(head(opalin_anno_cp[, c("seqnames", "start", "end", "annotation", "SYMBOL", "GENENAME", "peak_key")]))

# SAVE FOR DOWNSTREAM MOFA USE
save(
  bulk_anno_df,
  opalin_anno_cp,
  plekhg1_anno_cp,
  opc_anno_cp,
  file = file.path(mofa_out_dir, "chipseeker_annotations_all_datasets.RData")
)

cat("\nChIPseeker annotation complete and saved.\n")
cat("Saved file:", file.path(mofa_out_dir, "chipseeker_annotations_all_datasets.RData"), "\n")

