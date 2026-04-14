# MOFA ANALYSIS SCRIPT
# IMPORTANT: RUN CHIPSEEKER ANNOTATION SCRIPT BEFORE THIS
# CHIPSEEKER SCRIPT CREATES: opalin_anno_cp, plekhg1_anno_cp, opc_anno_cp, bulk_anno_df
# ALL WITH peak_key AND GENENAME.x COLUMNS ALREADY ADDED (dont need to worry about this)

.libPaths(c("/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))

library(MOFA2)
library(clusterProfiler)
library(rGREAT)
library(ggplot2)
library(org.Hs.eg.db)

# DO NOT RELOAD GRANGES RDATA HERE (IT WOULD OVERWRITE ANNOTATIONS FROM CHIPSEEKER SCRIPT)
# INSTEAD VERIFY THE REQUIRED OBJECTS ARE IN MEMORY
stopifnot(exists("opalin_anno_cp"))
stopifnot(exists("plekhg1_anno_cp"))
stopifnot(exists("opc_anno_cp"))
stopifnot(exists("bulk_anno_df"))
stopifnot(exists("opalin_vst"))
stopifnot(exists("plekhg1_vst"))
stopifnot(exists("opc_vst"))
stopifnot(exists("bulk_counts"))

cat("All required objects found in memory\n")

mofa_out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/new_sc_mofa"
dir.create(mofa_out_dir, recursive = TRUE, showWarnings = FALSE)

# SAMPLE METADATA FOR PD ASSOCIATION TESTING
sample_meta <- data.frame(
  sample    = c("PD1222", "PD1231", "PD726", "PD833", "PD936",
                "PDC167", "PDC184", "PDC197", "PDC87"),
  condition = c("PD", "PD", "PD", "PD", "PD",
                "Control", "Control", "Control", "Control")
)

# SELECT TOP 5000 HVFs PER VIEW
bulk_vars    <- apply(bulk_counts, 1, var)
bulk_top_idx <- head(order(bulk_vars, decreasing = TRUE), 5000)
bulk_hvf     <- bulk_counts[bulk_top_idx, ]

opalin_vars    <- apply(opalin_vst, 1, var)
opalin_top_idx <- head(order(opalin_vars, decreasing = TRUE), 5000)
opalin_hvf     <- opalin_vst[opalin_top_idx, ]

plekhg1_vars    <- apply(plekhg1_vst, 1, var)
plekhg1_top_idx <- head(order(plekhg1_vars, decreasing = TRUE), 5000)
plekhg1_hvf     <- plekhg1_vst[plekhg1_top_idx, ]

opc_vars    <- apply(opc_vst, 1, var)
opc_top_idx <- head(order(opc_vars, decreasing = TRUE), 5000)
opc_hvf     <- opc_vst[opc_top_idx, ]

# VERIFY HVF PEAKS MATCH ANNOTATION PEAK KEYS
cat("Opalin HVF in annotation:", sum(rownames(opalin_hvf) %in% opalin_anno_cp$peak_key), "/ 5000\n")
cat("Plekhg1 HVF in annotation:", sum(rownames(plekhg1_hvf) %in% plekhg1_anno_cp$peak_key), "/ 5000\n")
cat("OPC HVF in annotation:", sum(rownames(opc_hvf) %in% opc_anno_cp$peak_key), "/ 5000\n")
cat("Bulk HVF in annotation:", sum(rownames(bulk_hvf) %in% bulk_anno_df$peak_key), "/ 5000\n")

# LOAD BACKGROUND CONSENSUS PEAK SETS FOR GO ENRICHMENT UNIVERSE
sc_consensus   <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/single_cell_consensus_peaks.rds")
bulk_consensus <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/union_peaks/bulk_union_peaks.rds")

seqlevelsStyle(sc_consensus)   <- "UCSC"
seqlevelsStyle(bulk_consensus) <- "UCSC"

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

sc_consensus_anno   <- as.data.frame(annotatePeak(sc_consensus,   tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db"))
bulk_consensus_anno <- as.data.frame(annotatePeak(bulk_consensus, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db"))

sc_background_entrez   <- unique(sc_consensus_anno$geneId[!is.na(sc_consensus_anno$geneId)])
bulk_background_entrez <- unique(bulk_consensus_anno$geneId[!is.na(bulk_consensus_anno$geneId)])

cat("SC background genes:", length(sc_background_entrez), "\n")
cat("Bulk background genes:", length(bulk_background_entrez), "\n")

# FUNCTION TO CONVERT BULK PEAK NAMES TO GRANGES
bulk_peaks_to_gr <- function(peaks) {
  GRanges(
    seqnames = gsub("_.*", "", peaks),
    ranges   = IRanges(
      start = as.numeric(gsub(".*_(\\d+)_\\d+$", "\\1", peaks)),
      end   = as.numeric(gsub(".*_(\\d+)$", "\\1", peaks))
    )
  )
}

# FUNCTION TO CONVERT SC PEAK NAMES TO GRANGES
sc_peaks_to_gr <- function(peaks) {
  GRanges(
    seqnames = paste0("chr", sapply(strsplit(as.character(peaks), "-"), `[`, 1)),
    ranges   = IRanges(
      start = as.numeric(sapply(strsplit(as.character(peaks), "-"), `[`, 2)),
      end   = as.numeric(sapply(strsplit(as.character(peaks), "-"), `[`, 3))
    )
  )
}

# GO ENRICHMENT FUNCTION WITH SAFETY CHECKS AND DIAGNOSTICS
run_go_factor <- function(weights_df, anno_df, peak_key_col, geneid_col, label,
                          background_entrez,
                          n_top = 500, mode = c("abs", "positive", "negative")) {
  
  mode <- match.arg(mode)
  
  if (mode == "abs") {
    top_peaks <- weights_df[order(abs(weights_df$value), decreasing = TRUE), ]
  } else if (mode == "positive") {
    top_peaks <- weights_df[weights_df$value > 0, ]
    top_peaks <- top_peaks[order(top_peaks$value, decreasing = TRUE), ]
  } else if (mode == "negative") {
    top_peaks <- weights_df[weights_df$value < 0, ]
    top_peaks <- top_peaks[order(top_peaks$value, decreasing = FALSE), ]
  }
  
  top_peaks <- head(top_peaks, n_top)
  
  entrez <- unique(anno_df[[geneid_col]][anno_df[[peak_key_col]] %in% top_peaks$feature])
  entrez <- entrez[!is.na(entrez)]
  
  cat(label, "- top peaks:", nrow(top_peaks), "| mapped genes:", length(entrez), "\n")
  
  ego <- enrichGO(
    gene          = entrez,
    universe      = as.character(background_entrez),
    keyType       = "ENTREZID",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  ego_df <- if (is.null(ego)) data.frame() else as.data.frame(ego)
  
  write.csv(ego_df,
            file.path(mofa_out_dir, paste0("GO_mofa_", label, ".csv")),
            row.names = FALSE)
  
  if (!is.null(ego) && nrow(ego_df) > 0) {
    ggsave(
      file.path(mofa_out_dir, paste0("dotplot_GO_mofa_", label, ".png")),
      plot = dotplot(ego, showCategory = min(20, nrow(ego_df))) +
        labs(title = paste("GO -", label)) +
        theme(axis.text.y = element_text(size = 7),
              axis.text.x = element_text(size = 7),
              plot.title  = element_text(size = 10)),
      width = 10, height = 14, dpi = 150
    )
  } else {
    message("No significant GO terms for ", label)
  }
  
  return(ego)
}

# rGREAT ENRICHMENT FUNCTION
run_great_factor <- function(weights_df, view, factor_num, background_gr, label,
                             n_top = 500, mode = c("abs", "positive", "negative")) {
  
  mode <- match.arg(mode)
  
  if (mode == "abs") {
    top_peaks <- weights_df[order(abs(weights_df$value), decreasing = TRUE), ]
  } else if (mode == "positive") {
    top_peaks <- weights_df[weights_df$value > 0, ]
    top_peaks <- top_peaks[order(top_peaks$value, decreasing = TRUE), ]
  } else if (mode == "negative") {
    top_peaks <- weights_df[weights_df$value < 0, ]
    top_peaks <- top_peaks[order(top_peaks$value, decreasing = FALSE), ]
  }
  
  top_peaks <- head(top_peaks, n_top)
  
  if (view == "bulk") {
    peaks_gr <- bulk_peaks_to_gr(top_peaks$feature)
  } else {
    peaks_gr <- sc_peaks_to_gr(top_peaks$feature)
  }
  
  great_res <- great(peaks_gr,
                     gene_sets  = "GO:BP",
                     tss_source = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                     background = background_gr,
                     cores      = 8)
  
  tb     <- getEnrichmentTable(great_res)
  tb_sig <- tb[tb$p_adjust < 0.05, ]
  tb_sig <- tb_sig[order(tb_sig$p_adjust), ]
  
  cat(label, "- significant rGREAT terms:", nrow(tb_sig), "\n")
  
  write.csv(tb_sig[, c("id", "description", "fold_enrichment", "observed_region_hits", "p_adjust")],
            file.path(mofa_out_dir, paste0("rGREAT_", label, ".csv")),
            row.names = FALSE)
  
  return(tb_sig)
}

# FUNCTION TO GET TOP 20 ANNOTATED PEAKS PER FACTOR
get_top20 <- function(weights_df, anno_df, peak_key_col, gene_col, label, mode = "abs") {
  
  if (mode == "abs") {
    top_peaks <- weights_df[order(abs(weights_df$value), decreasing = TRUE), ]
  } else if (mode == "positive") {
    top_peaks <- weights_df[weights_df$value > 0, ]
    top_peaks <- top_peaks[order(top_peaks$value, decreasing = TRUE), ]
  } else if (mode == "negative") {
    top_peaks <- weights_df[weights_df$value < 0, ]
    top_peaks <- top_peaks[order(top_peaks$value, decreasing = FALSE), ]
  }
  
  top20 <- head(top_peaks, 20)
  
  annotated <- merge(top20,
                     anno_df[, c(peak_key_col, "annotation", gene_col, "distanceToTSS")],
                     by.x = "feature",
                     by.y = peak_key_col)
  
  annotated <- annotated[!duplicated(annotated$feature), ]
  
  write.csv(annotated[, c("feature", "value", "annotation", gene_col, "distanceToTSS")],
            file.path(mofa_out_dir, paste0("top20_", label, ".csv")),
            row.names = FALSE)
  
  cat("Saved top 20:", label, "\n")
  annotated
}

# MODEL 1: BULK VS OPALIN+
mofa_input_opalin <- list(bulk = as.matrix(bulk_hvf), opalin = as.matrix(opalin_hvf))
mofa_opalin_vst   <- create_mofa(mofa_input_opalin)

model_opts             <- get_default_model_options(mofa_opalin_vst)
model_opts$num_factors <- 9
train_opts             <- get_default_training_options(mofa_opalin_vst)
train_opts$convergence_mode <- "slow"
train_opts$seed        <- 42

mofa_opalin_vst <- prepare_mofa(mofa_opalin_vst,
                                data_options     = get_default_data_options(mofa_opalin_vst),
                                model_options    = model_opts,
                                training_options = train_opts)
mofa_opalin_vst <- run_mofa(mofa_opalin_vst,
                            outfile      = file.path(mofa_out_dir, "mofa_opalin_newsc.hdf5"),
                            use_basilisk = TRUE)

samples_metadata(mofa_opalin_vst) <- sample_meta
plot_variance_explained(mofa_opalin_vst, max_r2 = 15)
plot_factor_cor(mofa_opalin_vst)
plot_factor(mofa_opalin_vst, factors = 1:7, color_by = "condition")

# EXTRACT OPALIN MODEL WEIGHTS
# FACTOR 1 BULK, FACTORS 2+3 OPALIN, FACTOR 4 SHARED, FACTOR 5 BULK, FACTOR 6 OPALIN
bulk_weights_f1   <- get_weights(mofa_opalin_vst, views = "bulk",   factors = 1, as.data.frame = TRUE)
opalin_weights_f2 <- get_weights(mofa_opalin_vst, views = "opalin", factors = 2, as.data.frame = TRUE)
opalin_weights_f3 <- get_weights(mofa_opalin_vst, views = "opalin", factors = 3, as.data.frame = TRUE)
bulk_weights_f4   <- get_weights(mofa_opalin_vst, views = "bulk",   factors = 4, as.data.frame = TRUE)
opalin_weights_f4 <- get_weights(mofa_opalin_vst, views = "opalin", factors = 4, as.data.frame = TRUE)
bulk_weights_f5   <- get_weights(mofa_opalin_vst, views = "bulk",   factors = 5, as.data.frame = TRUE)
opalin_weights_f6 <- get_weights(mofa_opalin_vst, views = "opalin", factors = 6, as.data.frame = TRUE)

bulk_weights_f1$feature   <- as.character(bulk_weights_f1$feature)
opalin_weights_f2$feature <- as.character(opalin_weights_f2$feature)
opalin_weights_f3$feature <- as.character(opalin_weights_f3$feature)
bulk_weights_f4$feature   <- as.character(bulk_weights_f4$feature)
opalin_weights_f4$feature <- as.character(opalin_weights_f4$feature)
bulk_weights_f5$feature   <- as.character(bulk_weights_f5$feature)
opalin_weights_f6$feature <- as.character(opalin_weights_f6$feature)

# TOP 20 ANNOTATED PEAKS - OPALIN MODEL
get_top20(bulk_weights_f1,   bulk_anno_df,   "peak_key", "GENENAME",   "opalin_factor1_bulk")
get_top20(opalin_weights_f2, opalin_anno_cp, "peak_key", "GENENAME.x", "opalin_factor2_opalin")
get_top20(opalin_weights_f3, opalin_anno_cp, "peak_key", "GENENAME.x", "opalin_factor3_opalin")
get_top20(bulk_weights_f4,   bulk_anno_df,   "peak_key", "GENENAME",   "opalin_factor4_bulk")
get_top20(opalin_weights_f4, opalin_anno_cp, "peak_key", "GENENAME.x", "opalin_factor4_opalin")
get_top20(bulk_weights_f5,   bulk_anno_df,   "peak_key", "GENENAME",   "opalin_factor5_bulk")
get_top20(opalin_weights_f6, opalin_anno_cp, "peak_key", "GENENAME.x", "opalin_factor6_opalin")

# GO ENRICHMENT - OPALIN MODEL
ego_f1_bulk_pos   <- run_go_factor(bulk_weights_f1,   bulk_anno_df,   "peak_key", "geneId", "opalin_factor1_bulk_positive",   bulk_background_entrez, mode = "positive")
ego_f1_bulk_neg   <- run_go_factor(bulk_weights_f1,   bulk_anno_df,   "peak_key", "geneId", "opalin_factor1_bulk_negative",   bulk_background_entrez, mode = "negative")
ego_f2_opalin_pos <- run_go_factor(opalin_weights_f2, opalin_anno_cp, "peak_key", "geneId", "opalin_factor2_opalin_positive", sc_background_entrez,   mode = "positive")
ego_f2_opalin_neg <- run_go_factor(opalin_weights_f2, opalin_anno_cp, "peak_key", "geneId", "opalin_factor2_opalin_negative", sc_background_entrez,   mode = "negative")
ego_f3_opalin_pos <- run_go_factor(opalin_weights_f3, opalin_anno_cp, "peak_key", "geneId", "opalin_factor3_opalin_positive", sc_background_entrez,   mode = "positive")
ego_f3_opalin_neg <- run_go_factor(opalin_weights_f3, opalin_anno_cp, "peak_key", "geneId", "opalin_factor3_opalin_negative", sc_background_entrez,   mode = "negative")
ego_f4_bulk_pos   <- run_go_factor(bulk_weights_f4,   bulk_anno_df,   "peak_key", "geneId", "opalin_factor4_bulk_positive",   bulk_background_entrez, mode = "positive")
ego_f4_bulk_neg   <- run_go_factor(bulk_weights_f4,   bulk_anno_df,   "peak_key", "geneId", "opalin_factor4_bulk_negative",   bulk_background_entrez, mode = "negative")
ego_f4_opalin_pos <- run_go_factor(opalin_weights_f4, opalin_anno_cp, "peak_key", "geneId", "opalin_factor4_opalin_positive", sc_background_entrez,   mode = "positive")
ego_f4_opalin_neg <- run_go_factor(opalin_weights_f4, opalin_anno_cp, "peak_key", "geneId", "opalin_factor4_opalin_negative", sc_background_entrez,   mode = "negative")

cat("--- GO term counts - Opalin model ---\n")
cat("Factor 1 bulk positive:", nrow(data.frame(ego_f1_bulk_pos)), "terms\n")
cat("Factor 1 bulk negative:", nrow(data.frame(ego_f1_bulk_neg)), "terms\n")
cat("Factor 2 opalin positive:", nrow(data.frame(ego_f2_opalin_pos)), "terms\n")
cat("Factor 2 opalin negative:", nrow(data.frame(ego_f2_opalin_neg)), "terms\n")
cat("Factor 3 opalin positive:", nrow(data.frame(ego_f3_opalin_pos)), "terms\n")
cat("Factor 3 opalin negative:", nrow(data.frame(ego_f3_opalin_neg)), "terms\n")
cat("Factor 4 bulk positive:", nrow(data.frame(ego_f4_bulk_pos)), "terms\n")
cat("Factor 4 bulk negative:", nrow(data.frame(ego_f4_bulk_neg)), "terms\n")
cat("Factor 4 opalin positive:", nrow(data.frame(ego_f4_opalin_pos)), "terms\n")
cat("Factor 4 opalin negative:", nrow(data.frame(ego_f4_opalin_neg)), "terms\n")


opalin_weights_f2$feature <- as.character(opalin_weights_f2$feature)
opalin_weights_f3$feature <- as.character(opalin_weights_f3$feature)
opalin_weights_f4$feature <- as.character(opalin_weights_f4$feature)
bulk_weights_f1$feature   <- as.character(bulk_weights_f1$feature)
bulk_weights_f4$feature   <- as.character(bulk_weights_f4$feature)
bulk_weights_f5$feature   <- as.character(bulk_weights_f5$feature)
opalin_weights_f6$feature <- as.character(opalin_weights_f6$feature)

# rGREAT ENRICHMENT - OPALIN MODEL
great_f1_bulk_pos   <- run_great_factor(bulk_weights_f1,   "bulk",   1, bulk_consensus, "opalin_factor1_bulk_positive",   mode = "positive")
great_f2_opalin_neg <- run_great_factor(opalin_weights_f2, "opalin", 2, sc_consensus,   "opalin_factor2_opalin_negative", mode = "negative")
great_f3_opalin_pos <- run_great_factor(opalin_weights_f3, "opalin", 3, sc_consensus,   "opalin_factor3_opalin_positive", mode = "positive")
great_f4_bulk_pos   <- run_great_factor(bulk_weights_f4,   "bulk",   4, bulk_consensus, "opalin_factor4_bulk_positive",   mode = "positive")
great_f4_opalin_pos <- run_great_factor(opalin_weights_f4, "opalin", 4, sc_consensus,   "opalin_factor4_opalin_positive", mode = "positive")

saveRDS(mofa_opalin_vst, file.path(mofa_out_dir, "mofa_opalin_vst_object.rds"))

# MODEL 2: BULK VS PLEKHG1+
mofa_input_plekhg1 <- list(bulk = as.matrix(bulk_hvf), plekhg1 = as.matrix(plekhg1_hvf))
mofa_plekhg1_vst   <- create_mofa(mofa_input_plekhg1)

model_opts             <- get_default_model_options(mofa_plekhg1_vst)
model_opts$num_factors <- 9
train_opts             <- get_default_training_options(mofa_plekhg1_vst)
train_opts$convergence_mode <- "slow"
train_opts$seed        <- 42

mofa_plekhg1_vst <- prepare_mofa(mofa_plekhg1_vst,
                                 data_options     = get_default_data_options(mofa_plekhg1_vst),
                                 model_options    = model_opts,
                                 training_options = train_opts)
mofa_plekhg1_vst <- run_mofa(mofa_plekhg1_vst,
                             outfile      = file.path(mofa_out_dir, "mofa_plekhg1_newsc.hdf5"),
                             use_basilisk = TRUE)

samples_metadata(mofa_plekhg1_vst) <- sample_meta
plot_variance_explained(mofa_plekhg1_vst, max_r2 = 15)
plot_factor_cor(mofa_plekhg1_vst)
plot_factor(mofa_plekhg1_vst, factors = 1:7, color_by = "condition")

# EXTRACT PLEKHG1 MODEL WEIGHTS
# FACTOR 1 PLEKHG1, FACTORS 2+3+4 BULK, FACTOR 5 PLEKHG1
get_weights_plekhg1 <- function(factor_num, view) {
  w <- get_weights(mofa_plekhg1_vst, views = view, factors = factor_num, as.data.frame = TRUE)
  w$feature <- as.character(w$feature)
  w
}

plekhg1_weights_f1      <- get_weights_plekhg1(1, "plekhg1")
bulk_weights_plekhg1_f2 <- get_weights_plekhg1(2, "bulk")
bulk_weights_plekhg1_f3 <- get_weights_plekhg1(3, "bulk")
bulk_weights_plekhg1_f4 <- get_weights_plekhg1(4, "bulk")
plekhg1_weights_f5      <- get_weights_plekhg1(5, "plekhg1")

# TOP 20 ANNOTATED PEAKS - PLEKHG1 MODEL
get_top20(plekhg1_weights_f1,      plekhg1_anno_cp, "peak_key", "GENENAME.x", "plekhg1_factor1_plekhg1")
get_top20(bulk_weights_plekhg1_f2, bulk_anno_df,    "peak_key", "GENENAME",   "plekhg1_factor2_bulk")
get_top20(bulk_weights_plekhg1_f3, bulk_anno_df,    "peak_key", "GENENAME",   "plekhg1_factor3_bulk")
get_top20(bulk_weights_plekhg1_f4, bulk_anno_df,    "peak_key", "GENENAME",   "plekhg1_factor4_bulk")
get_top20(plekhg1_weights_f5,      plekhg1_anno_cp, "peak_key", "GENENAME.x", "plekhg1_factor5_plekhg1")

# GO ENRICHMENT - PLEKHG1 MODEL
ego_plekhg1_f1_pos <- run_go_factor(plekhg1_weights_f1,      plekhg1_anno_cp, "peak_key", "geneId", "plekhg1_factor1_plekhg1_positive", sc_background_entrez,   mode = "positive")
ego_plekhg1_f1_neg <- run_go_factor(plekhg1_weights_f1,      plekhg1_anno_cp, "peak_key", "geneId", "plekhg1_factor1_plekhg1_negative", sc_background_entrez,   mode = "negative")
ego_plekhg1_f2_pos <- run_go_factor(bulk_weights_plekhg1_f2, bulk_anno_df,    "peak_key", "geneId", "plekhg1_factor2_bulk_positive",    bulk_background_entrez, mode = "positive")
ego_plekhg1_f2_neg <- run_go_factor(bulk_weights_plekhg1_f2, bulk_anno_df,    "peak_key", "geneId", "plekhg1_factor2_bulk_negative",    bulk_background_entrez, mode = "negative")
ego_plekhg1_f3_pos <- run_go_factor(bulk_weights_plekhg1_f3, bulk_anno_df,    "peak_key", "geneId", "plekhg1_factor3_bulk_positive",    bulk_background_entrez, mode = "positive")
ego_plekhg1_f3_neg <- run_go_factor(bulk_weights_plekhg1_f3, bulk_anno_df,    "peak_key", "geneId", "plekhg1_factor3_bulk_negative",    bulk_background_entrez, mode = "negative")
ego_plekhg1_f4_pos <- run_go_factor(bulk_weights_plekhg1_f4, bulk_anno_df,    "peak_key", "geneId", "plekhg1_factor4_bulk_positive",    bulk_background_entrez, mode = "positive")
ego_plekhg1_f4_neg <- run_go_factor(bulk_weights_plekhg1_f4, bulk_anno_df,    "peak_key", "geneId", "plekhg1_factor4_bulk_negative",    bulk_background_entrez, mode = "negative")
ego_plekhg1_f5_pos <- run_go_factor(plekhg1_weights_f5,      plekhg1_anno_cp, "peak_key", "geneId", "plekhg1_factor5_plekhg1_positive", sc_background_entrez,   mode = "positive")
ego_plekhg1_f5_neg <- run_go_factor(plekhg1_weights_f5,      plekhg1_anno_cp, "peak_key", "geneId", "plekhg1_factor5_plekhg1_negative", sc_background_entrez,   mode = "negative")

cat("--- GO term counts - Plekhg1 model ---\n")
cat("Factor 1 plekhg1 positive:", nrow(data.frame(ego_plekhg1_f1_pos)), "terms\n")
cat("Factor 1 plekhg1 negative:", nrow(data.frame(ego_plekhg1_f1_neg)), "terms\n")
cat("Factor 2 bulk positive:", nrow(data.frame(ego_plekhg1_f2_pos)), "terms\n")
cat("Factor 2 bulk negative:", nrow(data.frame(ego_plekhg1_f2_neg)), "terms\n")
cat("Factor 3 bulk positive:", nrow(data.frame(ego_plekhg1_f3_pos)), "terms\n")
cat("Factor 3 bulk negative:", nrow(data.frame(ego_plekhg1_f3_neg)), "terms\n")
cat("Factor 4 bulk positive:", nrow(data.frame(ego_plekhg1_f4_pos)), "terms\n")
cat("Factor 4 bulk negative:", nrow(data.frame(ego_plekhg1_f4_neg)), "terms\n")
cat("Factor 5 plekhg1 positive:", nrow(data.frame(ego_plekhg1_f5_pos)), "terms\n")
cat("Factor 5 plekhg1 negative:", nrow(data.frame(ego_plekhg1_f5_neg)), "terms\n")

# rGREAT ENRICHMENT - PLEKHG1 MODEL
great_plekhg1_f1 <- run_great_factor(plekhg1_weights_f1,      "plekhg1", 1, sc_consensus,   "plekhg1_factor1_plekhg1_positive", mode = "positive")
great_plekhg1_f2 <- run_great_factor(bulk_weights_plekhg1_f2, "bulk",    2, bulk_consensus, "plekhg1_factor2_bulk_positive",    mode = "positive")
great_plekhg1_f3 <- run_great_factor(bulk_weights_plekhg1_f3, "bulk",    3, bulk_consensus, "plekhg1_factor3_bulk_positive",    mode = "positive")
great_plekhg1_f4 <- run_great_factor(bulk_weights_plekhg1_f4, "bulk",    4, bulk_consensus, "plekhg1_factor4_bulk_positive",    mode = "positive")
great_plekhg1_f5 <- run_great_factor(plekhg1_weights_f5,      "plekhg1", 5, sc_consensus,   "plekhg1_factor5_plekhg1_positive", mode = "positive")

saveRDS(mofa_plekhg1_vst, file.path(mofa_out_dir, "mofa_plekhg1_vst_object.rds"))

# MODEL 3: BULK VS OPCs
mofa_input_opc <- list(bulk = as.matrix(bulk_hvf), opc = as.matrix(opc_hvf))
mofa_opc_vst   <- create_mofa(mofa_input_opc)

model_opts             <- get_default_model_options(mofa_opc_vst)
model_opts$num_factors <- 9
train_opts             <- get_default_training_options(mofa_opc_vst)
train_opts$convergence_mode <- "slow"
train_opts$seed        <- 42

mofa_opc_vst <- prepare_mofa(mofa_opc_vst,
                             data_options     = get_default_data_options(mofa_opc_vst),
                             model_options    = model_opts,
                             training_options = train_opts)
mofa_opc_vst <- run_mofa(mofa_opc_vst,
                         outfile      = file.path(mofa_out_dir, "mofa_opc_newsc.hdf5"),
                         use_basilisk = TRUE)

samples_metadata(mofa_opc_vst) <- sample_meta
plot_variance_explained(mofa_opc_vst, max_r2 = 15)
plot_factor_cor(mofa_opc_vst)
plot_factor(mofa_opc_vst, factors = 1:7, color_by = "condition")

# EXTRACT OPC MODEL WEIGHTS
# FACTOR 1 BULK, FACTORS 2+3 OPC, FACTOR 4+5 BULK - UPDATE AFTER REVIEWING VARIANCE PLOT
get_weights_opc <- function(factor_num, view) {
  w <- get_weights(mofa_opc_vst, views = view, factors = factor_num, as.data.frame = TRUE)
  w$feature <- as.character(w$feature)
  w
}

bulk_weights_opc_f1 <- get_weights_opc(1, "bulk")
opc_weights_f2      <- get_weights_opc(2, "opc")
opc_weights_f3      <- get_weights_opc(3, "opc")
bulk_weights_opc_f4 <- get_weights_opc(4, "bulk")
bulk_weights_opc_f5 <- get_weights_opc(5, "bulk")

# VERIFY OPC WEIGHTS MATCH ANNOTATION
cat("OPC weights in annotation:", sum(opc_weights_f2$feature %in% opc_anno_cp$peak_key), "/ 5000\n")

opc_anno_cp$peak_key <- paste(gsub("chr", "", opc_anno_cp$seqnames),
                              opc_anno_cp$start,
                              opc_anno_cp$end, sep = "-")
# TOP 20 ANNOTATED PEAKS - OPC MODEL
get_top20(bulk_weights_opc_f1, bulk_anno_df, "peak_key", "GENENAME",   "opc_factor1_bulk")
get_top20(opc_weights_f2,      opc_anno_cp,  "peak_key", "GENENAME.x", "opc_factor2_opc")
get_top20(opc_weights_f3,      opc_anno_cp,  "peak_key", "GENENAME.x", "opc_factor3_opc")
get_top20(bulk_weights_opc_f4, bulk_anno_df, "peak_key", "GENENAME",   "opc_factor4_bulk")
get_top20(bulk_weights_opc_f5, bulk_anno_df, "peak_key", "GENENAME",   "opc_factor5_bulk")

# GO ENRICHMENT - OPC MODEL
ego_opc_f1_pos <- run_go_factor(bulk_weights_opc_f1, bulk_anno_df, "peak_key", "geneId", "opc_factor1_bulk_positive", bulk_background_entrez, mode = "positive")
ego_opc_f1_neg <- run_go_factor(bulk_weights_opc_f1, bulk_anno_df, "peak_key", "geneId", "opc_factor1_bulk_negative", bulk_background_entrez, mode = "negative")
ego_opc_f2_pos <- run_go_factor(opc_weights_f2,      opc_anno_cp,  "peak_key", "geneId", "opc_factor2_opc_positive",  sc_background_entrez,   mode = "positive")
ego_opc_f2_neg <- run_go_factor(opc_weights_f2,      opc_anno_cp,  "peak_key", "geneId", "opc_factor2_opc_negative",  sc_background_entrez,   mode = "negative")
ego_opc_f3_pos <- run_go_factor(opc_weights_f3,      opc_anno_cp,  "peak_key", "geneId", "opc_factor3_opc_positive",  sc_background_entrez,   mode = "positive")
ego_opc_f3_neg <- run_go_factor(opc_weights_f3,      opc_anno_cp,  "peak_key", "geneId", "opc_factor3_opc_negative",  sc_background_entrez,   mode = "negative")
ego_opc_f4_pos <- run_go_factor(bulk_weights_opc_f4, bulk_anno_df, "peak_key", "geneId", "opc_factor4_bulk_positive", bulk_background_entrez, mode = "positive")
ego_opc_f4_neg <- run_go_factor(bulk_weights_opc_f4, bulk_anno_df, "peak_key", "geneId", "opc_factor4_bulk_negative", bulk_background_entrez, mode = "negative")
ego_opc_f5_pos <- run_go_factor(bulk_weights_opc_f5, bulk_anno_df, "peak_key", "geneId", "opc_factor5_bulk_positive", bulk_background_entrez, mode = "positive")
ego_opc_f5_neg <- run_go_factor(bulk_weights_opc_f5, bulk_anno_df, "peak_key", "geneId", "opc_factor5_bulk_negative", bulk_background_entrez, mode = "negative")

cat("--- GO term counts - OPC model ---\n")
cat("Factor 1 bulk positive:", nrow(data.frame(ego_opc_f1_pos)), "terms\n")
cat("Factor 1 bulk negative:", nrow(data.frame(ego_opc_f1_neg)), "terms\n")
cat("Factor 2 opc positive:", nrow(data.frame(ego_opc_f2_pos)), "terms\n")
cat("Factor 2 opc negative:", nrow(data.frame(ego_opc_f2_neg)), "terms\n")
cat("Factor 3 opc positive:", nrow(data.frame(ego_opc_f3_pos)), "terms\n")
cat("Factor 3 opc negative:", nrow(data.frame(ego_opc_f3_neg)), "terms\n")
cat("Factor 4 bulk positive:", nrow(data.frame(ego_opc_f4_pos)), "terms\n")
cat("Factor 4 bulk negative:", nrow(data.frame(ego_opc_f4_neg)), "terms\n")
cat("Factor 5 bulk positive:", nrow(data.frame(ego_opc_f5_pos)), "terms\n")
cat("Factor 5 bulk negative:", nrow(data.frame(ego_opc_f5_neg)), "terms\n")

# rGREAT ENRICHMENT - OPC MODEL
great_opc_f1 <- run_great_factor(bulk_weights_opc_f1, "bulk", 1, bulk_consensus, "opc_factor1_bulk_positive", mode = "positive")
great_opc_f2 <- run_great_factor(opc_weights_f2,      "opc",  2, sc_consensus,   "opc_factor2_opc_positive",  mode = "positive")
great_opc_f3 <- run_great_factor(opc_weights_f3,      "opc",  3, sc_consensus,   "opc_factor3_opc_positive",  mode = "positive")
great_opc_f4 <- run_great_factor(bulk_weights_opc_f4, "bulk", 4, bulk_consensus, "opc_factor4_bulk_positive", mode = "positive")
great_opc_f5 <- run_great_factor(bulk_weights_opc_f5, "bulk", 5, bulk_consensus, "opc_factor5_bulk_positive", mode = "positive")

saveRDS(mofa_opc_vst, file.path(mofa_out_dir, "mofa_opc_vst_object.rds"))

cat("All results saved to:", mofa_out_dir, "\n")


# MODEL 4: 4-VIEW MOFA (BULK + ALL THREE SC SUBTYPES)
mofa_input_4view <- list(
  bulk    = as.matrix(bulk_hvf),
  opalin  = as.matrix(opalin_hvf),
  plekhg1 = as.matrix(plekhg1_hvf),
  opc     = as.matrix(opc_hvf)
)

mofa_4view_vst <- create_mofa(mofa_input_4view)

model_opts             <- get_default_model_options(mofa_4view_vst)
model_opts$num_factors <- 9
train_opts             <- get_default_training_options(mofa_4view_vst)
train_opts$convergence_mode <- "slow"
train_opts$seed        <- 42

mofa_4view_vst <- prepare_mofa(mofa_4view_vst,
                               data_options     = get_default_data_options(mofa_4view_vst),
                               model_options    = model_opts,
                               training_options = train_opts)

mofa_4view_vst <- run_mofa(mofa_4view_vst,
                           outfile      = file.path(mofa_out_dir, "mofa_4view_newsc.hdf5"),
                           use_basilisk = TRUE)

samples_metadata(mofa_4view_vst) <- sample_meta

plot_variance_explained(mofa_4view_vst, max_r2 = 15)
plot_factor_cor(mofa_4view_vst)
plot_factor(mofa_4view_vst, factors = 1:7, color_by = "condition")

# RECALCULATE OPC HVF FROM CURRENT opc_vst
opc_vars    <- apply(opc_vst, 1, var)
opc_top_idx <- head(order(opc_vars, decreasing = TRUE), 5000)
opc_hvf     <- opc_vst[opc_top_idx, ]

# VERIFY
cat("OPC HVF peaks in annotation:", sum(rownames(opc_hvf) %in% opc_anno_cp$peak_key), "\n")
cat("OPC HVF peaks in opc_vst:", sum(rownames(opc_hvf) %in% rownames(opc_vst)), "\n")

# EXTRACT 4-VIEW MODEL WEIGHTS
get_weights_4view <- function(factor_num, view) {
  w <- get_weights(mofa_4view_vst, views = view, factors = factor_num, as.data.frame = TRUE)
  w$feature <- as.character(w$feature)
  w
}

# FACTOR 1 - SC SHARED (ALL THREE SC SUBTYPES)
opalin_weights_4v_f1  <- get_weights_4view(1, "opalin")
plekhg1_weights_4v_f1 <- get_weights_4view(1, "plekhg1")
opc_weights_4v_f1     <- get_weights_4view(1, "opc")

# FACTOR 2 - BULK SPECIFIC
bulk_weights_4v_f2 <- get_weights_4view(2, "bulk")

# FACTOR 3 - OPALIN AND OPC SHARED
opalin_weights_4v_f3 <- get_weights_4view(3, "opalin")
opc_weights_4v_f3    <- get_weights_4view(3, "opc")

# FACTOR 4 - BULK SPECIFIC
bulk_weights_4v_f4 <- get_weights_4view(4, "bulk")

# TOP 20 ANNOTATED PEAKS
get_top20(opalin_weights_4v_f1,  opalin_anno_cp,  "peak_key", "GENENAME.x", "4view_factor1_opalin")
get_top20(plekhg1_weights_4v_f1, plekhg1_anno_cp, "peak_key", "GENENAME.x", "4view_factor1_plekhg1")
get_top20(opc_weights_4v_f1,     opc_anno_cp,     "peak_key", "GENENAME.x", "4view_factor1_opc")
get_top20(bulk_weights_4v_f2,    bulk_anno_df,    "peak_key", "GENENAME",   "4view_factor2_bulk")
get_top20(opalin_weights_4v_f3,  opalin_anno_cp,  "peak_key", "GENENAME.x", "4view_factor3_opalin")
get_top20(opc_weights_4v_f3,     opc_anno_cp,     "peak_key", "GENENAME.x", "4view_factor3_opc")
get_top20(bulk_weights_4v_f4,    bulk_anno_df,    "peak_key", "GENENAME",   "4view_factor4_bulk")

# GO ENRICHMENT
ego_4v_f1_opalin  <- run_go_factor(opalin_weights_4v_f1,  opalin_anno_cp,  "peak_key", "geneId", "4view_factor1_opalin_positive",  sc_background_entrez,   mode = "positive")
ego_4v_f1_plekhg1 <- run_go_factor(plekhg1_weights_4v_f1, plekhg1_anno_cp, "peak_key", "geneId", "4view_factor1_plekhg1_positive", sc_background_entrez,   mode = "positive")
ego_4v_f1_opc     <- run_go_factor(opc_weights_4v_f1,     opc_anno_cp,     "peak_key", "geneId", "4view_factor1_opc_positive",     sc_background_entrez,   mode = "positive")
ego_4v_f2_bulk    <- run_go_factor(bulk_weights_4v_f2,    bulk_anno_df,    "peak_key", "geneId", "4view_factor2_bulk_positive",    bulk_background_entrez, mode = "positive")
ego_4v_f3_opalin  <- run_go_factor(opalin_weights_4v_f3,  opalin_anno_cp,  "peak_key", "geneId", "4view_factor3_opalin_positive",  sc_background_entrez,   mode = "positive")
ego_4v_f3_opc     <- run_go_factor(opc_weights_4v_f3,     opc_anno_cp,     "peak_key", "geneId", "4view_factor3_opc_positive",     sc_background_entrez,   mode = "positive")
ego_4v_f4_bulk    <- run_go_factor(bulk_weights_4v_f4,    bulk_anno_df,    "peak_key", "geneId", "4view_factor4_bulk_positive",    bulk_background_entrez, mode = "positive")

cat("--- GO term counts - 4-view model ---\n")
cat("Factor 1 opalin positive:", nrow(data.frame(ego_4v_f1_opalin)), "terms\n")
cat("Factor 1 plekhg1 positive:", nrow(data.frame(ego_4v_f1_plekhg1)), "terms\n")
cat("Factor 1 opc positive:", nrow(data.frame(ego_4v_f1_opc)), "terms\n")
cat("Factor 2 bulk positive:", nrow(data.frame(ego_4v_f2_bulk)), "terms\n")
cat("Factor 3 opalin positive:", nrow(data.frame(ego_4v_f3_opalin)), "terms\n")
cat("Factor 3 opc positive:", nrow(data.frame(ego_4v_f3_opc)), "terms\n")
cat("Factor 4 bulk positive:", nrow(data.frame(ego_4v_f4_bulk)), "terms\n")

# rGREAT ENRICHMENT
great_4v_f1_opalin  <- run_great_factor(opalin_weights_4v_f1,  "opalin",  1, sc_consensus,   "4view_factor1_opalin_positive",  mode = "positive")
great_4v_f1_plekhg1 <- run_great_factor(plekhg1_weights_4v_f1, "plekhg1", 1, sc_consensus,   "4view_factor1_plekhg1_positive", mode = "positive")
great_4v_f1_opc     <- run_great_factor(opc_weights_4v_f1,     "opc",     1, sc_consensus,   "4view_factor1_opc_positive",     mode = "positive")
great_4v_f2_bulk    <- run_great_factor(bulk_weights_4v_f2,    "bulk",    2, bulk_consensus, "4view_factor2_bulk_positive",    mode = "positive")
great_4v_f3_opalin  <- run_great_factor(opalin_weights_4v_f3,  "opalin",  3, sc_consensus,   "4view_factor3_opalin_positive",  mode = "positive")
great_4v_f3_opc     <- run_great_factor(opc_weights_4v_f3,     "opc",     3, sc_consensus,   "4view_factor3_opc_positive",     mode = "positive")
great_4v_f4_bulk    <- run_great_factor(bulk_weights_4v_f4,    "bulk",    4, bulk_consensus, "4view_factor4_bulk_positive",    mode = "positive")

saveRDS(mofa_4view_vst, file.path(mofa_out_dir, "mofa_4view_vst_object.rds"))
cat("All 4-view results saved to:", mofa_out_dir, "\n")


# CHECK IF OPALIN AND OPC HVF PEAKS ARE DIFFERENT
cat("Opalin HVF peaks:", nrow(opalin_hvf), "\n")
cat("OPC HVF peaks:", nrow(opc_hvf), "\n")

# HOW MANY PEAKS ARE SHARED BETWEEN OPALIN AND OPC HVF
shared_hvf <- intersect(rownames(opalin_hvf), rownames(opc_hvf))
cat("Shared HVF peaks between Opalin and OPC:", length(shared_hvf), "\n")

# CHECK DIMENSIONS OF VST MATRICES
dim(opalin_vst)
dim(opc_vst)

# CHECK IF ROWNAMES ARE IDENTICAL
cat("Identical rownames:", identical(rownames(opalin_vst), rownames(opc_vst)), "\n")





