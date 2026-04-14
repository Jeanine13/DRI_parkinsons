# rGREAT ENRICHMENT SCRIPT
# RUN AFTER cluster_peaks_MOFA.R
# LOADS SAVED WEIGHTS AND RUNS rGREAT ENRICHMENT FOR ALL MODELS AND FACTORS
# USES POSITIVE AND NEGATIVE WEIGHTS FOR EACH FACTOR/VIEW COMBINATION
# BACKGROUND: SC CONSENSUS PEAKS FOR SC VIEWS, BULK UNION PEAKS FOR BULK VIEWS
# DOTPLOTS SAVED AS PDFs USING OKABE-ITO COLOUR BLIND FRIENDLY PALETTE

.libPaths(c(
  "/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5",
  .libPaths()
))

library(rGREAT)
library(GenomicRanges)
library(GenomeInfoDb)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)

mofa_out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/new_sc_mofa_clean_14042026"
rgreat_out_dir <- file.path(mofa_out_dir, "rgreat_results")
dir.create(rgreat_out_dir, recursive = TRUE, showWarnings = FALSE)

# LOAD SAVED WEIGHTS FROM MOFA SCRIPT
load(file.path(mofa_out_dir, "all_mofa_weights.RData"))
cat("All weights loaded\n")

# LOAD BACKGROUND CONSENSUS PEAK SETS
sc_consensus   <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/single_cell_consensus_peaks.rds")
bulk_consensus <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/union_peaks/bulk_union_peaks.rds")

seqlevelsStyle(sc_consensus)   <- "UCSC"
seqlevelsStyle(bulk_consensus) <- "UCSC"

cat("Background peak sets loaded\n")
cat("SC consensus peaks:", length(sc_consensus), "\n")
cat("Bulk consensus peaks:", length(bulk_consensus), "\n")

# FUNCTION TO CONVERT BULK PEAK NAMES TO GRANGES
bulk_peaks_to_gr <- function(peaks) {
  starts <- as.numeric(gsub(".*_(\\d+)_\\d+$", "\\1", peaks))
  ends   <- as.numeric(gsub(".*_(\\d+)$", "\\1", peaks))
  
  if (any(is.na(starts)) || any(is.na(ends))) {
    stop("Failed to parse one or more bulk peak names.")
  }
  
  GRanges(
    seqnames = gsub("_.*", "", peaks),
    ranges   = IRanges(start = starts, end = ends)
  )
}

# FUNCTION TO CONVERT SC PEAK NAMES TO GRANGES
sc_peaks_to_gr <- function(peaks) {
  split_peaks <- strsplit(as.character(peaks), "-")
  seqs   <- paste0("chr", sapply(split_peaks, `[`, 1))
  starts <- as.numeric(sapply(split_peaks, `[`, 2))
  ends   <- as.numeric(sapply(split_peaks, `[`, 3))
  
  if (any(is.na(starts)) || any(is.na(ends))) {
    stop("Failed to parse one or more single-cell peak names.")
  }
  
  GRanges(
    seqnames = seqs,
    ranges   = IRanges(start = starts, end = ends)
  )
}

# DOTPLOT FUNCTION USING OKABE-ITO COLOUR BLIND FRIENDLY PALETTE
# BLUE (#0072B2) FOR LOW P.ADJUST (MOST SIGNIFICANT)
# AMBER (#E69F00) FOR HIGH P.ADJUST (LESS SIGNIFICANT)
plot_great_dotplot <- function(tb_sig, label, mode, n_show = 20) {
  if (is.null(tb_sig) || nrow(tb_sig) == 0) {
    message("No significant terms for ", label, " - ", mode, " - skipping plot")
    return(NULL)
  }
  
  plot_df <- head(tb_sig, n_show)
  plot_df$description <- factor(plot_df$description, levels = rev(plot_df$description))
  
  p <- ggplot(
    plot_df,
    aes(
      x = fold_enrichment,
      y = description,
      size = observed_region_hits,
      colour = p_adjust
    )
  ) +
    geom_point() +
    scale_colour_gradient(
      low  = "#0072B2",
      high = "#E69F00",
      name = "p.adjust"
    ) +
    scale_size_continuous(name = "Region hits") +
    labs(
      title = paste("rGREAT GO:BP -", label, "-", mode),
      x     = "Fold enrichment",
      y     = ""
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title  = element_text(size = 10)
    )
  
  ggsave(
    file.path(rgreat_out_dir, paste0("dotplot_rGREAT_", label, "_", mode, ".pdf")),
    plot = p,
    width = 10,
    height = 8
  )
  
  return(p)
}

# rGREAT ENRICHMENT FUNCTION
run_great_factor <- function(weights_df, view, label, background_gr,
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
  
  if (nrow(top_peaks) == 0) {
    cat("No peaks for", label, "-", mode, "- skipping\n")
    return(NULL)
  }
  
  if (view == "bulk") {
    peaks_gr <- bulk_peaks_to_gr(top_peaks$feature)
  } else {
    peaks_gr <- sc_peaks_to_gr(top_peaks$feature)
  }
  
  cat(label, "-", mode, "- input peaks used:", length(peaks_gr), "\n")
  
  great_res <- great(
    peaks_gr,
    gene_sets  = "GO:BP",
    tss_source = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    background = background_gr,
    cores      = 4
  )
  
  tb <- getEnrichmentTable(great_res)
  tb_sig <- tb[tb$p_adjust < 0.05, ]
  tb_sig <- tb_sig[order(tb_sig$p_adjust), ]
  
  cat(label, "-", mode, "- significant terms:", nrow(tb_sig), "\n")
  
  write.csv(
    tb_sig[, c("id", "description", "fold_enrichment", "observed_region_hits", "p_adjust")],
    file.path(rgreat_out_dir, paste0("rGREAT_", label, "_", mode, ".csv")),
    row.names = FALSE
  )
  
  plot_great_dotplot(tb_sig, label, mode)
  
  return(tb_sig)
}

# HELPER TO RUN BOTH POSITIVE AND NEGATIVE FOR A GIVEN FACTOR/VIEW
run_great_both <- function(weights_df, view, label, background_gr, n_top = 500) {
  run_great_factor(weights_df, view, label, background_gr, n_top = n_top, mode = "positive")
  run_great_factor(weights_df, view, label, background_gr, n_top = n_top, mode = "negative")
}

# =========================
# MODEL 1: BULK VS OPALIN+
# =========================
cat("\n=== rGREAT - Opalin model ===\n")
run_great_both(bulk_weights_f1,   "bulk",   "opalin_factor1_bulk",   bulk_consensus)
run_great_both(opalin_weights_f2, "opalin", "opalin_factor2_opalin", sc_consensus)
run_great_both(opalin_weights_f3, "opalin", "opalin_factor3_opalin", sc_consensus)
run_great_both(bulk_weights_f4,   "bulk",   "opalin_factor4_bulk",   bulk_consensus)
run_great_both(opalin_weights_f4, "opalin", "opalin_factor4_opalin", sc_consensus)
run_great_both(bulk_weights_f5,   "bulk",   "opalin_factor5_bulk",   bulk_consensus)
run_great_both(opalin_weights_f6, "opalin", "opalin_factor6_opalin", sc_consensus)
run_great_both(bulk_weights_f7,   "bulk",   "opalin_factor7_bulk",   bulk_consensus)
run_great_both(opalin_weights_f7, "opalin", "opalin_factor7_opalin", sc_consensus)

# ==========================
# MODEL 2: BULK VS PLEKHG1+
# ==========================
cat("\n=== rGREAT - Plekhg1 model ===\n")
run_great_both(plekhg1_weights_f1,      "plekhg1", "plekhg1_factor1_plekhg1", sc_consensus)
run_great_both(bulk_weights_plekhg1_f2, "bulk",    "plekhg1_factor2_bulk",    bulk_consensus)
run_great_both(bulk_weights_plekhg1_f3, "bulk",    "plekhg1_factor3_bulk",    bulk_consensus)
run_great_both(bulk_weights_plekhg1_f4, "bulk",    "plekhg1_factor4_bulk",    bulk_consensus)
run_great_both(plekhg1_weights_f5,      "plekhg1", "plekhg1_factor5_plekhg1", sc_consensus)
run_great_both(bulk_weights_plekhg1_f6, "bulk",    "plekhg1_factor6_bulk",    bulk_consensus)
run_great_both(plekhg1_weights_f6,      "plekhg1", "plekhg1_factor6_plekhg1", sc_consensus)
run_great_both(bulk_weights_plekhg1_f7, "bulk",    "plekhg1_factor7_bulk",    bulk_consensus)
run_great_both(plekhg1_weights_f7,      "plekhg1", "plekhg1_factor7_plekhg1", sc_consensus)

# =====================
# MODEL 3: BULK VS OPC
# =====================
cat("\n=== rGREAT - OPC model ===\n")
run_great_both(bulk_weights_opc_f1, "bulk", "opc_factor1_bulk", bulk_consensus)
run_great_both(opc_weights_f2,      "opc",  "opc_factor2_opc",  sc_consensus)
run_great_both(opc_weights_f3,      "opc",  "opc_factor3_opc",  sc_consensus)
run_great_both(bulk_weights_opc_f4, "bulk", "opc_factor4_bulk", bulk_consensus)
run_great_both(bulk_weights_opc_f5, "bulk", "opc_factor5_bulk", bulk_consensus)
run_great_both(bulk_weights_opc_f6, "bulk", "opc_factor6_bulk", bulk_consensus)
run_great_both(opc_weights_f6,      "opc",  "opc_factor6_opc",  sc_consensus)
run_great_both(bulk_weights_opc_f7, "bulk", "opc_factor7_bulk", bulk_consensus)
run_great_both(opc_weights_f7,      "opc",  "opc_factor7_opc",  sc_consensus)

# ==========================================
# MODEL 4: 4-VIEW (BULK + ALL THREE SUBTYPES)
# ==========================================
cat("\n=== rGREAT - 4-view model ===\n")

# F1 SC SHARED - ALL THREE SC SUBTYPES
run_great_both(opalin_weights_4v_f1,  "opalin",  "4view_factor1_opalin",  sc_consensus)
run_great_both(plekhg1_weights_4v_f1, "plekhg1", "4view_factor1_plekhg1", sc_consensus)
run_great_both(opc_weights_4v_f1,     "opc",     "4view_factor1_opc",     sc_consensus)

# F2 BULK SPECIFIC
run_great_both(bulk_weights_4v_f2, "bulk", "4view_factor2_bulk", bulk_consensus)

# F3 OPALIN AND OPC SHARED, PLEKHG1 INCLUDED FOR COMPLETENESS
run_great_both(opalin_weights_4v_f3,  "opalin",  "4view_factor3_opalin",  sc_consensus)
run_great_both(opc_weights_4v_f3,     "opc",     "4view_factor3_opc",     sc_consensus)
run_great_both(plekhg1_weights_4v_f3, "plekhg1", "4view_factor3_plekhg1", sc_consensus)

# F4 BULK AND OPALIN
run_great_both(bulk_weights_4v_f4,   "bulk",   "4view_factor4_bulk",   bulk_consensus)
run_great_both(opalin_weights_4v_f4, "opalin", "4view_factor4_opalin", sc_consensus)

# F5 WEAK MIXED - ALL VIEWS
run_great_both(bulk_weights_4v_f5,    "bulk",    "4view_factor5_bulk",    bulk_consensus)
run_great_both(opalin_weights_4v_f5,  "opalin",  "4view_factor5_opalin",  sc_consensus)
run_great_both(plekhg1_weights_4v_f5, "plekhg1", "4view_factor5_plekhg1", sc_consensus)
run_great_both(opc_weights_4v_f5,     "opc",     "4view_factor5_opc",     sc_consensus)

# F6 WEAK - ALL VIEWS
run_great_both(bulk_weights_4v_f6,    "bulk",    "4view_factor6_bulk",    bulk_consensus)
run_great_both(opalin_weights_4v_f6,  "opalin",  "4view_factor6_opalin",  sc_consensus)
run_great_both(plekhg1_weights_4v_f6, "plekhg1", "4view_factor6_plekhg1", sc_consensus)
run_great_both(opc_weights_4v_f6,     "opc",     "4view_factor6_opc",     sc_consensus)

# F7 WEAK - ALL VIEWS
run_great_both(bulk_weights_4v_f7,    "bulk",    "4view_factor7_bulk",    bulk_consensus)
run_great_both(opalin_weights_4v_f7,  "opalin",  "4view_factor7_opalin",  sc_consensus)
run_great_both(plekhg1_weights_4v_f7, "plekhg1", "4view_factor7_plekhg1", sc_consensus)
run_great_both(opc_weights_4v_f7,     "opc",     "4view_factor7_opc",     sc_consensus)

cat("\nAll rGREAT enrichment complete\n")
cat("Results saved to:", rgreat_out_dir, "\n")