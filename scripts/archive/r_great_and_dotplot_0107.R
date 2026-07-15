# RGREAT FOR OPALIN FACTOR 2 BULK VIEW
# RUNS ENRICHMENT AND SAVES DOTPLOTS

.libPaths(c(
  "/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5",
  .libPaths()
))

library(rGREAT)
library(GenomicRanges)
library(GenomeInfoDb)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(stringr)

mofa_out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/new_sc_mofa_clean_14042026"
rgreat_dir   <- file.path(mofa_out_dir, "rgreat_results")
dotplot_dir  <- file.path(mofa_out_dir, "dotplots_top10")
dir.create(dotplot_dir, recursive = TRUE, showWarnings = FALSE)

# LOAD WEIGHTS
load(file.path(mofa_out_dir, "all_mofa_weights.RData"))
cat("Weights loaded\n")

# LOAD BACKGROUND
bulk_consensus <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/union_peaks/bulk_union_peaks.rds")
seqlevelsStyle(bulk_consensus) <- "UCSC"
cat("Background loaded\n")

# CONVERT BULK PEAK NAMES TO GRANGES
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

# DOTPLOT FUNCTION - TOP 10 TERMS
plot_great_dotplot <- function(tb_sig, label, mode, n_show = 10) {
  
  if (is.null(tb_sig) || nrow(tb_sig) == 0) {
    message("No significant terms for ", label, " - ", mode, " - skipping")
    return(NULL)
  }
  
  plot_df <- head(tb_sig, n_show)
  plot_df$description <- str_wrap(plot_df$description, width = 25)
  plot_df$description <- factor(plot_df$description,
                                levels = rev(plot_df$description))
  
  p <- ggplot(
    plot_df,
    aes(x      = fold_enrichment,
        y      = description,
        size   = observed_region_hits,
        colour = p_adjust)
  ) +
    geom_point() +
    scale_colour_gradient(
      low  = "#0072B2",
      high = "#E69F00",
      name = "p.adjust",
      guide = guide_colourbar(
        barwidth       = 0.4,
        barheight      = 3,
        direction      = "vertical",
        title.position = "top",
        title.hjust    = 0.5
      )
    ) +
    scale_size_continuous(
      name  = "Region hits",
      range = c(1, 4)
    ) +
    labs(
      title = paste(label, "-", mode),
      x     = "Fold enrichment",
      y     = NULL
    ) +
    guides(size = guide_legend(override.aes = list(size = 2.5))) +
    theme_classic(base_size = 7) +
    theme(
      axis.text.y  = element_text(size = 6),
      axis.text.x  = element_text(size = 6),
      axis.title.x = element_text(size = 7),
      plot.title   = element_text(size = 7, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 6),
      legend.text  = element_text(size = 5),
      legend.key.height = unit(0.2, "cm"),
      plot.margin  = margin(1, 1, 1, 1)
    )
  
  ggsave(
    filename = file.path(dotplot_dir,
                         paste0("dotplot_rGREAT_", label, "_", mode, ".png")),
    plot   = p,
    width  = 3.5,
    height = 2.8,
    units  = "in",
    dpi    = 300,
    bg     = "white"
  )
  
  cat("Saved dotplot:", label, "-", mode, "\n")
  return(p)
}

# RGREAT FUNCTION
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
  
  peaks_gr <- bulk_peaks_to_gr(top_peaks$feature)
  
  cat(label, "-", mode, "- input peaks:", length(peaks_gr), "\n")
  
  great_res <- great(
    peaks_gr,
    gene_sets  = "GO:BP",
    tss_source = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    background = background_gr,
    cores      = 4
  )
  
  tb     <- getEnrichmentTable(great_res)
  tb_sig <- tb[tb$p_adjust < 0.05, ]
  tb_sig <- tb_sig[order(tb_sig$p_adjust), ]
  
  cat(label, "-", mode, "- significant terms:", nrow(tb_sig), "\n")
  
  write.csv(
    tb_sig[, c("id", "description", "fold_enrichment", "observed_region_hits", "p_adjust")],
    file.path(rgreat_dir, paste0("rGREAT_", label, "_", mode, ".csv")),
    row.names = FALSE
  )
  
  plot_great_dotplot(tb_sig, label, mode)
  
  return(tb_sig)
}

# RUN FOR OPALIN FACTOR 2 BULK VIEW - BOTH DIRECTIONS
cat("\n=== rGREAT: Opalin Factor 2 bulk view ===\n")
run_great_factor(bulk_weights_f2, "bulk", "opalin_factor2_bulk", bulk_consensus, mode = "positive")
run_great_factor(bulk_weights_f2, "bulk", "opalin_factor2_bulk", bulk_consensus, mode = "negative")

cat("\nDone - results saved to:", rgreat_dir, "\n")
cat("Dotplots saved to:", dotplot_dir, "\n")
