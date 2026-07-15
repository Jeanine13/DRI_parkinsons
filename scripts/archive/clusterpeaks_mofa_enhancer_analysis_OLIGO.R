#In this analysis, were looking at how many features in each factor overlap with known enhancers for oligodendrocytes, to get some info on whether single-cell data is more enriched for cell type specific peaks


# load oligo enhancers from excel spreadsheet ( notts 2019 paper)
library(readxl)
enhancers <- read_excel("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/aay0793-nott-table-s5.xlsx")

library(readxl)
# see what sheets are available
excel_sheets("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/aay0793-nott-table-s5.xlsx")
enhancers <- read_excel("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/aay0793-nott-table-s5.xlsx",
                        sheet = "Oligo enhancers")

head(enhancers)
colnames(enhancers)
dim(enhancers)

#clean up the file rows and create granges object
# clean up the data frame
enhancers <- enhancers[-1, ]  # remove first NA row
colnames(enhancers) <- c("chr", "start", "end")

# convert to correct types
enhancers$start <- as.numeric(enhancers$start)
enhancers$end   <- as.numeric(enhancers$end)

# remove any remaining NA rows
enhancers <- enhancers[!is.na(enhancers$start), ]

head(enhancers)
dim(enhancers)

# create GRanges object
enhancers_gr <- GRanges(
  seqnames = enhancers$chr,
  ranges   = IRanges(start = enhancers$start, end = enhancers$end)
)
#granges object
enhancers_gr

head(seqlevels(enhancers_gr))
head(seqlevels(bulk_gr))
head(seqlevels(opalin_gr))

load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_granges_vst.RData")

# for SC peaks (1-start-end format)
peak_to_granges <- function(peaks) {
  parts <- do.call(rbind, strsplit(peaks, "-"))
  gr <- GRanges(
    seqnames = paste0("chr", parts[, 1]),
    ranges   = IRanges(
      start = as.numeric(parts[, 2]),
      end   = as.numeric(parts[, 3])
    )
  )
  names(gr) <- peaks
  gr
}

# for bulk peaks (chr1_start_end format)
bulk_peak_to_granges <- function(peaks) {
  gr <- GRanges(
    seqnames = gsub("_.*", "", peaks),
    ranges   = IRanges(
      start = as.numeric(gsub(".*_(\\d+)_\\d+$", "\\1", peaks)),
      end   = as.numeric(gsub(".*_(\\d+)$", "\\1", peaks))
    )
  )
  names(gr) <- peaks
  gr
}

#baxkground grnages
bulk_vars    <- apply(bulk_counts, 1, var)
bulk_top_idx <- head(order(bulk_vars, decreasing = TRUE), 5000)
bulk_hvf     <- bulk_counts[bulk_top_idx, ]

opalin_vars    <- apply(opalin_vst, 1, var)
opalin_top_idx <- head(order(opalin_vars, decreasing = TRUE), 5000)
opalin_hvf     <- opalin_vst[opalin_top_idx, ]

bulk_background_gr   <- bulk_peak_to_granges(rownames(bulk_hvf))
opalin_background_gr <- peak_to_granges(rownames(opalin_hvf))

#get mofa weights
mofa_opalin_vst <- load_model("mofa_opalin_vst.hdf5")

weights_df <- get_weights(mofa_opalin_vst,
                         views   = "all",
                         factors = "all",
                         as.data.frame = TRUE)
weights_df$feature <- as.character(weights_df$feature)
weights_df$factor  <- as.character(weights_df$factor)
weights_df$view    <- as.character(weights_df$view)

#functions

get_top_features <- function(df, factor_name, view_name, n = 200, direction = "positive") {
  sub <- df[df$factor == factor_name & df$view == view_name, ]
  
  if (direction == "positive") {
    sub <- sub[sub$value > 0, ]
    sub <- sub[order(sub$value, decreasing = TRUE), ]
  } else if (direction == "negative") {
    sub <- sub[sub$value < 0, ]
    sub <- sub[order(sub$value, decreasing = FALSE), ]
  } else {
    sub <- sub[order(abs(sub$value), decreasing = TRUE), ]
  }
  
  head(sub$feature, n)
}

test_enhancer_enrichment <- function(factor_peaks_gr, all_peaks_gr, enhancers_gr, label) {
  
  other_peaks_gr <- all_peaks_gr[!names(all_peaks_gr) %in% names(factor_peaks_gr)]
  
  factor_hits            <- suppressWarnings(findOverlaps(factor_peaks_gr, enhancers_gr, ignore.strand = TRUE))
  factor_enhancer_overlap <- length(unique(queryHits(factor_hits)))
  factor_no_overlap       <- length(factor_peaks_gr) - factor_enhancer_overlap
  
  bg_hits            <- suppressWarnings(findOverlaps(other_peaks_gr, enhancers_gr, ignore.strand = TRUE))
  bg_enhancer_overlap <- length(unique(queryHits(bg_hits)))
  bg_no_overlap       <- length(other_peaks_gr) - bg_enhancer_overlap
  
  mat <- matrix(
    c(factor_enhancer_overlap, factor_no_overlap,
      bg_enhancer_overlap,     bg_no_overlap),
    nrow = 2, byrow = TRUE
  )
  
  fisher <- fisher.test(mat, alternative = "greater")
  
  data.frame(
    label             = label,
    factor_overlap    = factor_enhancer_overlap,
    factor_total      = length(factor_peaks_gr),
    background_overlap = bg_enhancer_overlap,
    background_total  = length(other_peaks_gr),
    factor_prop       = factor_enhancer_overlap / length(factor_peaks_gr),
    background_prop   = bg_enhancer_overlap / length(other_peaks_gr),
    odds_ratio        = as.numeric(fisher$estimate),
    p_value           = fisher$p.value
  )
}

#enrichment analysis

all_results        <- list()
factors_to_test    <- unique(weights_df$factor)
directions_to_test <- c("positive", "negative")

for (view_name in c("bulk", "opalin")) {
  
  bg_gr <- if (view_name == "bulk") bulk_background_gr else opalin_background_gr
  
  for (factor_name in factors_to_test) {
    for (dir_name in directions_to_test) {
      
      top_peaks <- get_top_features(weights_df, factor_name, view_name, n = 200, direction = dir_name)
      if (length(top_peaks) == 0) next
      
      factor_gr <- if (view_name == "bulk") bulk_peak_to_granges(top_peaks) else peak_to_granges(top_peaks)
      
      res <- test_enhancer_enrichment(factor_gr, bg_gr, enhancers_gr,
                                      paste(view_name, factor_name, dir_name, sep = "_"))
      
      res$view      <- view_name
      res$factor    <- factor_name
      res$direction <- dir_name
      
      all_results[[length(all_results) + 1]] <- res
    }
  }
}

results_df       <- bind_rows(all_results)
results_df$padj  <- p.adjust(results_df$p_value, method = "BH")

print(results_df)

write.csv(results_df,
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/mofa_oligo_enhancer_enrichment_results.csv",
          row.names = FALSE)




