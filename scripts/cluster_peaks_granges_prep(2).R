
# GENOMIC OVERLAP ANALYSIS: BULK VS SC OLIGODENDROCYTE SUBTYPES
# OUTPUTS:
#   1. Pairwise Euler diagrams (Bulk vs each SC subtype)
#   2. 4-way genomic Venn diagram
#   3. 4-way genomic UpSet plot
# All overlaps are based on GenomicRanges coordinates, not raw peak names


.libPaths(c(
  "/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5",
  .libPaths()
))

load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_pseudobulked_vst.RData")

library(GenomicRanges)
library(GenomeInfoDb)
library(IRanges)
library(eulerr)
library(ggplot2)
library(ggVennDiagram)
library(UpSetR)

out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/plots/new_sc_data_plots_042026/descriptive_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Okabe-Ito colour-blind friendly palette
col_bulk    <- "#0072B2"
col_opalin  <- "#E69F00"
col_plekhg1 <- "#009E73"
col_opc     <- "#CC79A7"



# HELPER FUNCTIONS
# parsing 

# Parse bulk peak names: chr1_start_end
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

# Parse SC peak names: chr-start-end
sc_peaks_to_gr <- function(peaks) {
  split_peaks <- strsplit(as.character(peaks), "-")
  seqs   <- sapply(split_peaks, `[`, 1)
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

# Create string IDs from a GRanges object (used as set membership keys)
make_union_ids <- function(gr) {
  paste0(as.character(seqnames(gr)), "_", start(gr), "_", end(gr))
}


# SECTION 1: BUILD GRanges OBJECTS


bulk_peaks    <- rownames(bulk_counts)
opalin_peaks  <- rownames(opalin_vst)
plekhg1_peaks <- rownames(plekhg1_vst)
opc_peaks     <- rownames(opc_vst)

bulk_gr    <- bulk_peaks_to_gr(bulk_peaks)
opalin_gr  <- sc_peaks_to_gr(opalin_peaks)
plekhg1_gr <- sc_peaks_to_gr(plekhg1_peaks)
opc_gr     <- sc_peaks_to_gr(opc_peaks)

# Harmonise chromosome naming (bulk = chr1, SC = 1 -> both to UCSC style)
seqlevelsStyle(bulk_gr)    <- "UCSC"
seqlevelsStyle(opalin_gr)  <- "UCSC"
seqlevelsStyle(plekhg1_gr) <- "UCSC"
seqlevelsStyle(opc_gr)     <- "UCSC"

cat("Bulk oligodendrocytes:", length(bulk_peaks), "peaks\n")
cat("Opalin+:              ", length(opalin_peaks), "peaks\n")
cat("Plekhg1+:             ", length(plekhg1_peaks), "peaks\n")
cat("OPCs:                 ", length(opc_peaks), "peaks\n")


# =============================================================================
# SECTION 2: PAIRWISE OVERLAPS (BULK VS EACH SC SUBTYPE)
# =============================================================================

overlaps_opalin  <- suppressWarnings(findOverlaps(bulk_gr, opalin_gr,  ignore.strand = TRUE))
overlaps_plekhg1 <- suppressWarnings(findOverlaps(bulk_gr, plekhg1_gr, ignore.strand = TRUE))
overlaps_opc     <- suppressWarnings(findOverlaps(bulk_gr, opc_gr,     ignore.strand = TRUE))

# Opalin+
opalin_shared      <- length(unique(queryHits(overlaps_opalin)))
opalin_unique_bulk <- length(bulk_peaks) - opalin_shared
opalin_unique_sc   <- length(opalin_peaks) - length(unique(subjectHits(overlaps_opalin)))

# Plekhg1+
plekhg1_shared      <- length(unique(queryHits(overlaps_plekhg1)))
plekhg1_unique_bulk <- length(bulk_peaks) - plekhg1_shared
plekhg1_unique_sc   <- length(plekhg1_peaks) - length(unique(subjectHits(overlaps_plekhg1)))

# OPCs
opc_shared      <- length(unique(queryHits(overlaps_opc)))
opc_unique_bulk <- length(bulk_peaks) - opc_shared
opc_unique_sc   <- length(opc_peaks) - length(unique(subjectHits(overlaps_opc)))

cat("\nOpalin+  - shared:", opalin_shared,
    "| bulk unique:", opalin_unique_bulk,
    "| SC unique:", opalin_unique_sc, "\n")
cat("Plekhg1+ - shared:", plekhg1_shared,
    "| bulk unique:", plekhg1_unique_bulk,
    "| SC unique:", plekhg1_unique_sc, "\n")
cat("OPCs     - shared:", opc_shared,
    "| bulk unique:", opc_unique_bulk,
    "| SC unique:", opc_unique_sc, "\n")


# =============================================================================
# SECTION 3: PAIRWISE EULER DIAGRAMS
# =============================================================================

fit_opalin <- euler(c(
  "Bulk"         = opalin_unique_bulk,
  "Opalin+"      = opalin_unique_sc,
  "Bulk&Opalin+" = opalin_shared
))

fit_plekhg1 <- euler(c(
  "Bulk"          = plekhg1_unique_bulk,
  "Plekhg1+"      = plekhg1_unique_sc,
  "Bulk&Plekhg1+" = plekhg1_shared
))

fit_opc <- euler(c(
  "Bulk"      = opc_unique_bulk,
  "OPCs"      = opc_unique_sc,
  "Bulk&OPCs" = opc_shared
))

png(file.path(out_dir, "euler_bulk_vs_opalin.png"), width = 800, height = 700, res = 150)
plot(fit_opalin,
     quantities = TRUE,
     fills      = c(col_bulk, col_opalin),
     edges      = FALSE,
     main       = "Bulk vs Opalin+")
dev.off()

png(file.path(out_dir, "euler_bulk_vs_plekhg1.png"), width = 800, height = 700, res = 150)
plot(fit_plekhg1,
     quantities = TRUE,
     fills      = c(col_bulk, col_plekhg1),
     edges      = FALSE,
     main       = "Bulk vs Plekhg1+")
dev.off()

png(file.path(out_dir, "euler_bulk_vs_opc.png"), width = 800, height = 700, res = 150)
plot(fit_opc,
     quantities = TRUE,
     fills      = c(col_bulk, col_opc),
     edges      = FALSE,
     main       = "Bulk vs OPCs")
dev.off()

cat("\nPairwise Euler diagrams saved.\n")



# SECTION 4: 4-WAY GENOMIC VENN DIAGRAM
# -----------------------------------------------------------------------------
# Key fix: build set membership from the shared reduce() union, using
# countOverlaps() for each view — exactly the same logic as the UpSet plot.
# This ensures Venn and UpSet counts are consistent.


# Reduce all four peak sets into a single non-overlapping union of regions
all_union_gr <- reduce(
  c(bulk_gr, opalin_gr, plekhg1_gr, opc_gr),
  ignore.strand = TRUE
)
cat("\nTotal union genomic regions across all four datasets:", length(all_union_gr), "\n")

# Assign a unique string ID to each union region
union_ids <- make_union_ids(all_union_gr)

# For each view, collect the union IDs of regions that overlap that view's peaks
peak_sets_genomic <- list(
  Bulk    = union_ids[countOverlaps(all_union_gr, bulk_gr,    ignore.strand = TRUE) > 0],
  Opalin  = union_ids[countOverlaps(all_union_gr, opalin_gr,  ignore.strand = TRUE) > 0],
  Plekhg1 = union_ids[countOverlaps(all_union_gr, plekhg1_gr, ignore.strand = TRUE) > 0],
  OPCs    = union_ids[countOverlaps(all_union_gr, opc_gr,     ignore.strand = TRUE) > 0]
)

cat("Bulk union regions:    ", length(peak_sets_genomic$Bulk), "\n")
cat("Opalin union regions:  ", length(peak_sets_genomic$Opalin), "\n")
cat("Plekhg1 union regions: ", length(peak_sets_genomic$Plekhg1), "\n")
cat("OPCs union regions:    ", length(peak_sets_genomic$OPCs), "\n")

p_venn <- ggVennDiagram(
  peak_sets_genomic,
  label_alpha = 0,
  set_color   = c(col_bulk, col_opalin, col_plekhg1, col_opc)
) +
  scale_fill_gradient(low = "white", high = col_bulk) +
  labs(title = "Peak overlap across bulk and SC oligodendrocyte subtypes") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave(
  file.path(out_dir, "venn_all_datasets_genomic.png"),
  plot   = p_venn,
  width  = 10,
  height = 8,
  dpi    = 150
)

cat("4-way Venn diagram saved.\n")


# =============================================================================
# SECTION 5: 4-WAY GENOMIC UPSET PLOT
# -----------------------------------------------------------------------------
# Uses the same all_union_gr and countOverlaps() logic as the Venn above,
# so counts will match between the two figures.
# =============================================================================

upset_df <- data.frame(
  Bulk    = as.integer(countOverlaps(all_union_gr, bulk_gr,    ignore.strand = TRUE) > 0),
  Opalin  = as.integer(countOverlaps(all_union_gr, opalin_gr,  ignore.strand = TRUE) > 0),
  Plekhg1 = as.integer(countOverlaps(all_union_gr, plekhg1_gr, ignore.strand = TRUE) > 0),
  OPCs    = as.integer(countOverlaps(all_union_gr, opc_gr,     ignore.strand = TRUE) > 0)
)

png(file.path(out_dir, "upset_all_datasets_genomic.png"), width = 1200, height = 900, res = 150)
upset(
  upset_df,
  sets            = c("Bulk", "Opalin", "Plekhg1", "OPCs"),
  order.by        = "freq",
  sets.bar.color  = c(col_bulk, col_opalin, col_plekhg1, col_opc),
  text.scale      = 1.3,
  mainbar.y.label = "Number of genomic regions",
  sets.x.label    = "Total genomic regions per dataset"
)
dev.off()

cat("4-way UpSet plot saved.\n")



# SECTION 6: SAVE ALL OBJECTS FOR DOWNSTREAM USE


save(
  bulk_counts, oligo_bulk_data, bulk_meta, sample_mapping,
  opalin_vst, plekhg1_vst, opc_vst,
  bulk_peaks, opalin_peaks, plekhg1_peaks, opc_peaks,
  bulk_gr, opalin_gr, plekhg1_gr, opc_gr,
  overlaps_opalin, overlaps_plekhg1, overlaps_opc,
  opalin_shared, opalin_unique_bulk, opalin_unique_sc,
  plekhg1_shared, plekhg1_unique_bulk, plekhg1_unique_sc,
  opc_shared, opc_unique_bulk, opc_unique_sc,
  all_union_gr, union_ids, peak_sets_genomic, upset_df,
  file = "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_granges_vst.RData"
)

cat("\nAll objects saved to .RData file.\n")
cat("Done.\n")
