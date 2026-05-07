# GENOMIC OVERLAP ANALYSIS FOR BULK VS SC OLIGODENDROCYTE SUBTYPES
# CREATES
# 1. PAIRWISE EULER DIAGRAMS (BULK VS EACH SC SUBTYPE)
# 2. 4-WAY GENOMIC VENN DIAGRAM
# 3. 4-WAY GENOMIC UPSET PLOT
# ALL OVERLAPS ARE BASED ON GenomicRanges COORDINATES, NOT RAW PEAK NAMES

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

out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/plots/new_sc_data_plots_042026"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# COLOUR-BLIND FRIENDLY PALETTE (OKABE-ITO)
col_bulk    <- "#0072B2"
col_opalin  <- "#E69F00"
col_plekhg1 <- "#009E73"
col_opc     <- "#CC79A7"


#FUNCTIONS


# BULK PEAKS: chr1_start_end
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

# SC PEAKS: 1-start-end
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

# BUILD GENOMIC MEMBERSHIP TABLE FROM A UNION OF REGIONS
build_membership_df <- function(union_gr, bulk_gr, opalin_gr, plekhg1_gr, opc_gr) {
  data.frame(
    Bulk    = as.integer(countOverlaps(union_gr, bulk_gr,    ignore.strand = TRUE) > 0),
    Opalin  = as.integer(countOverlaps(union_gr, opalin_gr,  ignore.strand = TRUE) > 0),
    Plekhg1 = as.integer(countOverlaps(union_gr, plekhg1_gr, ignore.strand = TRUE) > 0),
    OPCs    = as.integer(countOverlaps(union_gr, opc_gr,     ignore.strand = TRUE) > 0)
  )
}

# CREATE UNION IDS FOR VENN INPUT
make_union_ids <- function(gr) {
  paste0(as.character(seqnames(gr)), "_", start(gr), "_", end(gr))
}


#1 CREATE GRANGES OBJECTS


bulk_peaks    <- rownames(bulk_counts)
opalin_peaks  <- rownames(opalin_vst)
plekhg1_peaks <- rownames(plekhg1_vst)
opc_peaks     <- rownames(opc_vst)

bulk_gr    <- bulk_peaks_to_gr(bulk_peaks)
opalin_gr  <- sc_peaks_to_gr(opalin_peaks)
plekhg1_gr <- sc_peaks_to_gr(plekhg1_peaks)
opc_gr     <- sc_peaks_to_gr(opc_peaks)

# FIX CHROMOSOME NAMING SO BULK (chr1) AND SC (1) MATCH
seqlevelsStyle(bulk_gr)    <- "UCSC"
seqlevelsStyle(opalin_gr)  <- "UCSC"
seqlevelsStyle(plekhg1_gr) <- "UCSC"
seqlevelsStyle(opc_gr)     <- "UCSC"

# PRINT PEAK COUNTS
cat("Bulk oligodendrocytes:", length(bulk_peaks), "peaks\n")
cat("Opalin+:", length(opalin_peaks), "peaks\n")
cat("Plekhg1+:", length(plekhg1_peaks), "peaks\n")
cat("OPCs:", length(opc_peaks), "peaks\n")

#2 PAIRWISE OVERLAPS: BULK VS EACH SC SUBTYPE


overlaps_opalin  <- suppressWarnings(findOverlaps(bulk_gr, opalin_gr,  ignore.strand = TRUE))
overlaps_plekhg1 <- suppressWarnings(findOverlaps(bulk_gr, plekhg1_gr, ignore.strand = TRUE))
overlaps_opc     <- suppressWarnings(findOverlaps(bulk_gr, opc_gr,     ignore.strand = TRUE))

opalin_shared      <- bulk_peaks[unique(queryHits(overlaps_opalin))]
opalin_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_opalin)]
opalin_unique_sc   <- opalin_peaks[!seq_along(opalin_peaks) %in% unique(subjectHits(overlaps_opalin))]

plekhg1_shared      <- bulk_peaks[unique(queryHits(overlaps_plekhg1))]
plekhg1_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_plekhg1)]
plekhg1_unique_sc   <- plekhg1_peaks[!seq_along(plekhg1_peaks) %in% unique(subjectHits(overlaps_plekhg1))]

opc_shared      <- bulk_peaks[unique(queryHits(overlaps_opc))]
opc_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_opc)]
opc_unique_sc   <- opc_peaks[!seq_along(opc_peaks) %in% unique(subjectHits(overlaps_opc))]

cat("Opalin+  - shared:", length(opalin_shared),
    "| bulk unique:", length(opalin_unique_bulk),
    "| sc unique:", length(opalin_unique_sc), "\n")
cat("Plekhg1+ - shared:", length(plekhg1_shared),
    "| bulk unique:", length(plekhg1_unique_bulk),
    "| sc unique:", length(plekhg1_unique_sc), "\n")
cat("OPCs     - shared:", length(opc_shared),
    "| bulk unique:", length(opc_unique_bulk),
    "| sc unique:", length(opc_unique_sc), "\n")

#3 PAIRWISE EULER DIAGRAMS


fit_opalin <- euler(c(
  "Bulk"         = length(opalin_unique_bulk),
  "Opalin+"      = length(opalin_unique_sc),
  "Bulk&Opalin+" = length(opalin_shared)
))

fit_plekhg1 <- euler(c(
  "Bulk"          = length(plekhg1_unique_bulk),
  "Plekhg1+"      = length(plekhg1_unique_sc),
  "Bulk&Plekhg1+" = length(plekhg1_shared)
))

fit_opc <- euler(c(
  "Bulk"      = length(opc_unique_bulk),
  "OPCs"      = length(opc_unique_sc),
  "Bulk&OPCs" = length(opc_shared)
))

png(file.path(out_dir, "euler_bulk_vs_opalin.png"), width = 800, height = 700, res = 150)
plot(fit_opalin, quantities = TRUE, fills = c(col_bulk, col_opalin), edges = FALSE, main = "Bulk vs Opalin+")
dev.off()

png(file.path(out_dir, "euler_bulk_vs_plekhg1.png"), width = 800, height = 700, res = 150)
plot(fit_plekhg1, quantities = TRUE, fills = c(col_bulk, col_plekhg1), edges = FALSE, main = "Bulk vs Plekhg1+")
dev.off()

png(file.path(out_dir, "euler_bulk_vs_opc.png"), width = 800, height = 700, res = 150)
plot(fit_opc, quantities = TRUE, fills = c(col_bulk, col_opc), edges = FALSE, main = "Bulk vs OPCs")
dev.off()


#4 4-WAY GENOMIC VENN DIAGRAM


# UNION OF ALL REGIONS ACROSS ALL FOUR DATASETS
all_union_gr <- reduce(c(bulk_gr, opalin_gr, plekhg1_gr, opc_gr), ignore.strand = TRUE)
cat("Total union genomic regions across all datasets:", length(all_union_gr), "\n")

union_ids <- make_union_ids(all_union_gr)

peak_sets_peaks <- list(
  Bulk    = bulk_peaks,
  Opalin  = opalin_peaks[unique(subjectHits(findOverlaps(bulk_gr, opalin_gr,  ignore.strand = TRUE)))],
  Plekhg1 = plekhg1_peaks[unique(subjectHits(findOverlaps(bulk_gr, plekhg1_gr, ignore.strand = TRUE)))],
  OPCs    = opc_peaks[unique(subjectHits(findOverlaps(bulk_gr, opc_gr,     ignore.strand = TRUE)))]
)

p_venn <- ggVennDiagram(
  peak_sets_peaks,
  label_alpha = 0,
  set_color   = c(col_bulk, col_opalin, col_plekhg1, col_opc)
) +
  scale_fill_gradient(low = "white", high = col_bulk) +
  labs(title = "Genomic peak overlap across bulk and SC oligodendrocyte subtypes") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

print(p_venn)

ggsave(
  file.path(out_dir, "venn_all_datasets_genomic.png"),
  plot = p_venn,
  width = 10,
  height = 8,
  dpi = 150
)


#5 4-WAY GENOMIC UPSET PLOT


upset_df <- build_membership_df(all_union_gr, bulk_gr, opalin_gr, plekhg1_gr, opc_gr)

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


#6 OPTIONAL (extra): PAIRWISE GENOMIC VENN FOR BULK VS OPALIN


bulk_opalin_union <- reduce(c(bulk_gr, opalin_gr), ignore.strand = TRUE)
bulk_opalin_ids <- make_union_ids(bulk_opalin_union)

pair_sets_opalin <- list(
  Bulk     = bulk_opalin_ids[countOverlaps(bulk_opalin_union, bulk_gr,   ignore.strand = TRUE) > 0],
  "Opalin+" = bulk_opalin_ids[countOverlaps(bulk_opalin_union, opalin_gr, ignore.strand = TRUE) > 0]
)

p_opalin <- ggVennDiagram(
  pair_sets_opalin,
  label_alpha = 0,
  set_color = c(col_bulk, col_opalin)
) +
  scale_fill_gradient(low = "white", high = col_opalin) +
  labs(title = "Bulk vs Opalin+ (genomic overlap)")

print(p_opalin)

ggsave(
  file.path(out_dir, "venn_bulk_vs_opalin_genomic.png"),
  plot = p_opalin,
  width = 8,
  height = 6,
  dpi = 150
)


# SAVE OBJECTS FOR DOWNSTREAM USE


save(
  bulk_counts, oligo_bulk_data, bulk_meta, sample_mapping,
  opalin_vst, plekhg1_vst, opc_vst,
  bulk_peaks, opalin_peaks, plekhg1_peaks, opc_peaks,
  bulk_gr, opalin_gr, plekhg1_gr, opc_gr,
  overlaps_opalin, overlaps_plekhg1, overlaps_opc,
  opalin_shared, opalin_unique_bulk, opalin_unique_sc,
  plekhg1_shared, plekhg1_unique_bulk, plekhg1_unique_sc,
  opc_shared, opc_unique_bulk, opc_unique_sc,
  all_union_gr, upset_df, peak_sets_genomic,
  file = "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_granges_vst.RData"
)

# QUICK CHECKS
dim(opalin_vst)
head(rownames(opalin_vst))



#confimring data

cat("Opalin+  - shared:", length(opalin_shared),
    "| bulk unique:", length(opalin_unique_bulk),
    "| sc unique:", length(opalin_unique_sc), "\n")
cat("Plekhg1+ - shared:", length(plekhg1_shared),
    "| bulk unique:", length(plekhg1_unique_bulk),
    "| sc unique:", length(plekhg1_unique_sc), "\n")
cat("OPCs     - shared:", length(opc_shared),
    "| bulk unique:", length(opc_unique_bulk),
    "| sc unique:", length(opc_unique_sc), "\n")
