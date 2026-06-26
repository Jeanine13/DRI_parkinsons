# GENOMIC OVERLAP ANALYSIS: BULK VS SC OLIGODENDROCYTE SUBTYPES
#
# Outputs:
# 1. Pairwise Euler diagrams
# 2. 4-way genomic Venn diagram
# 3. 4-way genomic UpSet plot
#
# Overlaps are based on genomic coordinates using GenomicRanges


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


out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/plots/new_sc_data_plots_042026/descriptive_analysis/da_2026"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


col_bulk <- "#0072B2"
col_opalin <- "#E69F00"
col_plekhg1 <- "#009E73"
col_opc <- "#CC79A7"



# Functions

bulk_peaks_to_gr <- function(peaks){
  
  starts <- as.numeric(
    gsub(".*_(\\d+)_\\d+$", "\\1", peaks)
  )
  
  ends <- as.numeric(
    gsub(".*_(\\d+)$", "\\1", peaks)
  )
  
  
  GRanges(
    seqnames = gsub("_.*", "", peaks),
    ranges = IRanges(
      start = starts,
      end = ends
    )
  )
}



sc_peaks_to_gr <- function(peaks){
  
  split_peaks <- strsplit(
    as.character(peaks),
    "-"
  )
  
  GRanges(
    seqnames = sapply(split_peaks, `[`, 1),
    ranges = IRanges(
      start = as.numeric(sapply(split_peaks, `[`, 2)),
      end = as.numeric(sapply(split_peaks, `[`, 3))
    )
  )
}



make_union_ids <- function(gr){
  
  paste0(
    as.character(seqnames(gr)),
    "_",
    start(gr),
    "_",
    end(gr)
  )
}



# Create GRanges objects

bulk_peaks <- rownames(bulk_counts)

opalin_peaks <- rownames(opalin_vst)

plekhg1_peaks <- rownames(plekhg1_vst)

opc_peaks <- rownames(opc_vst)


bulk_gr <- bulk_peaks_to_gr(bulk_peaks)

opalin_gr <- sc_peaks_to_gr(opalin_peaks)

plekhg1_gr <- sc_peaks_to_gr(plekhg1_peaks)

opc_gr <- sc_peaks_to_gr(opc_peaks)



# Harmonise chromosome naming

seqlevelsStyle(bulk_gr) <- "UCSC"
seqlevelsStyle(opalin_gr) <- "UCSC"
seqlevelsStyle(plekhg1_gr) <- "UCSC"
seqlevelsStyle(opc_gr) <- "UCSC"



cat("\nOriginal peak numbers\n")

bulk_n <- length(bulk_gr)
opalin_n <- length(opalin_gr)
plekhg1_n <- length(plekhg1_gr)
opc_n <- length(opc_gr)


cat("Bulk:", bulk_n, "\n")
cat("Opalin:", opalin_n, "\n")
cat("Plekhg1:", plekhg1_n, "\n")
cat("OPC:", opc_n, "\n")



# Pairwise overlaps

overlaps_opalin <- findOverlaps(
  bulk_gr,
  opalin_gr,
  ignore.strand = TRUE
)


overlaps_plekhg1 <- findOverlaps(
  bulk_gr,
  plekhg1_gr,
  ignore.strand = TRUE
)


overlaps_opc <- findOverlaps(
  bulk_gr,
  opc_gr,
  ignore.strand = TRUE
)



# Shared and unique peaks

opalin_shared <- length(unique(queryHits(overlaps_opalin)))

opalin_unique_bulk <- bulk_n - opalin_shared

opalin_unique_sc <- opalin_n -
  length(unique(subjectHits(overlaps_opalin)))



plekhg1_shared <- length(unique(queryHits(overlaps_plekhg1)))

plekhg1_unique_bulk <- bulk_n - plekhg1_shared

plekhg1_unique_sc <- plekhg1_n -
  length(unique(subjectHits(overlaps_plekhg1)))



opc_shared <- length(unique(queryHits(overlaps_opc)))

opc_unique_bulk <- bulk_n - opc_shared

opc_unique_sc <- opc_n -
  length(unique(subjectHits(overlaps_opc)))



cat("\nOverlap summary\n")

cat(
  "Opalin shared:",
  opalin_shared,
  "| bulk unique:",
  opalin_unique_bulk,
  "| SC unique:",
  opalin_unique_sc,
  "\n"
)


cat(
  "Plekhg1 shared:",
  plekhg1_shared,
  "| bulk unique:",
  plekhg1_unique_bulk,
  "| SC unique:",
  plekhg1_unique_sc,
  "\n"
)


cat(
  "OPC shared:",
  opc_shared,
  "| bulk unique:",
  opc_unique_bulk,
  "| SC unique:",
  opc_unique_sc,
  "\n"
)



# Euler diagrams

make_euler <- function(shared, bulk_unique, sc_unique, label, colour){
  
  fit <- euler(c(
    Bulk = bulk_unique,
    label = sc_unique,
    Shared = shared
  ))
  
  
  png(
    file.path(
      out_dir,
      paste0("euler_bulk_vs_", label, ".png")
    ),
    width = 800,
    height = 700,
    res = 150
  )
  
  
  plot(
    fit,
    quantities = TRUE,
    fills = c(col_bulk, colour),
    edges = FALSE,
    main = paste("Bulk vs", label)
  )
  
  
  dev.off()
  
}



make_euler(
  opalin_shared,
  opalin_unique_bulk,
  opalin_unique_sc,
  "Opalin",
  col_opalin
)


make_euler(
  plekhg1_shared,
  plekhg1_unique_bulk,
  plekhg1_unique_sc,
  "Plekhg1",
  col_plekhg1
)


make_euler(
  opc_shared,
  opc_unique_bulk,
  opc_unique_sc,
  "OPC",
  col_opc
)



# Genomic union for Venn and UpSet

all_union_gr <- reduce(
  c(
    bulk_gr,
    opalin_gr,
    plekhg1_gr,
    opc_gr
  ),
  ignore.strand = TRUE
)


union_ids <- make_union_ids(all_union_gr)



peak_sets_genomic <- list(
  
  Bulk =
    union_ids[
      countOverlaps(all_union_gr, bulk_gr) > 0
    ],
  
  Opalin =
    union_ids[
      countOverlaps(all_union_gr, opalin_gr) > 0
    ],
  
  Plekhg1 =
    union_ids[
      countOverlaps(all_union_gr, plekhg1_gr) > 0
    ],
  
  OPCs =
    union_ids[
      countOverlaps(all_union_gr, opc_gr) > 0
    ]
  
)



cat("\nTotal genomic regions:",
    length(all_union_gr),
    "\n")



# Venn diagram

p_venn <- ggVennDiagram(
  peak_sets_genomic,
  label_alpha = 0
)


ggsave(
  file.path(
    out_dir,
    "venn_all_datasets_genomic.png"
  ),
  p_venn,
  width = 10,
  height = 8,
  dpi = 150
)



# UpSet plot

upset_df <- data.frame(
  
  Bulk =
    as.integer(countOverlaps(all_union_gr, bulk_gr) > 0),
  
  Opalin =
    as.integer(countOverlaps(all_union_gr, opalin_gr) > 0),
  
  Plekhg1 =
    as.integer(countOverlaps(all_union_gr, plekhg1_gr) > 0),
  
  OPCs =
    as.integer(countOverlaps(all_union_gr, opc_gr) > 0)
  
)



png(
  file.path(
    out_dir,
    "upset_all_datasets_genomic.png"
  ),
  width = 1200,
  height = 900,
  res = 150
)


upset(
  upset_df,
  sets = c(
    "Bulk",
    "Opalin",
    "Plekhg1",
    "OPCs"
  ),
  order.by = "freq",
  mainbar.y.label = "Number of genomic regions",
  sets.x.label = "Total genomic regions"
)


dev.off()



# Save objects

save(
  bulk_gr,
  opalin_gr,
  plekhg1_gr,
  opc_gr,
  
  bulk_peaks,
  opalin_peaks,
  plekhg1_peaks,
  opc_peaks,
  
  overlaps_opalin,
  overlaps_plekhg1,
  overlaps_opc,
  
  opalin_shared,
  opalin_unique_bulk,
  opalin_unique_sc,
  
  plekhg1_shared,
  plekhg1_unique_bulk,
  plekhg1_unique_sc,
  
  opc_shared,
  opc_unique_bulk,
  opc_unique_sc,
  
  all_union_gr,
  peak_sets_genomic,
  upset_df,
  
  file =
    "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/genomic_overlap_analysis.RData"
)


cat("\nDone\n")
