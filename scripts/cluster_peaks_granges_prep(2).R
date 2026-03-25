#make granges objects to analyse peaks overlaps between subtypes

.libPaths(c("/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))

load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_pseudobulked_vst.RData")

library(GenomicRanges)
library(GenomeInfoDb)
#create Grnages objects

# bulk peaks (chr1_start_end format)
bulk_peaks <- rownames(bulk_counts)
bulk_gr <- GRanges(
  seqnames = gsub("_.*", "", bulk_peaks),
  ranges   = IRanges(
    start = as.numeric(gsub(".*_(\\d+)_\\d+$", "\\1", bulk_peaks)),
    end   = as.numeric(gsub(".*_(\\d+)$", "\\1", bulk_peaks))
  )
)

# function for SC peaks (1-start-end format)
make_gr <- function(peaks) {
  GRanges(
    seqnames = sapply(strsplit(peaks, "-"), `[`, 1),
    ranges   = IRanges(
      start = as.numeric(sapply(strsplit(peaks, "-"), `[`, 2)),
      end   = as.numeric(sapply(strsplit(peaks, "-"), `[`, 3))
    )
  )
}

opalin_peaks  <- rownames(opalin_vst)
plekhg1_peaks <- rownames(plekhg1_vst)
opc_peaks     <- rownames(opc_vst)

opalin_gr  <- make_gr(opalin_peaks)
plekhg1_gr <- make_gr(plekhg1_peaks)
opc_gr     <- make_gr(opc_peaks)

# fix chromosome naming
seqlevelsStyle(bulk_gr)    <- "UCSC"
seqlevelsStyle(opalin_gr)  <- "UCSC"
seqlevelsStyle(plekhg1_gr) <- "UCSC"
seqlevelsStyle(opc_gr)     <- "UCSC"

#find overlaps
overlaps_opalin  <- suppressWarnings(findOverlaps(bulk_gr, opalin_gr,  ignore.strand = TRUE))
overlaps_plekhg1 <- suppressWarnings(findOverlaps(bulk_gr, plekhg1_gr, ignore.strand = TRUE))
overlaps_opc     <- suppressWarnings(findOverlaps(bulk_gr, opc_gr,     ignore.strand = TRUE))

# classify peaks
opalin_shared      <- bulk_peaks[unique(queryHits(overlaps_opalin))]
opalin_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_opalin)]
opalin_unique_sc   <- opalin_peaks[!seq_along(opalin_peaks) %in% subjectHits(overlaps_opalin)]

plekhg1_shared      <- bulk_peaks[unique(queryHits(overlaps_plekhg1))]
plekhg1_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_plekhg1)]
plekhg1_unique_sc   <- plekhg1_peaks[!seq_along(plekhg1_peaks) %in% subjectHits(overlaps_plekhg1)]

opc_shared      <- bulk_peaks[unique(queryHits(overlaps_opc))]
opc_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_opc)]
opc_unique_sc   <- opc_peaks[!seq_along(opc_peaks) %in% subjectHits(overlaps_opc)]

# summarise
cat("Opalin+  - shared:", length(opalin_shared),
    "| bulk unique:", length(opalin_unique_bulk),
    "| sc unique:", length(opalin_unique_sc), "\n")
cat("Plekhg1+ - shared:", length(plekhg1_shared),
    "| bulk unique:", length(plekhg1_unique_bulk),
    "| sc unique:", length(plekhg1_unique_sc), "\n")
cat("OPCs     - shared:", length(opc_shared),
    "| bulk unique:", length(opc_unique_bulk),
    "| sc unique:", length(opc_unique_sc), "\n")

#save in saved_objects
save(bulk_counts, oligo_bulk_data, bulk_meta, sample_mapping,
     opalin_vst, plekhg1_vst, opc_vst,
     opalin_pseudo, plekhg1_pseudo, opc_pseudo,
     bulk_peaks, opalin_peaks, plekhg1_peaks, opc_peaks,
     bulk_gr, opalin_gr, plekhg1_gr, opc_gr,
     overlaps_opalin, overlaps_plekhg1, overlaps_opc,
     opalin_shared, opalin_unique_bulk, opalin_unique_sc,
     plekhg1_shared, plekhg1_unique_bulk, plekhg1_unique_sc,
     opc_shared, opc_unique_bulk, opc_unique_sc,
     file = "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_granges_vst.RData")


