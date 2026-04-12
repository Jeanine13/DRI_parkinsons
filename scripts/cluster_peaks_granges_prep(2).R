#make granges objects to analyse peaks overlaps between subtypes
#CREATE 4 WAY VENN DIAGRAM and seperate diagrams for each oligo subtype

.libPaths(c("/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))

load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_pseudobulked_vst.RData")

#packages for descriptive analyses and plots
library(GenomicRanges)
library(GenomeInfoDb)
library(eulerr)
library(ggplot2)
library(ggVennDiagram)

out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/plots/new_sc_data_plots_042026"

# colour-blind friendly palette (Okabe-Ito)
col_bulk    <- "#0072B2"
col_opalin  <- "#E69F00"
col_plekhg1 <- "#009E73"
col_opc     <- "#CC79A7"

# CREATE GRANGES OBJECTS FROM PEAK ROWNAMES
# BULK PEAKS ARE IN chr1_start_end FORMAT
bulk_peaks <- rownames(bulk_counts)
bulk_gr <- GRanges(
  seqnames = gsub("_.*", "", bulk_peaks),
  ranges   = IRanges(
    start = as.numeric(gsub(".*_(\\d+)_\\d+$", "\\1", bulk_peaks)),
    end   = as.numeric(gsub(".*_(\\d+)$", "\\1", bulk_peaks))
  )
)

# SC PEAKS ARE IN 1-start-end FORMAT SO NEED A DIFFERENT FUNCTION
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

# FIX CHROMOSOME NAMING SO BULK (chr1) AND SC (1) MATCH
seqlevelsStyle(bulk_gr)    <- "UCSC"
seqlevelsStyle(opalin_gr)  <- "UCSC"
seqlevelsStyle(plekhg1_gr) <- "UCSC"
seqlevelsStyle(opc_gr)     <- "UCSC"

# PRINT PEAK COUNTS FOR EACH DATASET
cat("Bulk oligodendrocytes:", length(bulk_peaks), "peaks\n")
cat("Opalin+:", length(opalin_peaks), "peaks\n")
cat("Plekhg1+:", length(plekhg1_peaks), "peaks\n")
cat("OPCs:", length(opc_peaks), "peaks\n")

# FIND COORDINATE-BASED OVERLAPS BETWEEN BULK AND EACH SC SUBTYPE
# ignore.strand = TRUE BECAUSE ATAC/CUT&TAG PEAKS ARE NOT STRAND SPECIFIC
overlaps_opalin  <- suppressWarnings(findOverlaps(bulk_gr, opalin_gr,  ignore.strand = TRUE))
overlaps_plekhg1 <- suppressWarnings(findOverlaps(bulk_gr, plekhg1_gr, ignore.strand = TRUE))
overlaps_opc     <- suppressWarnings(findOverlaps(bulk_gr, opc_gr,     ignore.strand = TRUE))

# CLASSIFY PEAKS AS SHARED OR UNIQUE TO BULK OR SC
opalin_shared      <- bulk_peaks[unique(queryHits(overlaps_opalin))]
opalin_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_opalin)]
opalin_unique_sc   <- opalin_peaks[!seq_along(opalin_peaks) %in% subjectHits(overlaps_opalin)]

plekhg1_shared      <- bulk_peaks[unique(queryHits(overlaps_plekhg1))]
plekhg1_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_plekhg1)]
plekhg1_unique_sc   <- plekhg1_peaks[!seq_along(plekhg1_peaks) %in% subjectHits(overlaps_plekhg1)]

opc_shared      <- bulk_peaks[unique(queryHits(overlaps_opc))]
opc_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_opc)]
opc_unique_sc   <- opc_peaks[!seq_along(opc_peaks) %in% subjectHits(overlaps_opc)]

# PRINT OVERLAP SUMMARY
cat("Opalin+  - shared:", length(opalin_shared),
    "| bulk unique:", length(opalin_unique_bulk),
    "| sc unique:", length(opalin_unique_sc), "\n")
cat("Plekhg1+ - shared:", length(plekhg1_shared),
    "| bulk unique:", length(plekhg1_unique_bulk),
    "| sc unique:", length(plekhg1_unique_sc), "\n")
cat("OPCs     - shared:", length(opc_shared),
    "| bulk unique:", length(opc_unique_bulk),
    "| sc unique:", length(opc_unique_sc), "\n")

# EULER DIAGRAMS - PROPORTIONAL VENN DIAGRAMS FOR EACH PAIRWISE COMPARISON
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

# DISPLAY IN RSTUDIO THEN SAVE TO FILE
plot(fit_opalin,  quantities = TRUE, fills = c(col_bulk, col_opalin),  edges = FALSE, main = "Bulk vs Opalin+")
png(file.path(out_dir, "euler_bulk_vs_opalin.png"), width = 800, height = 700, res = 150)
plot(fit_opalin,  quantities = TRUE, fills = c(col_bulk, col_opalin),  edges = FALSE, main = "Bulk vs Opalin+")
dev.off()

plot(fit_plekhg1, quantities = TRUE, fills = c(col_bulk, col_plekhg1), edges = FALSE, main = "Bulk vs Plekhg1+")
png(file.path(out_dir, "euler_bulk_vs_plekhg1.png"), width = 800, height = 700, res = 150)
plot(fit_plekhg1, quantities = TRUE, fills = c(col_bulk, col_plekhg1), edges = FALSE, main = "Bulk vs Plekhg1+")
dev.off()

plot(fit_opc,     quantities = TRUE, fills = c(col_bulk, col_opc),     edges = FALSE, main = "Bulk vs OPCs")
png(file.path(out_dir, "euler_bulk_vs_opc.png"), width = 800, height = 700, res = 150)
plot(fit_opc,     quantities = TRUE, fills = c(col_bulk, col_opc),     edges = FALSE, main = "Bulk vs OPCs")
dev.off()

# 4-WAY VENN DIAGRAM SHOWING OVERLAP ACROSS ALL FOUR DATASETS TOGETHER
# USING ggVennDiagram AS eulerr STRUGGLES WITH 4-WAY PROPORTIONAL DIAGRAMS
peak_sets <- list(
  Bulk    = bulk_peaks,
  Opalin  = opalin_peaks,
  Plekhg1 = plekhg1_peaks,
  OPCs    = opc_peaks
)

p_venn <- ggVennDiagram(peak_sets,
                        label_alpha = 0,
                        set_color   = c(col_bulk, col_opalin, col_plekhg1, col_opc)) +
  scale_fill_gradient(low = "white", high = col_bulk) +
  labs(title = "Peak overlap across bulk and SC oligodendrocyte subtypes") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

print(p_venn)
ggsave(file.path(out_dir, "venn_all_datasets.png"),
       plot = p_venn, width = 10, height = 8, dpi = 150)

# SAVE ALL OBJECTS FOR USE IN DOWNSTREAM ANALYSES
save(bulk_counts, oligo_bulk_data, bulk_meta, sample_mapping,
     opalin_vst, plekhg1_vst, opc_vst,
     bulk_peaks, opalin_peaks, plekhg1_peaks, opc_peaks,
     bulk_gr, opalin_gr, plekhg1_gr, opc_gr,
     overlaps_opalin, overlaps_plekhg1, overlaps_opc,
     opalin_shared, opalin_unique_bulk, opalin_unique_sc,
     plekhg1_shared, plekhg1_unique_bulk, plekhg1_unique_sc,
     opc_shared, opc_unique_bulk, opc_unique_sc,
     file = "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_granges_vst.RData")



dim(opalin_vst)
head(rownames(opalin_vst))




#use colourblind palettte ( i need to verify this ?)
# pairwise using ggVennDiagram
p_opalin <- ggVennDiagram(
  list(Bulk = bulk_peaks, "Opalin+" = opalin_peaks),
  label_alpha = 0,
  set_color = c(col_bulk, col_opalin)
) +
  scale_fill_gradient(low = "white", high = col_opalin) +
  labs(title = "Bulk vs Opalin+")

print(p_opalin)


#TRYING AN UPSET PLOT

library(UpSetR)
# use the overlap results we already calculated
# create binary membership based on coordinate overlaps


# first create GRanges for all unique coordinate-based peaks

# get SC peaks that overlap with bulk
opalin_in_bulk  <- opalin_peaks[unique(subjectHits(overlaps_opalin))]
plekhg1_in_bulk <- plekhg1_peaks[unique(subjectHits(overlaps_plekhg1))]
opc_in_bulk     <- opc_peaks[unique(subjectHits(overlaps_opc))]

# get SC unique peaks (not in bulk)
opalin_only  <- opalin_peaks[!seq_along(opalin_peaks) %in% subjectHits(overlaps_opalin)]
plekhg1_only <- plekhg1_peaks[!seq_along(plekhg1_peaks) %in% subjectHits(overlaps_plekhg1)]
opc_only     <- opc_peaks[!seq_along(opc_peaks) %in% subjectHits(overlaps_opc)]

# build upset data frame
# for bulk peaks - check which SC subtypes they overlap with
upset_bulk <- data.frame(
  Bulk    = rep(1, length(bulk_peaks)),
  Opalin  = as.integer(bulk_peaks %in% opalin_shared),
  Plekhg1 = as.integer(bulk_peaks %in% plekhg1_shared),
  OPCs    = as.integer(bulk_peaks %in% opc_shared)
)

# for SC unique peaks - they are only in their respective dataset
upset_opalin_only <- data.frame(
  Bulk    = rep(0, length(opalin_only)),
  Opalin  = rep(1, length(opalin_only)),
  Plekhg1 = rep(0, length(opalin_only)),
  OPCs    = rep(0, length(opalin_only))
)

upset_plekhg1_only <- data.frame(
  Bulk    = rep(0, length(plekhg1_only)),
  Opalin  = rep(0, length(plekhg1_only)),
  Plekhg1 = rep(1, length(plekhg1_only)),
  OPCs    = rep(0, length(plekhg1_only))
)

upset_opc_only <- data.frame(
  Bulk    = rep(0, length(opc_only)),
  Opalin  = rep(0, length(opc_only)),
  Plekhg1 = rep(0, length(opc_only)),
  OPCs    = rep(1, length(opc_only))
)

# combine all rows
upset_df_full <- rbind(upset_bulk, upset_opalin_only, upset_plekhg1_only, upset_opc_only)

# plot
upset(upset_df_full,
      sets            = c("Bulk", "Opalin", "Plekhg1", "OPCs"),
      order.by        = "freq",
      sets.bar.color  = c(col_bulk, col_opalin, col_plekhg1, col_opc),
      text.scale      = 1.3,
      mainbar.y.label = "Number of peaks",
      sets.x.label    = "Total peaks per dataset")





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
