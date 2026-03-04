load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_pseudobulked.RData")
.libPaths(c("/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
.libPaths()

library(GenomicRanges)

#convert to granges
bulk_peaks <- rownames(bulk_counts)
bulk_gr <- GRanges(
  seqnames = gsub("_.*", "", gsub("chr", "", bulk_peaks)),
  ranges   = IRanges(
    start = as.numeric(gsub(".*_(\\d+)_\\d+$", "\\1", bulk_peaks)),
    end   = as.numeric(gsub(".*_(\\d+)$", "\\1", bulk_peaks))
  )
)

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

#find overlaps between peaks for each oligo bulk vs olgio subtype
overlaps_opalin  <- findOverlaps(bulk_gr, opalin_gr)
overlaps_plekhg1 <- findOverlaps(bulk_gr, plekhg1_gr)
overlaps_opc     <- findOverlaps(bulk_gr, opc_gr)

#classification of peaks
opalin_shared      <- bulk_peaks[unique(queryHits(overlaps_opalin))]
opalin_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_opalin)]
opalin_unique_sc   <- opalin_peaks[!seq_along(opalin_peaks) %in% subjectHits(overlaps_opalin)]

plekhg1_shared      <- bulk_peaks[unique(queryHits(overlaps_plekhg1))]
plekhg1_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_plekhg1)]
plekhg1_unique_sc   <- plekhg1_peaks[!seq_along(plekhg1_peaks) %in% subjectHits(overlaps_plekhg1)]

opc_shared      <- bulk_peaks[unique(queryHits(overlaps_opc))]
opc_unique_bulk <- bulk_peaks[!seq_along(bulk_peaks) %in% queryHits(overlaps_opc)]
opc_unique_sc   <- opc_peaks[!seq_along(opc_peaks) %in% subjectHits(overlaps_opc)]

#summarise peak figures
cat("Opalin+  - shared:", length(opalin_shared), 
    "| bulk unique:", length(opalin_unique_bulk), 
    "| sc unique:", length(opalin_unique_sc), "\n")
cat("Plekhg1+ - shared:", length(plekhg1_shared), 
    "| bulk unique:", length(plekhg1_unique_bulk), 
    "| sc unique:", length(plekhg1_unique_sc), "\n")
cat("OPCs     - shared:", length(opc_shared), 
    "| bulk unique:", length(opc_unique_bulk), 
    "| sc unique:", length(opc_unique_sc), "\n")

#visualisation (venn diagram)
library(ggVennDiagram)
library(ggplot2)
library(eulerr)

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

plot(fit_opalin,  quantities = TRUE, fills = c("steelblue", "#E69F00"), edges = FALSE, main = "Bulk vs Opalin+")
plot(fit_plekhg1, quantities = TRUE, fills = c("steelblue", "#009E73"), edges = FALSE, main = "Bulk vs Plekhg1+")
plot(fit_opc,     quantities = TRUE, fills = c("steelblue", "#CC79A7"),edges = FALSE, main = "Bulk vs OPCs")

#ANALYSIS, LOOKING AT PEAK SIZE

#calculate the width of the peaks (GRanges)
bulk_width_df <- data.frame(
  peak     = bulk_peaks,
  width    = width(bulk_gr),
  category = ifelse(bulk_peaks %in% opalin_shared, "Shared", "Bulk unique")
)

opalin_width_df <- data.frame(
  peak     = opalin_peaks,
  width    = width(opalin_gr),
  category = ifelse(seq_along(opalin_peaks) %in% subjectHits(overlaps_opalin), "Shared", "Opalin+ unique")
)

# combine
width_df <- rbind(bulk_width_df, opalin_width_df)

# plot
width_df_unique <- width_df[width_df$category != "Shared", ]

ggplot(width_df_unique, aes(x = category, y = width, fill = category)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Bulk unique" = "steelblue", "Opalin+ unique" = "#E69F00")) +
  coord_cartesian(ylim = c(0, 6000)) +
  labs(title = "Peak width distribution: Bulk unique vs Opalin+ unique",
       x = "", y = "Peak width (bp)") +
  theme_classic() +
  theme(legend.position = "none")

#do the same for remaning subtypes

# plekhg1 width df
plekhg1_width_df <- data.frame(
  peak     = plekhg1_peaks,
  width    = width(plekhg1_gr),
  category = ifelse(seq_along(plekhg1_peaks) %in% subjectHits(overlaps_plekhg1), "Shared", "Plekhg1+ unique")
)

bulk_width_plekhg1 <- data.frame(
  peak     = bulk_peaks,
  width    = width(bulk_gr),
  category = ifelse(bulk_peaks %in% plekhg1_shared, "Shared", "Bulk unique")
)

# opc width df
opc_width_df <- data.frame(
  peak     = opc_peaks,
  width    = width(opc_gr),
  category = ifelse(seq_along(opc_peaks) %in% subjectHits(overlaps_opc), "Shared", "OPC unique")
)

bulk_width_opc <- data.frame(
  peak     = bulk_peaks,
  width    = width(bulk_gr),
  category = ifelse(bulk_peaks %in% opc_shared, "Shared", "Bulk unique")
)

# combine and filter out shared
plekhg1_df_unique <- rbind(bulk_width_plekhg1, plekhg1_width_df)
plekhg1_df_unique <- plekhg1_df_unique[plekhg1_df_unique$category != "Shared", ]

opc_df_unique <- rbind(bulk_width_opc, opc_width_df)
opc_df_unique <- opc_df_unique[opc_df_unique$category != "Shared", ]

# plots
ggplot(plekhg1_df_unique, aes(x = category, y = width, fill = category)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Bulk unique" = "steelblue", "Plekhg1+ unique" = "#009E73")) +
  coord_cartesian(ylim = c(0, 10000)) +
  labs(title = "Peak width distribution: Bulk unique vs Plekhg1+ unique",
       x = "", y = "Peak width (bp)") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(opc_df_unique, aes(x = category, y = width, fill = category)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Bulk unique" = "steelblue", "OPC unique" = "#CC79A7")) +
  coord_cartesian(ylim = c(0, 10000)) +
  labs(title = "Peak width distribution: Bulk unique vs OPC unique",
       x = "", y = "Peak width (bp)") +
  theme_classic() +
  theme(legend.position = "none")

#check rnage of peak szies (graphs need fixing)
summary(opc_df_unique$width[opc_df_unique$category == "OPC unique"])


#ANALYSIS 3 look at genomic context and annotate peaks
#use ChIPseeker
#issue intstalling these packages
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")


