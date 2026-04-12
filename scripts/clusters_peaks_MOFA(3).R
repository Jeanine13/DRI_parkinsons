
#MOFA ON NEW SC VST NORMALISED DATA 

#make sure to run chipseeker annotation script before doing GO enrichment/rGREAT
library(MOFA2)

load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_granges_vst.RData")

#adding output directory for new sc data models
mofa_out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/new_sc_mofa"
dir.create(mofa_out_dir, recursive = TRUE, showWarnings = FALSE)

# SAMPLE METADATA FOR PD ASSOCIATION TESTING
sample_meta <- data.frame(
  sample    = c("PD1222", "PD1231", "PD726", "PD833", "PD936",
                "PDC167", "PDC184", "PDC197", "PDC87"),
  condition = c("PD", "PD", "PD", "PD", "PD",
                "Control", "Control", "Control", "Control")
)

# select top 5000 HVFs per view
bulk_vars    <- apply(bulk_counts, 1, var)
bulk_top_idx <- head(order(bulk_vars, decreasing = TRUE), 5000)
bulk_hvf     <- bulk_counts[bulk_top_idx, ]

opalin_vars    <- apply(opalin_vst, 1, var)
opalin_top_idx <- head(order(opalin_vars, decreasing = TRUE), 5000)
opalin_hvf     <- opalin_vst[opalin_top_idx, ]

plekhg1_vars    <- apply(plekhg1_vst, 1, var)
plekhg1_top_idx <- head(order(plekhg1_vars, decreasing = TRUE), 5000)
plekhg1_hvf     <- plekhg1_vst[plekhg1_top_idx, ]

# CHECK DIMENSIONS to verify top 5000 HVF have been selected
dim(bulk_hvf)
dim(opalin_hvf)
dim(plekhg1_hvf)
dim(opc_hvf)



# create MOFA object 1 (bulk vs opalin)
mofa_input_opalin <- list(
  bulk   = as.matrix(bulk_hvf),
  opalin = as.matrix(opalin_hvf)
)

mofa_opalin_vst <- create_mofa(mofa_input_opalin)

data_opts  <- get_default_data_options(mofa_opalin_vst)
model_opts <- get_default_model_options(mofa_opalin_vst)
model_opts$num_factors <- 9

train_opts <- get_default_training_options(mofa_opalin_vst)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

mofa_opalin_vst <- prepare_mofa(mofa_opalin_vst,
                                data_options     = data_opts,
                                model_options    = model_opts,
                                training_options = train_opts)

mofa_opalin_vst <- run_mofa(mofa_opalin_vst,
                            outfile      = "mofa_opalin_newsc.hdf5",
                            use_basilisk = TRUE)
 


sample_meta <- data.frame(
  sample    = c("PD1222", "PD1231", "PD726", "PD833", "PD936",
                "PDC167", "PDC184", "PDC197", "PDC87"),
  condition = c("PD", "PD", "PD", "PD", "PD",
                "Control", "Control", "Control", "Control")
)

samples_metadata(mofa_opalin_vst) <- sample_meta

plot_variance_explained(mofa_opalin_vst, max_r2 = 15)
plot_factor_cor(mofa_opalin_vst)



#check if there is PD association 

plot_factor(mofa_opalin_vst, factors = 1:7, color_by = "condition")


#FACTOR WEIGHTS PLOTTING ANF IN A TABLE

#adding peak key,need to help match peaks between the two dataframes
bulk_anno_df$peak_key <- paste(bulk_anno_df$seqnames, bulk_anno_df$start, bulk_anno_df$end, sep = "_")

# FACTOR 1 - BULK VIEW (EXPLAINS MOST VARIANCE IN BULK)
bulk_weights_f1 <- get_weights(mofa_opalin_vst, views = "bulk", factors = 1, as.data.frame = TRUE)
bulk_weights_f1 <- bulk_weights_f1[order(abs(bulk_weights_f1$value), decreasing = TRUE), ]
top_peaks_f1_annotated <- merge(head(bulk_weights_f1, 20),
                                bulk_anno_df[, c("peak_key", "annotation", "GENENAME", "distanceToTSS")],
                                by.x = "feature", by.y = "peak_key")
top_peaks_f1_annotated <- top_peaks_f1_annotated[!duplicated(top_peaks_f1_annotated$feature), ]
top_peaks_f1_annotated$view <- "bulk"
write.csv(top_peaks_f1_annotated[, c("feature", "value", "view", "annotation", "GENENAME", "distanceToTSS")],
          file.path(mofa_out_dir, "mofa_opalin_factor1_bulk_top20.csv"), row.names = FALSE)

#quantify promoters, intronic, intergenic
table(bulk_anno_df$annotation[bulk_anno_df$peak_key %in% top_bulk_f1_peaks$feature])
# FACTOR 2 - OPALIN VIEW (EXPLAINS MOST VARIANCE IN OPALIN)



opalin_weights_f2 <- get_weights(mofa_opalin_vst, views = "opalin", factors = 2, as.data.frame = TRUE)
opalin_weights_f2 <- opalin_weights_f2[order(abs(opalin_weights_f2$value), decreasing = TRUE), ]
top_peaks_f2_annotated <- merge(head(opalin_weights_f2, 20),
                                opalin_anno_cp[, c("peak_key", "annotation", "GENENAME.x", "distanceToTSS")],
                                by.x = "feature", by.y = "peak_key")
top_peaks_f2_annotated <- top_peaks_f2_annotated[!duplicated(top_peaks_f2_annotated$feature), ]
top_peaks_f2_annotated$view <- "opalin"
write.csv(top_peaks_f2_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")],
          file.path(mofa_out_dir, "mofa_opalin_factor2_opalin_top20.csv"), row.names = FALSE)

# FACTOR 3 - OPALIN VIEW (EXPLAINS MOST VARIANCE IN OPALIN)
opalin_weights_f3 <- get_weights(mofa_opalin_vst, views = "opalin", factors = 3, as.data.frame = TRUE)
opalin_weights_f3 <- opalin_weights_f3[order(abs(opalin_weights_f3$value), decreasing = TRUE), ]
top_peaks_f3_annotated <- merge(head(opalin_weights_f3, 20),
                                opalin_anno_cp[, c("peak_key", "annotation", "GENENAME.x", "distanceToTSS")],
                                by.x = "feature", by.y = "peak_key")
top_peaks_f3_annotated <- top_peaks_f3_annotated[!duplicated(top_peaks_f3_annotated$feature), ]
top_peaks_f3_annotated$view <- "opalin"
write.csv(top_peaks_f3_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")],
          file.path(mofa_out_dir, "mofa_opalin_factor3_opalin_top20.csv"), row.names = FALSE)

# FACTOR 4 - BOTH VIEWS (SHARED FACTOR - MOST BIOLOGICALLY INTERESTING)
bulk_weights_f4 <- get_weights(mofa_opalin_vst, views = "bulk", factors = 4, as.data.frame = TRUE)
bulk_weights_f4 <- bulk_weights_f4[order(abs(bulk_weights_f4$value), decreasing = TRUE), ]
top_peaks_f4_bulk_annotated <- merge(head(bulk_weights_f4, 20),
                                     bulk_anno_df[, c("peak_key", "annotation", "GENENAME", "distanceToTSS")],
                                     by.x = "feature", by.y = "peak_key")
top_peaks_f4_bulk_annotated <- top_peaks_f4_bulk_annotated[!duplicated(top_peaks_f4_bulk_annotated$feature), ]
top_peaks_f4_bulk_annotated$view <- "bulk"
write.csv(top_peaks_f4_bulk_annotated[, c("feature", "value", "view", "annotation", "GENENAME", "distanceToTSS")],
          file.path(mofa_out_dir, "mofa_opalin_factor4_bulk_top20.csv"), row.names = FALSE)

opalin_weights_f4 <- get_weights(mofa_opalin_vst, views = "opalin", factors = 4, as.data.frame = TRUE)
opalin_weights_f4 <- opalin_weights_f4[order(abs(opalin_weights_f4$value), decreasing = TRUE), ]
top_peaks_f4_opalin_annotated <- merge(head(opalin_weights_f4, 20),
                                       opalin_anno_cp[, c("peak_key", "annotation", "GENENAME.x", "distanceToTSS")],
                                       by.x = "feature", by.y = "peak_key")
top_peaks_f4_opalin_annotated <- top_peaks_f4_opalin_annotated[!duplicated(top_peaks_f4_opalin_annotated$feature), ]
top_peaks_f4_opalin_annotated$view <- "opalin"
write.csv(top_peaks_f4_opalin_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")],
          file.path(mofa_out_dir, "mofa_opalin_factor4_opalin_top20.csv"), row.names = FALSE)

# FACTOR 5 - BULK VIEW
bulk_weights_f5 <- get_weights(mofa_opalin_vst, views = "bulk", factors = 5, as.data.frame = TRUE)
bulk_weights_f5 <- bulk_weights_f5[order(abs(bulk_weights_f5$value), decreasing = TRUE), ]
top_peaks_f5_annotated <- merge(head(bulk_weights_f5, 20),
                                bulk_anno_df[, c("peak_key", "annotation", "GENENAME", "distanceToTSS")],
                                by.x = "feature", by.y = "peak_key")
top_peaks_f5_annotated <- top_peaks_f5_annotated[!duplicated(top_peaks_f5_annotated$feature), ]
top_peaks_f5_annotated$view <- "bulk"
write.csv(top_peaks_f5_annotated[, c("feature", "value", "view", "annotation", "GENENAME", "distanceToTSS")],
          file.path(mofa_out_dir, "mofa_opalin_factor5_bulk_top20.csv"), row.names = FALSE)

# FACTOR 6 - OPALIN VIEW
opalin_weights_f6 <- get_weights(mofa_opalin_vst, views = "opalin", factors = 6, as.data.frame = TRUE)
opalin_weights_f6 <- opalin_weights_f6[order(abs(opalin_weights_f6$value), decreasing = TRUE), ]
top_peaks_f6_annotated <- merge(head(opalin_weights_f6, 20),
                                opalin_anno_cp[, c("peak_key", "annotation", "GENENAME.x", "distanceToTSS")],
                                by.x = "feature", by.y = "peak_key")
top_peaks_f6_annotated <- top_peaks_f6_annotated[!duplicated(top_peaks_f6_annotated$feature), ]
top_peaks_f6_annotated$view <- "opalin"
write.csv(top_peaks_f6_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")],
          file.path(mofa_out_dir, "mofa_opalin_factor6_opalin_top20.csv"), row.names = FALSE)

# SAVE MOFA OBJECT
saveRDS(mofa_opalin_vst,
        file.path(mofa_out_dir, "mofa_opalin_vst_object.rds"))

cat("All factor weight CSVs saved to:", mofa_out_dir, "\n")



#load saved mofa model
mofa_opalin_vst <- readRDS(
  file.path(mofa_out_dir, "mofa_opalin_vst_object.rds")
)


## BACKGROUND UNIVERSE PEAK SETS FOR GO ENRICHMENT
# LOAD BACKGROUND PEAK SETS - THESE ARE USED AS THE UNIVERSE IN GO ENRICHMENT
# USING THE FULL CONSENSUS PEAK SET AS BACKGROUND IS MORE APPROPRIATE THAN
# USING ALL HUMAN GENES AS IT RESTRICTS THE COMPARISON TO ACCESSIBLE REGIONS



sc_consensus  <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/single_cell_consensus_peaks.rds")
bulk_consensus <- readRDS("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/student_data_package/data_in/union_peaks/bulk_union_peaks.rds")

# FIX CHROMOSOME NAMING TO MATCH TXDB
seqlevelsStyle(sc_consensus)   <- "UCSC"
seqlevelsStyle(bulk_consensus) <- "UCSC"

# ANNOTATE CONSENSUS PEAKS TO GET BACKGROUND GENE LISTS
sc_consensus_anno   <- as.data.frame(annotatePeak(sc_consensus,   tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db"))
bulk_consensus_anno <- as.data.frame(annotatePeak(bulk_consensus, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db"))

# EXTRACT UNIQUE ENTREZ IDs AS BACKGROUND GENE LISTS
sc_background_entrez   <- unique(sc_consensus_anno$geneId[!is.na(sc_consensus_anno$geneId)])
bulk_background_entrez <- unique(bulk_consensus_anno$geneId[!is.na(bulk_consensus_anno$geneId)])

cat("SC background genes:", length(sc_background_entrez), "\n")
cat("Bulk background genes:", length(bulk_background_entrez), "\n")

# GO ENRICHMENT FUNCTION
# background_entrez IS NOW A REQUIRED ARGUMENT TO ENSURE CORRECT UNIVERSE IS USED
# mode CONTROLS WHETHER TO USE POSITIVE, NEGATIVE OR ABSOLUTE WEIGHTS
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
  
  write.csv(
    data.frame(ego),
    file.path(mofa_out_dir, paste0("GO_mofa_opalin_", label, ".csv")),
    row.names = FALSE
  )
  
  ggsave(
    file.path(mofa_out_dir, paste0("dotplot_GO_mofa_opalin_", label, ".png")),
    plot = dotplot(ego, showCategory = 30) + labs(title = paste("GO -", label)),
    width = 10, height = 14, dpi = 150
  )
  
  return(ego)
}

# RUN GO FOR FACTORS 2 AND 3 USING SC BACKGROUND
# FACTORS 2 AND 3 EXPLAIN MOST VARIANCE IN OPALIN VIEW
ego_f2_opalin_pos <- run_go_factor(opalin_weights_f2, opalin_anno_cp, "peak_key", "geneId",
                                   "factor2_opalin_positive",
                                   background_entrez = sc_background_entrez,
                                   mode = "positive")

ego_f2_opalin_neg <- run_go_factor(opalin_weights_f2, opalin_anno_cp, "peak_key", "geneId",
                                   "factor2_opalin_negative",
                                   background_entrez = sc_background_entrez,
                                   mode = "negative")

ego_f3_opalin_pos <- run_go_factor(opalin_weights_f3, opalin_anno_cp, "peak_key", "geneId",
                                   "factor3_opalin_positive",
                                   background_entrez = sc_background_entrez,
                                   mode = "positive")

ego_f3_opalin_neg <- run_go_factor(opalin_weights_f3, opalin_anno_cp, "peak_key", "geneId",
                                   "factor3_opalin_negative",
                                   background_entrez = sc_background_entrez,
                                   mode = "negative")

# RUN GO FOR FACTOR 1 USING BULK BACKGROUND
# FACTOR 1 EXPLAINS MOST VARIANCE IN BULK VIEW
ego_f1_bulk_pos <- run_go_factor(bulk_weights_f1, bulk_anno_df, "peak_key", "geneId",
                                 "factor1_bulk_positive",
                                 background_entrez = bulk_background_entrez,
                                 mode = "positive")

ego_f1_bulk_neg <- run_go_factor(bulk_weights_f1, bulk_anno_df, "peak_key", "geneId",
                                 "factor1_bulk_negative",
                                 background_entrez = bulk_background_entrez,
                                 mode = "negative")

# RUN GO FOR FACTOR 4 USING BOTH BACKGROUNDS
# FACTOR 4 IS THE SHARED FACTOR 
ego_f4_bulk_pos <- run_go_factor(bulk_weights_f4, bulk_anno_df, "peak_key", "geneId",
                                 "factor4_bulk_positive",
                                 background_entrez = bulk_background_entrez,
                                 mode = "positive")

ego_f4_bulk_neg <- run_go_factor(bulk_weights_f4, bulk_anno_df, "peak_key", "geneId",
                                 "factor4_bulk_negative",
                                 background_entrez = bulk_background_entrez,
                                 mode = "negative")

ego_f4_opalin_pos <- run_go_factor(opalin_weights_f4, opalin_anno_cp, "peak_key", "geneId",
                                   "factor4_opalin_positive",
                                   background_entrez = sc_background_entrez,
                                   mode = "positive")

ego_f4_opalin_neg <- run_go_factor(opalin_weights_f4, opalin_anno_cp, "peak_key", "geneId",
                                   "factor4_opalin_negative",
                                   background_entrez = sc_background_entrez,
                                   mode = "negative")

# CHECK HOW MANY SIGNIFICANT TERMS WERE FOUND PER ENRICHMENT
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

head(data.frame(ego_f1_bulk_pos)[, c("Description", "GeneRatio", "p.adjust")])
head(data.frame(ego_f3_opalin_pos)[, c("Description", "GeneRatio", "p.adjust")])
head(data.frame(ego_f4_bulk_pos)[, c("Description", "GeneRatio", "p.adjust")])


print(dotplot(ego_f1_bulk_pos, showCategory = 20) + 
        labs(title = "GO - Factor 1 bulk positive weights") +
        theme(axis.text.y = element_text(size = 7),
              axis.text.x = element_text(size = 7),
              plot.title  = element_text(size = 10)))

print(dotplot(ego_f2_opalin_neg, showCategory = 4) + 
        labs(title = "GO - Factor 2 Opalin+ negative weights") +
        theme(axis.text.y = element_text(size = 7),
              axis.text.x = element_text(size = 7),
              plot.title  = element_text(size = 10)) +
        scale_size_continuous(range = c(2, 8)))



print(dotplot(ego_f3_opalin_pos, showCategory = 13) + 
        labs(title = "GO - Factor 3 Opalin+ positive weights"))

print(dotplot(ego_f4_bulk_pos, showCategory = 6) + 
        labs(title = "GO - Factor 4 bulk positive weights"))

#TRYING RGREAT INSTEAD OF CHIPSEEKER

library(rGREAT)

# rGREAT TAKES GRANGES OBJECTS DIRECTLY AS INPUT
# WE NEED TO CONVERT OUR TOP WEIGHTED PEAK NAMES BACK TO GRANGES

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
    seqnames = paste0("chr", sapply(strsplit(peaks, "-"), `[`, 1)),
    ranges   = IRanges(
      start = as.numeric(sapply(strsplit(peaks, "-"), `[`, 2)),
      end   = as.numeric(sapply(strsplit(peaks, "-"), `[`, 3))
    )
  )
}

# TEST WITH FACTOR 1 BULK POSITIVE WEIGHTS
top_bulk_f1_peaks <- bulk_weights_f1[bulk_weights_f1$value > 0, ]
top_bulk_f1_peaks <- head(top_bulk_f1_peaks[order(top_bulk_f1_peaks$value, decreasing = TRUE), ], 500)

top_bulk_f1_gr <- bulk_peaks_to_gr(top_bulk_f1_peaks$feature)

# RUN GREAT ENRICHMENT
great_f1_bulk <- great(top_bulk_f1_gr,
                       gene_sets   = "GO:BP",
                       tss_source  = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       background  = bulk_consensus,
                       cores       = 1)

# VIEW RESULTS
tb <- getEnrichmentTable(great_f1_bulk)
head(tb[, c("id", "description", "p_adjust", "fold_enrichment")])


tb[tb$p_adjust < 0.05, c("id", "description", "p_adjust", "fold_enrichment", "observed_region_hits")]

#CREATE TABLE
tb_sig <- tb[tb$p_adjust < 0.05, ]
tb_sig <- tb_sig[order(tb_sig$p_adjust), ]

# CHECK HOW MANY SIGNIFICANT TERMS
nrow(tb_sig)
write.csv(tb_sig[, c("id", "description", "fold_enrichment", "observed_region_hits", "p_adjust")],
          file.path(mofa_out_dir, "rGREAT_factor1_bulk_positive.csv"),
          row.names = FALSE)


#factor 2 - 4-great
# CONVERT FEATURE COLUMN TO CHARACTER FIRST
opalin_weights_f2$feature <- as.character(opalin_weights_f2$feature)
opalin_weights_f3$feature <- as.character(opalin_weights_f3$feature)
opalin_weights_f4$feature <- as.character(opalin_weights_f4$feature)
bulk_weights_f1$feature   <- as.character(bulk_weights_f1$feature)
bulk_weights_f4$feature   <- as.character(bulk_weights_f4$feature)

# FUNCTION TO RUN rGREAT FOR A GIVEN FACTOR AND VIEW
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
  
  # CONVERT PEAK NAMES TO GRANGES DEPENDING ON VIEW
  if (view == "bulk") {
    peaks_gr <- bulk_peaks_to_gr(top_peaks$feature)
  } else {
    peaks_gr <- sc_peaks_to_gr(top_peaks$feature)
  }
  
  
  # RUN GREAT
  great_res <- great(peaks_gr,
                     gene_sets  = "GO:BP",
                     tss_source = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                     background = background_gr,
                     cores      = 1)
  
  tb <- getEnrichmentTable(great_res)
  tb_sig <- tb[tb$p_adjust < 0.05, ]
  tb_sig <- tb_sig[order(tb_sig$p_adjust), ]
  
  cat(label, "- significant terms:", nrow(tb_sig), "\n")
  
  write.csv(tb_sig[, c("id", "description", "fold_enrichment", "observed_region_hits", "p_adjust")],
            file.path(mofa_out_dir, paste0("rGREAT_", label, ".csv")),
            row.names = FALSE)
  
  return(tb_sig)
}

# FACTOR 2 OPALIN NEGATIVE (4 clusterProfiler terms)
great_f2_opalin_neg <- run_great_factor(opalin_weights_f2, "opalin", 2,
                                        sc_consensus, "factor2_opalin_negative",
                                        mode = "negative")

# FACTOR 3 OPALIN POSITIVE (13 clusterProfiler terms)
great_f3_opalin_pos <- run_great_factor(opalin_weights_f3, "opalin", 3,
                                        sc_consensus, "factor3_opalin_positive",
                                        mode = "positive")

# FACTOR 4 BULK POSITIVE (6 clusterProfiler terms)
great_f4_bulk_pos <- run_great_factor(bulk_weights_f4, "bulk", 4,
                                      bulk_consensus, "factor4_bulk_positive",
                                      mode = "positive")

# FACTOR 4 OPALIN POSITIVE (1 clusterProfiler term)
great_f4_opalin_pos <- run_great_factor(opalin_weights_f4, "opalin", 4,
                                        sc_consensus, "factor4_opalin_positive",
                                        mode = "positive")

#ran successful






#MOFA MODEL 2
#MOFA FOR PLEKHG1 VS BULK
# SELECT TOP 5000 HVFs FOR PLEKHG1
plekhg1_vars    <- apply(plekhg1_vst, 1, var)
plekhg1_top_idx <- head(order(plekhg1_vars, decreasing = TRUE), 5000)
plekhg1_hvf     <- plekhg1_vst[plekhg1_top_idx, ]

# CREATE 2-VIEW MOFA OBJECT: BULK VS PLEKHG1+
mofa_input_plekhg1 <- list(
  bulk    = as.matrix(bulk_hvf),
  plekhg1 = as.matrix(plekhg1_hvf)
)

mofa_plekhg1_vst <- create_mofa(mofa_input_plekhg1)

data_opts  <- get_default_data_options(mofa_plekhg1_vst)
model_opts <- get_default_model_options(mofa_plekhg1_vst)
model_opts$num_factors <- 9

train_opts <- get_default_training_options(mofa_plekhg1_vst)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

mofa_plekhg1_vst <- prepare_mofa(mofa_plekhg1_vst,
                                 data_options     = data_opts,
                                 model_options    = model_opts,
                                 training_options = train_opts)

mofa_plekhg1_vst <- run_mofa(mofa_plekhg1_vst,
                             outfile      = file.path(mofa_out_dir, "mofa_plekhg1_newsc.hdf5"),
                             use_basilisk = TRUE)


samples_metadata(mofa_plekhg1_vst) <- sample_meta

plot_variance_explained(mofa_plekhg1_vst, max_r2 = 15)
plot_factor_cor(mofa_plekhg1_vst)

# CONVERT FEATURE COLUMN TO CHARACTER TO AVOID FACTOR LEVEL ISSUES
get_weights_plekhg1 <- function(factor_num, view) {
  w <- get_weights(mofa_plekhg1_vst, views = view, factors = factor_num, as.data.frame = TRUE)
  w$feature <- as.character(w$feature)
  w
}

# FACTOR 1 - PLEKHG1 VIEW
plekhg1_weights_f1 <- get_weights_plekhg1(1, "plekhg1")

# FACTOR 2 - BULK VIEW
bulk_weights_plekhg1_f2 <- get_weights_plekhg1(2, "bulk")

# FACTOR 3 - BULK VIEW
bulk_weights_plekhg1_f3 <- get_weights_plekhg1(3, "bulk")

# FACTOR 4 - BULK VIEW
bulk_weights_plekhg1_f4 <- get_weights_plekhg1(4, "bulk")
bulk_weights_plekhg1_f4$feature <- as.character(bulk_weights_plekhg1_f4$feature)

# FACTOR 5 - PLEKHG1 VIEW
plekhg1_weights_f5 <- get_weights_plekhg1(5, "plekhg1")

# CHECK DIMENSIONS
cat("Factor 1 plekhg1:", nrow(plekhg1_weights_f1), "\n")
cat("Factor 2 bulk:", nrow(bulk_weights_plekhg1_f2), "\n")
cat("Factor 3 bulk:", nrow(bulk_weights_plekhg1_f3), "\n")
cat("Factor 5 plekhg1:", nrow(plekhg1_weights_f5), "\n")

#redefine function 

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
  
  write.csv(
    data.frame(ego),
    file.path(mofa_out_dir, paste0("GO_mofa_", label, ".csv")),
    row.names = FALSE
  )
  
  ggsave(
    file.path(mofa_out_dir, paste0("dotplot_GO_mofa_", label, ".png")),
    plot = dotplot(ego, showCategory = 20) + 
      labs(title = paste("GO -", label)) +
      theme(axis.text.y = element_text(size = 7),
            axis.text.x = element_text(size = 7),
            plot.title  = element_text(size = 10)),
    width = 10, height = 14, dpi = 150
  )
  
  return(ego)
}




# ANNOTATE PLEKHG1 PEAKS WITH PEAK KEY FOR GO ENRICHMENT
plekhg1_anno_cp$peak_key <- paste(gsub("chr", "", plekhg1_anno_cp$seqnames),
                                  plekhg1_anno_cp$start,
                                  plekhg1_anno_cp$end, sep = "-")

# GO ENRICHMENT
ego_plekhg1_f1 <- run_go_factor(plekhg1_weights_f1, plekhg1_anno_cp, "peak_key", "geneId",
                                "plekhg1_factor1_plekhg1_positive",
                                background_entrez = sc_background_entrez,
                                mode = "positive")

ego_plekhg1_f2 <- run_go_factor(bulk_weights_plekhg1_f2, bulk_anno_df, "peak_key", "geneId",
                                "plekhg1_factor2_bulk_positive",
                                background_entrez = bulk_background_entrez,
                                mode = "positive")

ego_plekhg1_f3 <- run_go_factor(bulk_weights_plekhg1_f3, bulk_anno_df, "peak_key", "geneId",
                                "plekhg1_factor3_bulk_positive",
                                background_entrez = bulk_background_entrez,
                                mode = "positive")


ego_plekhg1_f4 <- run_go_factor(bulk_weights_plekhg1_f4, bulk_anno_df, "peak_key", "geneId",
                                "plekhg1_factor4_bulk_positive",
                                background_entrez = bulk_background_entrez,
                                mode = "positive")

ego_plekhg1_f5 <- run_go_factor(plekhg1_weights_f5, plekhg1_anno_cp, "peak_key", "geneId",
                                "plekhg1_factor5_plekhg1_positive",
                                background_entrez = sc_background_entrez,
                                mode = "positive")

# CHECK RESULTS
cat("Factor 1 plekhg1 positive:", nrow(data.frame(ego_plekhg1_f1)), "terms\n")
cat("Factor 2 bulk positive:", nrow(data.frame(ego_plekhg1_f2)), "terms\n")
cat("Factor 3 bulk positive:", nrow(data.frame(ego_plekhg1_f3)), "terms\n")
cat("Factor 4 bulk positive:", nrow(data.frame(ego_plekhg1_f4)), "terms\n")
cat("Factor 5 plekhg1 positive:", nrow(data.frame(ego_plekhg1_f5)), "terms\n")



#define great fucntion
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
                     cores      = 1)
  
  tb <- getEnrichmentTable(great_res)
  tb_sig <- tb[tb$p_adjust < 0.05, ]
  tb_sig <- tb_sig[order(tb_sig$p_adjust), ]
  
  cat(label, "- significant terms:", nrow(tb_sig), "\n")
  
  write.csv(tb_sig[, c("id", "description", "fold_enrichment", "observed_region_hits", "p_adjust")],
            file.path(mofa_out_dir, paste0("rGREAT_", label, ".csv")),
            row.names = FALSE)
  
  return(tb_sig)
}
# rGREAT ENRICHMENT
great_plekhg1_f1 <- run_great_factor(plekhg1_weights_f1, "plekhg1", 1,
                                     sc_consensus, "plekhg1_factor1_plekhg1_positive",
                                     mode = "positive")

great_plekhg1_f2 <- run_great_factor(bulk_weights_plekhg1_f2, "bulk", 2,
                                     bulk_consensus, "plekhg1_factor2_bulk_positive",
                                     mode = "positive")

great_plekhg1_f3 <- run_great_factor(bulk_weights_plekhg1_f3, "bulk", 3,
                                     bulk_consensus, "plekhg1_factor3_bulk_positive",


great_plekhg1_f4 <- run_great_factor(bulk_weights_plekhg1_f4, "bulk", 4,
                                      bulk_consensus, "plekhg1_factor4_bulk_positive",
                                      mode = "positive")                                                                       mode = "positive")

great_plekhg1_f5 <- run_great_factor(plekhg1_weights_f5, "plekhg1", 5,
                                     sc_consensus, "plekhg1_factor5_plekhg1_positive",
                                     mode = "positive")


#chipseeker annotations (top 20 genes)

# FUNCTION TO GET TOP 20 ANNOTATED PEAKS FOR PLEKHG1 MODEL
get_top20_plekhg1 <- function(weights_df, anno_df, peak_key_col, gene_col, label, mode = "abs") {
  
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
            file.path(mofa_out_dir, paste0("mofa_plekhg1_", label, "_top20.csv")),
            row.names = FALSE)
  
  cat("Saved:", label, "\n")
  annotated
}

# FACTOR 1 - PLEKHG1 VIEW
top20_plekhg1_f1 <- get_top20_plekhg1(plekhg1_weights_f1, plekhg1_anno_cp,
                                      "peak_key", "GENENAME.x",
                                      "factor1_plekhg1", mode = "abs")

# FACTOR 2 - BULK VIEW
top20_plekhg1_f2 <- get_top20_plekhg1(bulk_weights_plekhg1_f2, bulk_anno_df,
                                      "peak_key", "GENENAME",
                                      "factor2_bulk", mode = "abs")

# FACTOR 3 - BULK VIEW
top20_plekhg1_f3 <- get_top20_plekhg1(bulk_weights_plekhg1_f3, bulk_anno_df,
                                      "peak_key", "GENENAME",
                                      "factor3_bulk", mode = "abs")

# TOP 20 ANNOTATION
top20_plekhg1_f4 <- get_top20_plekhg1(bulk_weights_plekhg1_f4, bulk_anno_df,
                                      "peak_key", "GENENAME",
                                      "factor4_bulk", mode = "abs")
# FACTOR 5 - PLEKHG1 VIEW
top20_plekhg1_f5 <- get_top20_plekhg1(plekhg1_weights_f5, plekhg1_anno_cp,
                                      "peak_key", "GENENAME.x",
                                      "factor5_plekhg1", mode = "abs")


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


# EXTRACT OPC MODEL WEIGHTS BASED ON VARIANCE DECOMPOSITION
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

opc_anno_cp$peak_key <- paste(gsub("chr", "", opc_anno_cp$seqnames),
                              opc_anno_cp$start,
                              opc_anno_cp$end, sep = "-")

opc_anno_cp <- as.data.frame(annotatePeak(opc_gr, tssRegion = c(-3000, 3000), 
                                          TxDb = txdb, annoDb = "org.Hs.eg.db"))

opc_anno_cp <- add_gene_symbols(opc_anno_cp)

opc_anno_cp$peak_key <- paste(gsub("chr", "", opc_anno_cp$seqnames),
                              opc_anno_cp$start,
                              opc_anno_cp$end, sep = "-")

colnames(opc_anno_cp)

# TOP 20 ANNOTATED PEAKS
get_top20(bulk_weights_opc_f1, bulk_anno_df, "peak_key", "GENENAME",   "opc_factor1_bulk")
get_top20(opc_weights_f2,      opc_anno_cp,  "peak_key", "GENENAME.x", "opc_factor2_opc")
get_top20(opc_weights_f3,      opc_anno_cp,  "peak_key", "GENENAME.x", "opc_factor3_opc")
get_top20(bulk_weights_opc_f4, bulk_anno_df, "peak_key", "GENENAME",   "opc_factor4_bulk")
get_top20(bulk_weights_opc_f5, bulk_anno_df, "peak_key", "GENENAME",   "opc_factor5_bulk")

# GO ENRICHMENT
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











