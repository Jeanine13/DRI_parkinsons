

#MOFA ON NEW SC VST NORMALISED DATA 

library(MOFA2)

# select top 5000 HVFs per view
bulk_vars    <- apply(bulk_counts, 1, var)
bulk_top_idx <- head(order(bulk_vars, decreasing = TRUE), 5000)
bulk_hvf     <- bulk_counts[bulk_top_idx, ]

opalin_vars    <- apply(opalin_vst, 1, var)
opalin_top_idx <- head(order(opalin_vars, decreasing = TRUE), 5000)
opalin_hvf     <- opalin_vst[opalin_top_idx, ]

# create MOFA object (bulk vs opalin)
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
                            outfile      = "mofa_opalin_vst.hdf5",
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

bulk_weights_f1 <- get_weights(mofa_opalin_vst,
                               views = "bulk",
                               factors = 1,
                               as.data.frame = TRUE)

bulk_weights_f1 <- bulk_weights_f1[order(abs(bulk_weights_f1$value), decreasing = TRUE), ]
top_peaks_f1 <- head(bulk_weights_f1, 20)

# annotate
bulk_anno_df$peak_key <- paste(bulk_anno_df$seqnames, bulk_anno_df$start, bulk_anno_df$end, sep = "_")

top_peaks_f1_annotated <- merge(top_peaks_f1,
                                bulk_anno_df[, c("peak_key", "annotation", "GENENAME", "distanceToTSS")],
                                by.x = "feature",
                                by.y = "peak_key")

top_peaks_f1_annotated <- top_peaks_f1_annotated[!duplicated(top_peaks_f1_annotated$feature), ]

top_peaks_f1_annotated[, c("feature", "value", "annotation", "GENENAME", "distanceToTSS")]
top_peaks_f1_annotated$view <- "bulk"

write.csv(top_peaks_f1_annotated[, c("feature", "value", "view", "annotation", "GENENAME", "distanceToTSS")],
          file.path(out_dir, "mofa_opalin_vst_factor1_bulk_top20_weights.csv"),
          row.names = FALSE)

#look at opalin view for factor 1 as this factor also explain some of the variance in opalin
#rerun peak_key
opalin_anno_cp$peak_key <- paste(gsub("chr", "", opalin_anno_cp$seqnames),
                                 opalin_anno_cp$start,
                                 opalin_anno_cp$end, sep = "-")


opalin_weights_f1 <- get_weights(mofa_opalin_vst,
                                 views = "opalin",
                                 factors = 1,
                                 as.data.frame = TRUE)

opalin_weights_f1 <- opalin_weights_f1[order(abs(opalin_weights_f1$value), decreasing = TRUE), ]
top_peaks_f1_opalin <- head(opalin_weights_f1, 20)

top_peaks_f1_opalin_annotated <- merge(top_peaks_f1_opalin,
                                       opalin_anno_cp[, c("peak_key", "annotation", "GENENAME.x", "distanceToTSS")],
                                       by.x = "feature",
                                       by.y = "peak_key")

top_peaks_f1_opalin_annotated <- top_peaks_f1_opalin_annotated[!duplicated(top_peaks_f1_opalin_annotated$feature), ]
top_peaks_f1_opalin_annotated$view <- "opalin"

write.csv(top_peaks_f1_opalin_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")],
          file.path(out_dir, "mofa_opalin_vst_factor1_opalin_top20_weights.csv"),
          row.names = FALSE)

top_peaks_f1_opalin_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")]


#look at factor 2 in the oplain view, where it explain most of varaince in opalin rather than bulk
opalin_weights_f2 <- get_weights(mofa_opalin_vst,
                                 views = "opalin",
                                 factors = 2,
                                 as.data.frame = TRUE)

opalin_weights_f2 <- opalin_weights_f2[order(abs(opalin_weights_f2$value), decreasing = TRUE), ]
top_peaks_f2_opalin <- head(opalin_weights_f2, 20)

top_peaks_f2_opalin_annotated <- merge(top_peaks_f2_opalin,
                                       opalin_anno_cp[, c("peak_key", "annotation", "GENENAME.x", "distanceToTSS")],
                                       by.x = "feature",
                                       by.y = "peak_key")

top_peaks_f2_opalin_annotated <- top_peaks_f2_opalin_annotated[!duplicated(top_peaks_f2_opalin_annotated$feature), ]
top_peaks_f2_opalin_annotated$view <- "opalin"

top_peaks_f2_opalin_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")]

write.csv(top_peaks_f2_opalin_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")],
          file.path(out_dir, "mofa_opalin_vst_factor2_opalin_top20_weights.csv"),
          row.names = FALSE)

#FACTOR 2 (OPLAIN VIEW )

opalin_weights_f3 <- get_weights(mofa_opalin_vst,
                                 views = "opalin",
                                 factors = 3,
                                 as.data.frame = TRUE)

opalin_weights_f3 <- opalin_weights_f3[order(abs(opalin_weights_f3$value), decreasing = TRUE), ]
top_peaks_f3_opalin <- head(opalin_weights_f3, 20)

top_peaks_f3_opalin_annotated <- merge(top_peaks_f3_opalin,
                                       opalin_anno_cp[, c("peak_key", "annotation", "GENENAME.x", "distanceToTSS")],
                                       by.x = "feature",
                                       by.y = "peak_key")

top_peaks_f3_opalin_annotated <- top_peaks_f3_opalin_annotated[!duplicated(top_peaks_f3_opalin_annotated$feature), ]
top_peaks_f3_opalin_annotated$view <- "opalin"

top_peaks_f3_opalin_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")]

write.csv(top_peaks_f3_opalin_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")],
          file.path(out_dir, "mofa_opalin_vst_factor3_opalin_top20_weights.csv"),
          row.names = FALSE)



#GO ENRICHMENT OF TOP WEIGHTED PEAKS PER FACTOR (OPALIN VS BULK MODEL)
# function to run GO on top weighted peaks

run_go_factor <- function(weights_df, anno_df, peak_key_col, geneid_col, label,
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
    keyType       = "ENTREZID",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  write.csv(
    data.frame(ego),
    file.path(out_dir, paste0("GO_mofa_opalin_vst_", label, ".csv")),
    row.names = FALSE
  )
  
  ggsave(
    file.path(out_dir, paste0("dotplot_GO_mofa_opalin_vst_", label, ".png")),
    plot = dotplot(ego, showCategory = 30) + labs(title = paste("GO -", label)),
    width = 10, height = 14, dpi = 150
  )
  
  return(ego)
}

#running enrichment calls on factors 2 and 3 as they explain mostly opalin variance 

# Factor 2 opalin - positive and negative
ego_f2_opalin_pos <- run_go_factor(opalin_weights_f2, opalin_anno_cp, "peak_key", "geneId",
                                   "factor2_opalin_positive", mode = "positive")
ego_f2_opalin_neg <- run_go_factor(opalin_weights_f2, opalin_anno_cp, "peak_key", "geneId",
                                   "factor2_opalin_negative", mode = "negative")

# Factor 3 opalin - positive and negative
ego_f3_opalin_pos <- run_go_factor(opalin_weights_f3, opalin_anno_cp, "peak_key", "geneId",
                                   "factor3_opalin_positive", mode = "positive")
ego_f3_opalin_neg <- run_go_factor(opalin_weights_f3, opalin_anno_cp, "peak_key", "geneId",
                                   "factor3_opalin_negative", mode = "negative")


#----------------------
#plots are empty, debugging why ??

nrow(data.frame(ego_f2_opalin_pos))
nrow(data.frame(ego_f2_opalin_neg))
nrow(data.frame(ego_f3_opalin_pos))
nrow(data.frame(ego_f3_opalin_neg))

head(opalin_weights_f2$feature)
head(opalin_anno_cp$peak_key)
sum(opalin_anno_cp$peak_key %in% opalin_weights_f2$feature)


top_peaks <- opalin_weights_f2[opalin_weights_f2$value > 0, ]
top_peaks <- top_peaks[order(top_peaks$value, decreasing = TRUE), ]
top_peaks <- head(top_peaks, 300)

entrez <- unique(opalin_anno_cp$geneId[opalin_anno_cp$peak_key %in% top_peaks$feature])
entrez <- entrez[!is.na(entrez)]
length(entrez)

ego_test <- enrichGO(
  gene          = entrez,
  keyType       = "ENTREZID",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

nrow(as.data.frame(ego_test))
nrow(data.frame(ego_f2_opalin_neg))
head(data.frame(ego_f2_opalin_neg))
#-------------------------------------



#MOFA FOR PLEKHG1 VS BULK

# select HVFs for plekhg1
plekhg1_vars    <- apply(plekhg1_vst, 1, var)
plekhg1_top_idx <- head(order(plekhg1_vars, decreasing = TRUE), 5000)
plekhg1_hvf     <- plekhg1_vst[plekhg1_top_idx, ]

# create 2-view MOFA object bulk vs plekhg1
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
                             outfile      = "mofa_plekhg1_vst.hdf5",
                             use_basilisk = TRUE)


samples_metadata(mofa_plekhg1_vst) <- sample_meta

plot_variance_explained(mofa_plekhg1_vst, max_r2 = 15)
plot_factor_cor(mofa_plekhg1_vst)

plot_factor(mofa_plekhg1_vst, factors = 1:7, color_by = "condition")

#plot feature weights for FACTOR 1 (plekhg1 view)

plekhg1_weights_f1 <- get_weights(mofa_plekhg1_vst,
                                  views = "plekhg1",
                                  factors = 1,
                                  as.data.frame = TRUE)

plekhg1_weights_f1 <- plekhg1_weights_f1[order(abs(plekhg1_weights_f1$value), decreasing = TRUE), ]
top_peaks_f1_plekhg1 <- head(plekhg1_weights_f1, 20)

# create peak key for plekhg1 annotation
plekhg1_anno_cp$peak_key <- paste(gsub("chr", "", plekhg1_anno_cp$seqnames),
                                  plekhg1_anno_cp$start,
                                  plekhg1_anno_cp$end, sep = "-")

top_peaks_f1_plekhg1_annotated <- merge(top_peaks_f1_plekhg1,
                                        plekhg1_anno_cp[, c("peak_key", "annotation", "GENENAME.x", "distanceToTSS")],
                                        by.x = "feature",
                                        by.y = "peak_key")

top_peaks_f1_plekhg1_annotated <- top_peaks_f1_plekhg1_annotated[!duplicated(top_peaks_f1_plekhg1_annotated$feature), ]
top_peaks_f1_plekhg1_annotated$view <- "plekhg1"

top_peaks_f1_plekhg1_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")]

write.csv(top_peaks_f1_plekhg1_annotated[, c("feature", "value", "view", "annotation", "GENENAME.x", "distanceToTSS")],
          file.path(out_dir, "mofa_plekhg1_vst_factor1_plekhg1_top20_weights.csv"),
          row.names = FALSE)

#RUN GO ON factor 1 (plekhg1 vs bulk model)

ego_plekhg1_f1 <- run_go_factor(plekhg1_weights_f1, plekhg1_anno_cp, "peak_key", "geneId",
                                "factor1_plekhg1", mode = "abs")

nrow(data.frame(ego_plekhg1_f1))





