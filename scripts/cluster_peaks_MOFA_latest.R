# MOFA TRAINING AND WEIGHT EXTRACTION SCRIPT
# RUN AFTER cluster_peaks_chipseeker.R
# TRAINS ALL FOUR MOFA MODELS
# SAVES MODEL PLOTS + VARIANCE TABLES
# EXTRACTS TOP 20 ANNOTATED PEAKS FOR SELECTED FACTOR/VIEW COMBINATIONS
# SAVES ABSOLUTE, POSITIVE, AND NEGATIVE TOP FEATURES

.libPaths(c(
  "/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5",
  .libPaths()
))

library(MOFA2)
library(ggplot2)

# DEFINE OUTPUT DIRECTORY
mofa_out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/new_sc_mofa_clean_14042026"
dir.create(mofa_out_dir, recursive = TRUE, showWarnings = FALSE)

# LOAD SAVED CHIPSEEKER ANNOTATIONS AND INPUT DATA
load(file.path(mofa_out_dir, "chipseeker_annotations_all_datasets.RData"))
load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_granges_vst.RData")

cat("All required objects loaded\n")

# SAMPLE METADATA FOR PD ASSOCIATION TESTING
sample_meta <- data.frame(
  sample    = c("PD1222", "PD1231", "PD726", "PD833", "PD936",
                "PDC167", "PDC184", "PDC197", "PDC87"),
  condition = c("PD", "PD", "PD", "PD", "PD",
                "Control", "Control", "Control", "Control"),
  stringsAsFactors = FALSE
)

# SELECT TOP 5000 HVFs PER VIEW
# bulk_counts is used here as the bulk input matrix for MOFA
bulk_vars    <- apply(bulk_counts, 1, var)
bulk_top_idx <- head(order(bulk_vars, decreasing = TRUE), 5000)
bulk_hvf     <- bulk_counts[bulk_top_idx, ]

opalin_vars    <- apply(opalin_vst, 1, var)
opalin_top_idx <- head(order(opalin_vars, decreasing = TRUE), 5000)
opalin_hvf     <- opalin_vst[opalin_top_idx, ]

plekhg1_vars    <- apply(plekhg1_vst, 1, var)
plekhg1_top_idx <- head(order(plekhg1_vars, decreasing = TRUE), 5000)
plekhg1_hvf     <- plekhg1_vst[plekhg1_top_idx, ]

opc_vars    <- apply(opc_vst, 1, var)
opc_top_idx <- head(order(opc_vars, decreasing = TRUE), 5000)
opc_hvf     <- opc_vst[opc_top_idx, ]

# CHECK SAMPLE ORDER CONSISTENCY ACROSS VIEWS
cat("Bulk vs Opalin sample order identical: ", identical(colnames(bulk_hvf), colnames(opalin_hvf)), "\n")
cat("Bulk vs Plekhg1 sample order identical: ", identical(colnames(bulk_hvf), colnames(plekhg1_hvf)), "\n")
cat("Bulk vs OPC sample order identical: ", identical(colnames(bulk_hvf), colnames(opc_hvf)), "\n")

# ALIGN SAMPLE METADATA TO BULK_HVF COLUMN ORDER
sample_meta <- sample_meta[match(colnames(bulk_hvf), sample_meta$sample), ]

if (any(is.na(sample_meta$sample))) {
  stop("Sample metadata alignment failed: at least one bulk_hvf sample was not found in sample_meta.")
}

cat("Sample metadata aligned to bulk_hvf columns\n")
print(sample_meta)

# VERIFY HVF PEAKS MATCH ANNOTATION
cat("Bulk HVF in annotation:", sum(rownames(bulk_hvf) %in% bulk_anno_df$peak_key), "/ 5000\n")
cat("Opalin HVF in annotation:", sum(rownames(opalin_hvf) %in% opalin_anno_cp$peak_key), "/ 5000\n")
cat("Plekhg1 HVF in annotation:", sum(rownames(plekhg1_hvf) %in% plekhg1_anno_cp$peak_key), "/ 5000\n")
cat("OPC HVF in annotation:", sum(rownames(opc_hvf) %in% opc_anno_cp$peak_key), "/ 5000\n")

# HELPER: TRAIN MOFA
train_mofa <- function(input_list, outfile, n_factors = 7, seed = 42) {
  mofa_obj <- create_mofa(input_list)
  
  model_opts <- get_default_model_options(mofa_obj)
  model_opts$num_factors <- n_factors
  
  train_opts <- get_default_training_options(mofa_obj)
  train_opts$convergence_mode <- "slow"
  train_opts$seed <- seed
  
  mofa_obj <- prepare_mofa(
    mofa_obj,
    data_options     = get_default_data_options(mofa_obj),
    model_options    = model_opts,
    training_options = train_opts
  )
  
  mofa_obj <- run_mofa(
    mofa_obj,
    outfile      = file.path(mofa_out_dir, outfile),
    use_basilisk = TRUE
  )
  
  samples_metadata(mofa_obj) <- sample_meta
  return(mofa_obj)
}

# EXTRACT WEIGHTS
get_w <- function(mofa_obj, factor_num, view) {
  w <- get_weights(mofa_obj, views = view, factors = factor_num, as.data.frame = TRUE)
  w$feature <- as.character(w$feature)
  w
}

# SAVE MOFA PLOTS + VARIANCE TABLES
# Uses png() instead of ggsave()
save_mofa_outputs <- function(mofa_obj, prefix) {
  png(
    filename = file.path(mofa_out_dir, paste0(prefix, "_variance_explained.png")),
    width = 3000, height = 2400, res = 300
  )
  print(plot_variance_explained(mofa_obj, max_r2 = 15))
  dev.off()
  
  png(
    filename = file.path(mofa_out_dir, paste0(prefix, "_factor_cor.png")),
    width = 3000, height = 2400, res = 300
  )
  print(plot_factor_cor(mofa_obj))
  dev.off()
  
  png(
    filename = file.path(mofa_out_dir, paste0(prefix, "_factors.png")),
    width = 3000, height = 2400, res = 300
  )
  print(plot_factor(
    mofa_obj,
    factors = 1:min(7, get_dimensions(mofa_obj)$K),
    color_by = "condition"
  ))
  dev.off()
  
  ve <- get_variance_explained(mofa_obj)
  
  write.csv(
    ve$r2_per_factor,
    file.path(mofa_out_dir, paste0(prefix, "_variance_explained_by_factor.csv")),
    row.names = FALSE
  )
  
  if (!is.null(ve$r2_total)) {
    write.csv(
      ve$r2_total,
      file.path(mofa_out_dir, paste0(prefix, "_variance_explained_total.csv")),
      row.names = FALSE
    )
  }
}

# HELPER: SAVE TOP FEATURES
get_top_features <- function(weights_df, anno_df, label, mode = c("abs", "positive", "negative"), n_top = 20) {
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
  
  if (nrow(top_peaks) == 0) {
    cat("No peaks found for", label, "-", mode, "\n")
    return(NULL)
  }
  
  top_peaks <- head(top_peaks, n_top)
  
  annotated <- merge(
    top_peaks,
    anno_df[, c("peak_key", "annotation", "SYMBOL", "GENENAME", "distanceToTSS")],
    by.x = "feature",
    by.y = "peak_key",
    all.x = TRUE
  )
  
  annotated <- annotated[!duplicated(annotated$feature), ]
  
  write.csv(
    annotated,
    file.path(mofa_out_dir, paste0("top_", n_top, "_", label, "_", mode, ".csv")),
    row.names = FALSE
  )
  
  cat("Saved:", label, "-", mode, "- n =", nrow(annotated), "\n")
  annotated
}

# HELPER: RUN ALL THREE MODES
save_all_feature_modes <- function(weights_df, anno_df, label, n_top = 20) {
  get_top_features(weights_df, anno_df, label, mode = "abs", n_top = n_top)
  get_top_features(weights_df, anno_df, label, mode = "positive", n_top = n_top)
  get_top_features(weights_df, anno_df, label, mode = "negative", n_top = n_top)
}

# MODEL 1: BULK VS OPALIN+

mofa_opalin_vst <- train_mofa(
  list(bulk = as.matrix(bulk_hvf), opalin = as.matrix(opalin_hvf)),
  "mofa_opalin_newsc.hdf5",
  n_factors = 7
)

save_mofa_outputs(mofa_opalin_vst, "opalin")

# Extract weights for all 7 factors
bulk_weights_f1   <- get_w(mofa_opalin_vst, 1, "bulk")
opalin_weights_f1 <- get_w(mofa_opalin_vst, 1, "opalin")

bulk_weights_f2   <- get_w(mofa_opalin_vst, 2, "bulk")
opalin_weights_f2 <- get_w(mofa_opalin_vst, 2, "opalin")

bulk_weights_f3   <- get_w(mofa_opalin_vst, 3, "bulk")
opalin_weights_f3 <- get_w(mofa_opalin_vst, 3, "opalin")

bulk_weights_f4   <- get_w(mofa_opalin_vst, 4, "bulk")
opalin_weights_f4 <- get_w(mofa_opalin_vst, 4, "opalin")

bulk_weights_f5   <- get_w(mofa_opalin_vst, 5, "bulk")
opalin_weights_f5 <- get_w(mofa_opalin_vst, 5, "opalin")

bulk_weights_f6   <- get_w(mofa_opalin_vst, 6, "bulk")
opalin_weights_f6 <- get_w(mofa_opalin_vst, 6, "opalin")

bulk_weights_f7   <- get_w(mofa_opalin_vst, 7, "bulk")
opalin_weights_f7 <- get_w(mofa_opalin_vst, 7, "opalin")

# Opalin model extraction plan
save_all_feature_modes(bulk_weights_f1, bulk_anno_df, "opalin_factor1_bulk")
save_all_feature_modes(opalin_weights_f2, opalin_anno_cp, "opalin_factor2_opalin")
save_all_feature_modes(opalin_weights_f3, opalin_anno_cp, "opalin_factor3_opalin")
save_all_feature_modes(bulk_weights_f4, bulk_anno_df, "opalin_factor4_bulk")
save_all_feature_modes(opalin_weights_f4, opalin_anno_cp, "opalin_factor4_opalin")
save_all_feature_modes(bulk_weights_f5, bulk_anno_df, "opalin_factor5_bulk")
save_all_feature_modes(opalin_weights_f6, opalin_anno_cp, "opalin_factor6_opalin")
save_all_feature_modes(bulk_weights_f7, bulk_anno_df, "opalin_factor7_bulk")
save_all_feature_modes(opalin_weights_f7, opalin_anno_cp, "opalin_factor7_opalin")

saveRDS(mofa_opalin_vst, file.path(mofa_out_dir, "mofa_opalin_vst_object.rds"))
cat("Opalin model complete\n")


# MODEL 2: BULK VS PLEKHG1+

mofa_plekhg1_vst <- train_mofa(
  list(bulk = as.matrix(bulk_hvf), plekhg1 = as.matrix(plekhg1_hvf)),
  "mofa_plekhg1_newsc.hdf5",
  n_factors = 7
)

save_mofa_outputs(mofa_plekhg1_vst, "plekhg1")

# Extract weights for all 7 factors
bulk_weights_plekhg1_f1 <- get_w(mofa_plekhg1_vst, 1, "bulk")
plekhg1_weights_f1      <- get_w(mofa_plekhg1_vst, 1, "plekhg1")

bulk_weights_plekhg1_f2 <- get_w(mofa_plekhg1_vst, 2, "bulk")
plekhg1_weights_f2      <- get_w(mofa_plekhg1_vst, 2, "plekhg1")

bulk_weights_plekhg1_f3 <- get_w(mofa_plekhg1_vst, 3, "bulk")
plekhg1_weights_f3      <- get_w(mofa_plekhg1_vst, 3, "plekhg1")

bulk_weights_plekhg1_f4 <- get_w(mofa_plekhg1_vst, 4, "bulk")
plekhg1_weights_f4      <- get_w(mofa_plekhg1_vst, 4, "plekhg1")

bulk_weights_plekhg1_f5 <- get_w(mofa_plekhg1_vst, 5, "bulk")
plekhg1_weights_f5      <- get_w(mofa_plekhg1_vst, 5, "plekhg1")

bulk_weights_plekhg1_f6 <- get_w(mofa_plekhg1_vst, 6, "bulk")
plekhg1_weights_f6      <- get_w(mofa_plekhg1_vst, 6, "plekhg1")

bulk_weights_plekhg1_f7 <- get_w(mofa_plekhg1_vst, 7, "bulk")
plekhg1_weights_f7      <- get_w(mofa_plekhg1_vst, 7, "plekhg1")

# Plekhg1 model extraction plan
save_all_feature_modes(plekhg1_weights_f1, plekhg1_anno_cp, "plekhg1_factor1_plekhg1")
save_all_feature_modes(bulk_weights_plekhg1_f2, bulk_anno_df, "plekhg1_factor2_bulk")
save_all_feature_modes(bulk_weights_plekhg1_f3, bulk_anno_df, "plekhg1_factor3_bulk")
save_all_feature_modes(bulk_weights_plekhg1_f4, bulk_anno_df, "plekhg1_factor4_bulk")
save_all_feature_modes(plekhg1_weights_f5, plekhg1_anno_cp, "plekhg1_factor5_plekhg1")
save_all_feature_modes(bulk_weights_plekhg1_f6, bulk_anno_df, "plekhg1_factor6_bulk")
save_all_feature_modes(plekhg1_weights_f6, plekhg1_anno_cp, "plekhg1_factor6_plekhg1")
save_all_feature_modes(bulk_weights_plekhg1_f7, bulk_anno_df, "plekhg1_factor7_bulk")
save_all_feature_modes(plekhg1_weights_f7, plekhg1_anno_cp, "plekhg1_factor7_plekhg1")

saveRDS(mofa_plekhg1_vst, file.path(mofa_out_dir, "mofa_plekhg1_vst_object.rds"))
cat("Plekhg1 model complete\n")


# MODEL 3: BULK VS OPC

mofa_opc_vst <- train_mofa(
  list(bulk = as.matrix(bulk_hvf), opc = as.matrix(opc_hvf)),
  "mofa_opc_newsc.hdf5",
  n_factors = 7
)

save_mofa_outputs(mofa_opc_vst, "opc")

# Extract weights for all 7 factors
bulk_weights_opc_f1 <- get_w(mofa_opc_vst, 1, "bulk")
opc_weights_f1      <- get_w(mofa_opc_vst, 1, "opc")

bulk_weights_opc_f2 <- get_w(mofa_opc_vst, 2, "bulk")
opc_weights_f2      <- get_w(mofa_opc_vst, 2, "opc")

bulk_weights_opc_f3 <- get_w(mofa_opc_vst, 3, "bulk")
opc_weights_f3      <- get_w(mofa_opc_vst, 3, "opc")

bulk_weights_opc_f4 <- get_w(mofa_opc_vst, 4, "bulk")
opc_weights_f4      <- get_w(mofa_opc_vst, 4, "opc")

bulk_weights_opc_f5 <- get_w(mofa_opc_vst, 5, "bulk")
opc_weights_f5      <- get_w(mofa_opc_vst, 5, "opc")

bulk_weights_opc_f6 <- get_w(mofa_opc_vst, 6, "bulk")
opc_weights_f6      <- get_w(mofa_opc_vst, 6, "opc")

bulk_weights_opc_f7 <- get_w(mofa_opc_vst, 7, "bulk")
opc_weights_f7      <- get_w(mofa_opc_vst, 7, "opc")

# OPC model extraction plan
save_all_feature_modes(bulk_weights_opc_f1, bulk_anno_df, "opc_factor1_bulk")
save_all_feature_modes(opc_weights_f2, opc_anno_cp, "opc_factor2_opc")
save_all_feature_modes(opc_weights_f3, opc_anno_cp, "opc_factor3_opc")
save_all_feature_modes(bulk_weights_opc_f4, bulk_anno_df, "opc_factor4_bulk")
save_all_feature_modes(bulk_weights_opc_f5, bulk_anno_df, "opc_factor5_bulk")
save_all_feature_modes(bulk_weights_opc_f6, bulk_anno_df, "opc_factor6_bulk")
save_all_feature_modes(opc_weights_f6, opc_anno_cp, "opc_factor6_opc")
save_all_feature_modes(bulk_weights_opc_f7, bulk_anno_df, "opc_factor7_bulk")
save_all_feature_modes(opc_weights_f7, opc_anno_cp, "opc_factor7_opc")

saveRDS(mofa_opc_vst, file.path(mofa_out_dir, "mofa_opc_vst_object.rds"))
cat("OPC model complete\n")


# MODEL 4: 4-VIEW (BULK + ALL THREE SUBTYPES)

mofa_4view_vst <- train_mofa(
  list(
    bulk    = as.matrix(bulk_hvf),
    opalin  = as.matrix(opalin_hvf),
    plekhg1 = as.matrix(plekhg1_hvf),
    opc     = as.matrix(opc_hvf)
  ),
  "mofa_4view_newsc.hdf5",
  n_factors = 7
)

save_mofa_outputs(mofa_4view_vst, "4view")

# Extract weights for all 7 factors across all 4 views
bulk_weights_4v_f1    <- get_w(mofa_4view_vst, 1, "bulk")
opalin_weights_4v_f1  <- get_w(mofa_4view_vst, 1, "opalin")
plekhg1_weights_4v_f1 <- get_w(mofa_4view_vst, 1, "plekhg1")
opc_weights_4v_f1     <- get_w(mofa_4view_vst, 1, "opc")

bulk_weights_4v_f2    <- get_w(mofa_4view_vst, 2, "bulk")
opalin_weights_4v_f2  <- get_w(mofa_4view_vst, 2, "opalin")
plekhg1_weights_4v_f2 <- get_w(mofa_4view_vst, 2, "plekhg1")
opc_weights_4v_f2     <- get_w(mofa_4view_vst, 2, "opc")

bulk_weights_4v_f3    <- get_w(mofa_4view_vst, 3, "bulk")
opalin_weights_4v_f3  <- get_w(mofa_4view_vst, 3, "opalin")
plekhg1_weights_4v_f3 <- get_w(mofa_4view_vst, 3, "plekhg1")
opc_weights_4v_f3     <- get_w(mofa_4view_vst, 3, "opc")

bulk_weights_4v_f4    <- get_w(mofa_4view_vst, 4, "bulk")
opalin_weights_4v_f4  <- get_w(mofa_4view_vst, 4, "opalin")
plekhg1_weights_4v_f4 <- get_w(mofa_4view_vst, 4, "plekhg1")
opc_weights_4v_f4     <- get_w(mofa_4view_vst, 4, "opc")

bulk_weights_4v_f5    <- get_w(mofa_4view_vst, 5, "bulk")
opalin_weights_4v_f5  <- get_w(mofa_4view_vst, 5, "opalin")
plekhg1_weights_4v_f5 <- get_w(mofa_4view_vst, 5, "plekhg1")
opc_weights_4v_f5     <- get_w(mofa_4view_vst, 5, "opc")

bulk_weights_4v_f6    <- get_w(mofa_4view_vst, 6, "bulk")
opalin_weights_4v_f6  <- get_w(mofa_4view_vst, 6, "opalin")
plekhg1_weights_4v_f6 <- get_w(mofa_4view_vst, 6, "plekhg1")
opc_weights_4v_f6     <- get_w(mofa_4view_vst, 6, "opc")

bulk_weights_4v_f7    <- get_w(mofa_4view_vst, 7, "bulk")
opalin_weights_4v_f7  <- get_w(mofa_4view_vst, 7, "opalin")
plekhg1_weights_4v_f7 <- get_w(mofa_4view_vst, 7, "plekhg1")
opc_weights_4v_f7     <- get_w(mofa_4view_vst, 7, "opc")

# 4-view extraction plan
save_all_feature_modes(opalin_weights_4v_f1, opalin_anno_cp, "4view_factor1_opalin")
save_all_feature_modes(plekhg1_weights_4v_f1, plekhg1_anno_cp, "4view_factor1_plekhg1")
save_all_feature_modes(opc_weights_4v_f1, opc_anno_cp, "4view_factor1_opc")

save_all_feature_modes(bulk_weights_4v_f2, bulk_anno_df, "4view_factor2_bulk")

save_all_feature_modes(opalin_weights_4v_f3, opalin_anno_cp, "4view_factor3_opalin")
save_all_feature_modes(opc_weights_4v_f3, opc_anno_cp, "4view_factor3_opc")
save_all_feature_modes(plekhg1_weights_4v_f3, plekhg1_anno_cp, "4view_factor3_plekhg1")

save_all_feature_modes(bulk_weights_4v_f4, bulk_anno_df, "4view_factor4_bulk")
save_all_feature_modes(opalin_weights_4v_f4, opalin_anno_cp, "4view_factor4_opalin")

save_all_feature_modes(bulk_weights_4v_f5, bulk_anno_df, "4view_factor5_bulk")
save_all_feature_modes(opalin_weights_4v_f5, opalin_anno_cp, "4view_factor5_opalin")
save_all_feature_modes(plekhg1_weights_4v_f5, plekhg1_anno_cp, "4view_factor5_plekhg1")
save_all_feature_modes(opc_weights_4v_f5, opc_anno_cp, "4view_factor5_opc")

save_all_feature_modes(bulk_weights_4v_f6, bulk_anno_df, "4view_factor6_bulk")
save_all_feature_modes(opalin_weights_4v_f6, opalin_anno_cp, "4view_factor6_opalin")
save_all_feature_modes(plekhg1_weights_4v_f6, plekhg1_anno_cp, "4view_factor6_plekhg1")
save_all_feature_modes(opc_weights_4v_f6, opc_anno_cp, "4view_factor6_opc")

save_all_feature_modes(bulk_weights_4v_f7, bulk_anno_df, "4view_factor7_bulk")
save_all_feature_modes(opalin_weights_4v_f7, opalin_anno_cp, "4view_factor7_opalin")
save_all_feature_modes(plekhg1_weights_4v_f7, plekhg1_anno_cp, "4view_factor7_plekhg1")
save_all_feature_modes(opc_weights_4v_f7, opc_anno_cp, "4view_factor7_opc")

saveRDS(mofa_4view_vst, file.path(mofa_out_dir, "mofa_4view_vst_object.rds"))

# SAVE ALL WEIGHTS FOR rGREAT SCRIPT
save(
  bulk_weights_f1, opalin_weights_f1,
  bulk_weights_f2, opalin_weights_f2,
  bulk_weights_f3, opalin_weights_f3,
  bulk_weights_f4, opalin_weights_f4,
  bulk_weights_f5, opalin_weights_f5,
  bulk_weights_f6, opalin_weights_f6,
  bulk_weights_f7, opalin_weights_f7,
  
  bulk_weights_plekhg1_f1, plekhg1_weights_f1,
  bulk_weights_plekhg1_f2, plekhg1_weights_f2,
  bulk_weights_plekhg1_f3, plekhg1_weights_f3,
  bulk_weights_plekhg1_f4, plekhg1_weights_f4,
  bulk_weights_plekhg1_f5, plekhg1_weights_f5,
  bulk_weights_plekhg1_f6, plekhg1_weights_f6,
  bulk_weights_plekhg1_f7, plekhg1_weights_f7,
  
  bulk_weights_opc_f1, opc_weights_f1,
  bulk_weights_opc_f2, opc_weights_f2,
  bulk_weights_opc_f3, opc_weights_f3,
  bulk_weights_opc_f4, opc_weights_f4,
  bulk_weights_opc_f5, opc_weights_f5,
  bulk_weights_opc_f6, opc_weights_f6,
  bulk_weights_opc_f7, opc_weights_f7,
  
  bulk_weights_4v_f1, opalin_weights_4v_f1, plekhg1_weights_4v_f1, opc_weights_4v_f1,
  bulk_weights_4v_f2, opalin_weights_4v_f2, plekhg1_weights_4v_f2, opc_weights_4v_f2,
  bulk_weights_4v_f3, opalin_weights_4v_f3, plekhg1_weights_4v_f3, opc_weights_4v_f3,
  bulk_weights_4v_f4, opalin_weights_4v_f4, plekhg1_weights_4v_f4, opc_weights_4v_f4,
  bulk_weights_4v_f5, opalin_weights_4v_f5, plekhg1_weights_4v_f5, opc_weights_4v_f5,
  bulk_weights_4v_f6, opalin_weights_4v_f6, plekhg1_weights_4v_f6, opc_weights_4v_f6,
  bulk_weights_4v_f7, opalin_weights_4v_f7, plekhg1_weights_4v_f7, opc_weights_4v_f7,
  
  file = file.path(mofa_out_dir, "all_mofa_weights.RData")
)

cat("All models complete - weights saved for rGREAT script\n")
cat("All results saved to:", mofa_out_dir, "\n")



