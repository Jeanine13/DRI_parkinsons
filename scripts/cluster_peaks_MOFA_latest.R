# MOFA TRAINING AND WEIGHT EXTRACTION SCRIPT
# RUN AFTER cluster_peaks_chipseeker.R
# TRAINS ALL FOUR MOFA MODELS AND EXTRACTS TOP 20 ANNOTATED PEAKS PER FACTOR

.libPaths(c(
  "/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5",
  .libPaths()
))

library(MOFA2)
library(ggplot2)

# LOAD SAVED CHIPSEEKER ANNOTATIONS
mofa_out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/new_sc_mofa_clean_14042026"

load(file.path(mofa_out_dir, "chipseeker_annotations_all_datasets.RData"))
load("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/mofa_OLIGO_clusterpeaks_granges_vst.RData")

cat("All required objects loaded\n")

# SAMPLE METADATA FOR PD ASSOCIATION TESTING
# stringsAsFactors = FALSE PREVENTS FACTOR CONVERSION ISSUES
# ALIGNED TO BULK_HVF COLUMN ORDER TO ENSURE CORRECT SAMPLE MATCHING
sample_meta <- data.frame(
  sample    = c("PD1222", "PD1231", "PD726", "PD833", "PD936",
                "PDC167", "PDC184", "PDC197", "PDC87"),
  condition = c("PD", "PD", "PD", "PD", "PD",
                "Control", "Control", "Control", "Control"),
  stringsAsFactors = FALSE
)

# SELECT TOP 5000 HVFs PER VIEW
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

# ALIGN SAMPLE METADATA TO BULK_HVF COLUMN ORDER
sample_meta <- sample_meta[match(colnames(bulk_hvf), sample_meta$sample), ]
cat("Sample metadata aligned to bulk_hvf columns\n")
print(sample_meta)

# VERIFY HVF PEAKS MATCH ANNOTATION
cat("Bulk HVF in annotation:", sum(rownames(bulk_hvf) %in% bulk_anno_df$peak_key), "/ 5000\n")
cat("Opalin HVF in annotation:", sum(rownames(opalin_hvf) %in% opalin_anno_cp$peak_key), "/ 5000\n")
cat("Plekhg1 HVF in annotation:", sum(rownames(plekhg1_hvf) %in% plekhg1_anno_cp$peak_key), "/ 5000\n")
cat("OPC HVF in annotation:", sum(rownames(opc_hvf) %in% opc_anno_cp$peak_key), "/ 5000\n")

# MOFA TRAINING HELPER FUNCTION
train_mofa <- function(input_list, outfile, n_factors = 9, seed = 42) {
  mofa_obj <- create_mofa(input_list)
  model_opts             <- get_default_model_options(mofa_obj)
  model_opts$num_factors <- n_factors
  train_opts             <- get_default_training_options(mofa_obj)
  train_opts$convergence_mode <- "slow"
  train_opts$seed        <- seed
  mofa_obj <- prepare_mofa(mofa_obj,
                           data_options     = get_default_data_options(mofa_obj),
                           model_options    = model_opts,
                           training_options = train_opts)
  mofa_obj <- run_mofa(mofa_obj,
                       outfile      = file.path(mofa_out_dir, outfile),
                       use_basilisk = TRUE)
  samples_metadata(mofa_obj) <- sample_meta
  return(mofa_obj)
}

# TOP 20 ANNOTATION FUNCTION
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
  cat("Annotated peaks returned:", nrow(annotated), "\n")
  annotated
}

# HELPER TO EXTRACT WEIGHTS
get_w <- function(mofa_obj, factor_num, view) {
  w <- get_weights(mofa_obj, views = view, factors = factor_num, as.data.frame = TRUE)
  w$feature <- as.character(w$feature)
  w
}

# MODEL 1: BULK VS OPALIN+
mofa_opalin_vst <- train_mofa(
  list(bulk = as.matrix(bulk_hvf), opalin = as.matrix(opalin_hvf)),
  "mofa_opalin_newsc.hdf5"
)

pdf(file.path(mofa_out_dir, "opalin_variance_explained.pdf"))
print(plot_variance_explained(mofa_opalin_vst, max_r2 = 15))
dev.off()

pdf(file.path(mofa_out_dir, "opalin_factor_cor.pdf"))
print(plot_factor_cor(mofa_opalin_vst))
dev.off()

pdf(file.path(mofa_out_dir, "opalin_factors.pdf"))
print(plot_factor(mofa_opalin_vst, factors = 1:get_dimensions(mofa_opalin_vst)$K, color_by = "condition"))
dev.off()

# EXTRACT OPALIN WEIGHTS
# FACTOR 1 BULK, FACTORS 2+3 OPALIN, FACTOR 4 SHARED, FACTOR 5 BULK, FACTOR 6 OPALIN
bulk_weights_f1   <- get_w(mofa_opalin_vst, 1, "bulk")
opalin_weights_f2 <- get_w(mofa_opalin_vst, 2, "opalin")
opalin_weights_f3 <- get_w(mofa_opalin_vst, 3, "opalin")
bulk_weights_f4   <- get_w(mofa_opalin_vst, 4, "bulk")
opalin_weights_f4 <- get_w(mofa_opalin_vst, 4, "opalin")
bulk_weights_f5   <- get_w(mofa_opalin_vst, 5, "bulk")
opalin_weights_f6 <- get_w(mofa_opalin_vst, 6, "opalin")

# TOP 20 - OPALIN MODEL
get_top20(bulk_weights_f1,   bulk_anno_df,   "peak_key", "GENENAME", "opalin_factor1_bulk")
get_top20(opalin_weights_f2, opalin_anno_cp, "peak_key", "GENENAME", "opalin_factor2_opalin")
get_top20(opalin_weights_f3, opalin_anno_cp, "peak_key", "GENENAME", "opalin_factor3_opalin")
get_top20(bulk_weights_f4,   bulk_anno_df,   "peak_key", "GENENAME", "opalin_factor4_bulk")
get_top20(opalin_weights_f4, opalin_anno_cp, "peak_key", "GENENAME", "opalin_factor4_opalin")
get_top20(bulk_weights_f5,   bulk_anno_df,   "peak_key", "GENENAME", "opalin_factor5_bulk")
get_top20(opalin_weights_f6, opalin_anno_cp, "peak_key", "GENENAME", "opalin_factor6_opalin")

saveRDS(mofa_opalin_vst, file.path(mofa_out_dir, "mofa_opalin_vst_object.rds"))
cat("Opalin model complete\n")

# MODEL 2: BULK VS PLEKHG1+
mofa_plekhg1_vst <- train_mofa(
  list(bulk = as.matrix(bulk_hvf), plekhg1 = as.matrix(plekhg1_hvf)),
  "mofa_plekhg1_newsc.hdf5"
)

pdf(file.path(mofa_out_dir, "plekhg1_variance_explained.pdf"))
print(plot_variance_explained(mofa_plekhg1_vst, max_r2 = 15))
dev.off()

pdf(file.path(mofa_out_dir, "plekhg1_factor_cor.pdf"))
print(plot_factor_cor(mofa_plekhg1_vst))
dev.off()

pdf(file.path(mofa_out_dir, "plekhg1_factors.pdf"))
print(plot_factor(mofa_plekhg1_vst, factors = 1:get_dimensions(mofa_plekhg1_vst)$K, color_by = "condition"))
dev.off()

# EXTRACT PLEKHG1 WEIGHTS
# FACTOR 1 PLEKHG1, FACTORS 2+3+4 BULK, FACTOR 5 PLEKHG1
plekhg1_weights_f1      <- get_w(mofa_plekhg1_vst, 1, "plekhg1")
bulk_weights_plekhg1_f2 <- get_w(mofa_plekhg1_vst, 2, "bulk")
bulk_weights_plekhg1_f3 <- get_w(mofa_plekhg1_vst, 3, "bulk")
bulk_weights_plekhg1_f4 <- get_w(mofa_plekhg1_vst, 4, "bulk")
plekhg1_weights_f5      <- get_w(mofa_plekhg1_vst, 5, "plekhg1")

# TOP 20 - PLEKHG1 MODEL
get_top20(plekhg1_weights_f1,      plekhg1_anno_cp, "peak_key", "GENENAME", "plekhg1_factor1_plekhg1")
get_top20(bulk_weights_plekhg1_f2, bulk_anno_df,    "peak_key", "GENENAME", "plekhg1_factor2_bulk")
get_top20(bulk_weights_plekhg1_f3, bulk_anno_df,    "peak_key", "GENENAME", "plekhg1_factor3_bulk")
get_top20(bulk_weights_plekhg1_f4, bulk_anno_df,    "peak_key", "GENENAME", "plekhg1_factor4_bulk")
get_top20(plekhg1_weights_f5,      plekhg1_anno_cp, "peak_key", "GENENAME", "plekhg1_factor5_plekhg1")

saveRDS(mofa_plekhg1_vst, file.path(mofa_out_dir, "mofa_plekhg1_vst_object.rds"))
cat("Plekhg1 model complete\n")

# MODEL 3: BULK VS OPCs
mofa_opc_vst <- train_mofa(
  list(bulk = as.matrix(bulk_hvf), opc = as.matrix(opc_hvf)),
  "mofa_opc_newsc.hdf5"
)

pdf(file.path(mofa_out_dir, "opc_variance_explained.pdf"))
print(plot_variance_explained(mofa_opc_vst, max_r2 = 15))
dev.off()

pdf(file.path(mofa_out_dir, "opc_factor_cor.pdf"))
print(plot_factor_cor(mofa_opc_vst))
dev.off()

pdf(file.path(mofa_out_dir, "opc_factors.pdf"))
print(plot_factor(mofa_opc_vst, factors = 1:get_dimensions(mofa_opc_vst)$K, color_by = "condition"))
dev.off()

# EXTRACT OPC WEIGHTS
# UPDATE FACTOR/VIEW ASSIGNMENTS AFTER REVIEWING VARIANCE PLOT
bulk_weights_opc_f1 <- get_w(mofa_opc_vst, 1, "bulk")
opc_weights_f2      <- get_w(mofa_opc_vst, 2, "opc")
opc_weights_f3      <- get_w(mofa_opc_vst, 3, "opc")
bulk_weights_opc_f4 <- get_w(mofa_opc_vst, 4, "bulk")
bulk_weights_opc_f5 <- get_w(mofa_opc_vst, 5, "bulk")

# TOP 20 - OPC MODEL
get_top20(bulk_weights_opc_f1, bulk_anno_df, "peak_key", "GENENAME", "opc_factor1_bulk")
get_top20(opc_weights_f2,      opc_anno_cp,  "peak_key", "GENENAME", "opc_factor2_opc")
get_top20(opc_weights_f3,      opc_anno_cp,  "peak_key", "GENENAME", "opc_factor3_opc")
get_top20(bulk_weights_opc_f4, bulk_anno_df, "peak_key", "GENENAME", "opc_factor4_bulk")
get_top20(bulk_weights_opc_f5, bulk_anno_df, "peak_key", "GENENAME", "opc_factor5_bulk")

saveRDS(mofa_opc_vst, file.path(mofa_out_dir, "mofa_opc_vst_object.rds"))
cat("OPC model complete\n")

# MODEL 4: 4-VIEW (BULK + ALL THREE SC SUBTYPES)
mofa_4view_vst <- train_mofa(
  list(bulk    = as.matrix(bulk_hvf),
       opalin  = as.matrix(opalin_hvf),
       plekhg1 = as.matrix(plekhg1_hvf),
       opc     = as.matrix(opc_hvf)),
  "mofa_4view_newsc.hdf5"
)

pdf(file.path(mofa_out_dir, "4view_variance_explained.pdf"))
print(plot_variance_explained(mofa_4view_vst, max_r2 = 15))
dev.off()

pdf(file.path(mofa_out_dir, "4view_factor_cor.pdf"))
print(plot_factor_cor(mofa_4view_vst))
dev.off()

pdf(file.path(mofa_out_dir, "4view_factors.pdf"))
print(plot_factor(mofa_4view_vst, factors = 1:get_dimensions(mofa_4view_vst)$K, color_by = "condition"))
dev.off()

# EXTRACT 4-VIEW WEIGHTS
# UPDATE FACTOR/VIEW ASSIGNMENTS AFTER REVIEWING VARIANCE PLOT
# FACTOR 1 SC SHARED, FACTOR 2 BULK, FACTOR 3 OPALIN+OPC SHARED, FACTOR 4 BULK
opalin_weights_4v_f1  <- get_w(mofa_4view_vst, 1, "opalin")
plekhg1_weights_4v_f1 <- get_w(mofa_4view_vst, 1, "plekhg1")
opc_weights_4v_f1     <- get_w(mofa_4view_vst, 1, "opc")
bulk_weights_4v_f2    <- get_w(mofa_4view_vst, 2, "bulk")
opalin_weights_4v_f3  <- get_w(mofa_4view_vst, 3, "opalin")
opc_weights_4v_f3     <- get_w(mofa_4view_vst, 3, "opc")
bulk_weights_4v_f4    <- get_w(mofa_4view_vst, 4, "bulk")

# TOP 20 - 4-VIEW MODEL
get_top20(opalin_weights_4v_f1,  opalin_anno_cp,  "peak_key", "GENENAME", "4view_factor1_opalin")
get_top20(plekhg1_weights_4v_f1, plekhg1_anno_cp, "peak_key", "GENENAME", "4view_factor1_plekhg1")
get_top20(opc_weights_4v_f1,     opc_anno_cp,     "peak_key", "GENENAME", "4view_factor1_opc")
get_top20(bulk_weights_4v_f2,    bulk_anno_df,    "peak_key", "GENENAME", "4view_factor2_bulk")
get_top20(opalin_weights_4v_f3,  opalin_anno_cp,  "peak_key", "GENENAME", "4view_factor3_opalin")
get_top20(opc_weights_4v_f3,     opc_anno_cp,     "peak_key", "GENENAME", "4view_factor3_opc")
get_top20(bulk_weights_4v_f4,    bulk_anno_df,    "peak_key", "GENENAME", "4view_factor4_bulk")

saveRDS(mofa_4view_vst, file.path(mofa_out_dir, "mofa_4view_vst_object.rds"))

# SAVE ALL WEIGHTS FOR rGREAT SCRIPT
save(bulk_weights_f1, opalin_weights_f2, opalin_weights_f3,
     bulk_weights_f4, opalin_weights_f4, bulk_weights_f5, opalin_weights_f6,
     plekhg1_weights_f1, bulk_weights_plekhg1_f2, bulk_weights_plekhg1_f3,
     bulk_weights_plekhg1_f4, plekhg1_weights_f5,
     bulk_weights_opc_f1, opc_weights_f2, opc_weights_f3,
     bulk_weights_opc_f4, bulk_weights_opc_f5,
     opalin_weights_4v_f1, plekhg1_weights_4v_f1, opc_weights_4v_f1,
     bulk_weights_4v_f2, opalin_weights_4v_f3, opc_weights_4v_f3,
     bulk_weights_4v_f4,
     file = file.path(mofa_out_dir, "all_mofa_weights.RData"))

cat("All models complete - weights saved for rGREAT script\n")
cat("All results saved to:", mofa_out_dir, "\n")





