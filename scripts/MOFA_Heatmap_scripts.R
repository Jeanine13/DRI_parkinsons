.libPaths(c(
  "/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5",
  .libPaths()
))
library(MOFA2)


# TOP-FEATURES DATA HEATMAPS (4-VIEW MODEL)

# Standalone script — loads the already-trained 4-view MOFA object directly

# steps. plot_data_heatmap() selects its own top features and pulls sample
# metadata (condition) straight from the saved object

# mofa_out_dir = where the trained model (.rds) already lives — do not change
mofa_out_dir <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/new_sc_mofa_clean_14042026"

# heatmap_dir = separate subfolder where the new heatmap PDFs get written
heatmap_dir <- file.path(mofa_out_dir, "mofa_heatmaps")
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

mofa_4view_vst <- readRDS(file.path(mofa_out_dir, "mofa_4view_vst_object.rds"))
cat("4-view MOFA object loaded\n")

save_data_heatmap <- function(mofa_obj, view, factor_num, prefix, n_features = 20) {
  pdf(
    file = file.path(heatmap_dir, paste0(prefix, "_factor", factor_num, "_", view, "_heatmap.pdf")),
    width = 8, height = 10
  )
  print(plot_data_heatmap(
    mofa_obj,
    view                = view,
    factor              = factor_num,
    features            = n_features,
    annotation_samples  = "condition",
    cluster_cols        = TRUE,
    cluster_rows        = TRUE,
    show_rownames       = TRUE,
    show_colnames       = TRUE,
    fontsize_row        = 6
  ))
  dev.off()
  cat("Saved heatmap:", prefix, "- factor", factor_num, "-", view, "\n")
}

# Mirrors the 4-view feature-extraction plan from the main script:
heatmap_plan_4view <- list(
  list(factor = 1, views = c("opalin", "plekhg1", "opc")),
  list(factor = 2, views = "bulk"),
  list(factor = 3, views = c("opalin", "opc", "plekhg1")),
  list(factor = 4, views = c("bulk", "opalin")),
  list(factor = 5, views = c("bulk", "opalin", "plekhg1", "opc")),
  list(factor = 6, views = c("bulk", "opalin", "plekhg1", "opc")),
  list(factor = 7, views = c("bulk", "opalin", "plekhg1", "opc"))
)

for (plan in heatmap_plan_4view) {
  for (v in plan$views) {
    save_data_heatmap(
      mofa_4view_vst,
      view       = v,
      factor_num = plan$factor,
      prefix     = "4view",
      n_features = 20
    )
  }
}

cat("All 4-view heatmaps saved to:", heatmap_dir, "\n")
