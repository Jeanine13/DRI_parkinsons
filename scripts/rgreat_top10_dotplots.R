# RE-GENERATE TOP 10 DOTPLOTS FROM SAVED RGREAT CSVS
# ALL FOUR MODELS - NO RE-RUNNING RGREAT
# SAVES DOTPLOTS ONLY - NO NEW CSVS

.libPaths(c(
  "/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5",
  .libPaths()
))

library(ggplot2)

# DIRECTORIES
mofa_out_dir   <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/new_sc_mofa_clean_14042026"
rgreat_dir     <- file.path(mofa_out_dir, "rgreat_results")
dotplot_dir    <- file.path(mofa_out_dir, "dotplots_top20")
dir.create(dotplot_dir, recursive = TRUE, showWarnings = FALSE)

# DOTPLOT FUNCTION - TOP 10 TERMS
# BLUE (#0072B2) = MOST SIGNIFICANT, AMBER (#E69F00) = LESS SIGNIFICANT
plot_great_dotplot <- function(tb_sig, label, mode, n_show = 20) {
  
  if (is.null(tb_sig) || nrow(tb_sig) == 0) {
    message("No significant terms for ", label, " - ", mode, " - skipping")
    return(NULL)
  }
  
  plot_df <- head(tb_sig, n_show)
  plot_df$description <- factor(plot_df$description,
                                levels = rev(plot_df$description))
  
  p <- ggplot(
    plot_df,
    aes(x      = fold_enrichment,
        y      = description,
        size   = observed_region_hits,
        colour = p_adjust)
  ) +
    geom_point() +
    scale_colour_gradient(
      low  = "#0072B2",
      high = "#E69F00",
      name = "p.adjust",
      guide = guide_colourbar(
        barwidth       = 0.4,
        barheight      = 3,
        direction      = "vertical",
        title.position = "top",
        title.hjust    = 0.5
      )
    ) +
    scale_size_continuous(name = "Region hits", range = c(1, 4)) +
    labs(title = paste(label, "-", mode),
         x     = "Fold enrichment",
         y     = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.y     = element_text(size = 9),
      axis.text.x     = element_text(size = 9),
      axis.title.x    = element_text(size = 10),
      plot.title      = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title    = element_text(size = 8),
      legend.text     = element_text(size = 7),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(dotplot_dir,
                         paste0("dotplot_rGREAT_", label, "_", mode, ".png")),
    plot   = p,
    width  = 5,
    height = 4,
    units  = "in",
    dpi    = 300,
    bg     = "white"
  )
  
  cat("Saved dotplot:", label, "-", mode, "\n")
  return(p)
}
  

# FUNCTION TO LOAD CSV AND PLOT
plot_from_csv <- function(csv_path, label, mode) {
  
  if (!file.exists(csv_path)) {
    message("CSV not found - skipping: ", csv_path)
    return(NULL)
  }
  
  tb <- read.csv(csv_path)
  
  if (nrow(tb) == 0) {
    message("No significant terms in CSV - skipping: ", label, " - ", mode)
    return(NULL)
  }
  
  plot_great_dotplot(tb, label, mode)
}

# HELPER TO PLOT BOTH DIRECTIONS FOR A GIVEN LABEL
plot_both <- function(label) {
  plot_from_csv(
    file.path(rgreat_dir, paste0("rGREAT_", label, "_positive.csv")),
    label, "positive"
  )
  plot_from_csv(
    file.path(rgreat_dir, paste0("rGREAT_", label, "_negative.csv")),
    label, "negative"
  )
}

# MODEL 1: BULK VS OPALIN+

cat("\n=== Opalin model ===\n")
plot_both("opalin_factor1_bulk")
plot_both("opalin_factor2_opalin")
plot_both("opalin_factor3_opalin")
plot_both("opalin_factor4_bulk")
plot_both("opalin_factor4_opalin")
plot_both("opalin_factor5_bulk")
plot_both("opalin_factor6_opalin")
plot_both("opalin_factor7_bulk")
plot_both("opalin_factor7_opalin")


# MODEL 2: BULK VS PLEKHG1+

cat("\n=== Plekhg1 model ===\n")
plot_both("plekhg1_factor1_plekhg1")
plot_both("plekhg1_factor2_bulk")
plot_both("plekhg1_factor3_bulk")
plot_both("plekhg1_factor4_bulk")
plot_both("plekhg1_factor5_plekhg1")
plot_both("plekhg1_factor6_bulk")
plot_both("plekhg1_factor6_plekhg1")
plot_both("plekhg1_factor7_bulk")
plot_both("plekhg1_factor7_plekhg1")


# MODEL 3: BULK VS OPC

cat("\n=== OPC model ===\n")
plot_both("opc_factor1_bulk")
plot_both("opc_factor2_opc")
plot_both("opc_factor3_opc")
plot_both("opc_factor4_bulk")
plot_both("opc_factor5_bulk")
plot_both("opc_factor6_bulk")
plot_both("opc_factor6_opc")
plot_both("opc_factor7_bulk")
plot_both("opc_factor7_opc")


# MODEL 4: 4-VIEW

cat("\n=== 4-view model ===\n")

# USE 4-VIEW RGREAT DIRECTORY
rgreat_dir <- file.path(mofa_out_dir, "rgreat_results_4view")

plot_both("4view_factor1_opalin")
plot_both("4view_factor1_plekhg1")
plot_both("4view_factor1_opc")
plot_both("4view_factor1_bulk")
plot_both("4view_factor2_bulk")
plot_both("4view_factor2_opalin")
plot_both("4view_factor2_plekhg1")
plot_both("4view_factor2_opc")
plot_both("4view_factor3_opalin")
plot_both("4view_factor3_opc")
plot_both("4view_factor3_plekhg1")
plot_both("4view_factor3_bulk")
plot_both("4view_factor4_bulk")
plot_both("4view_factor4_opalin")
plot_both("4view_factor4_plekhg1")
plot_both("4view_factor4_opc")
plot_both("4view_factor5_bulk")
plot_both("4view_factor5_opalin")
plot_both("4view_factor5_plekhg1")
plot_both("4view_factor5_opc")
plot_both("4view_factor6_bulk")
plot_both("4view_factor6_opalin")
plot_both("4view_factor6_plekhg1")
plot_both("4view_factor6_opc")
plot_both("4view_factor7_bulk")
plot_both("4view_factor7_opalin")
plot_both("4view_factor7_plekhg1")
plot_both("4view_factor7_opc")

cat("\nAll dotplots saved to:", dotplot_dir, "\n")
