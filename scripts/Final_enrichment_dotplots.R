# REGENERATE DOTPLOTS FROM SAVED RGREAT CSVS
# TOP 10 TERMS - NARROW FORMAT - LARGER DOTS - CLOSER TERMS
# ALL FOUR MODELS

.libPaths(c(
  "/cephfs/volumes/hpc_data_usr/k25093549/eabe5dc4-1fa9-4cdc-b2af-6a4d37d00142/R/R/x86_64-pc-linux-gnu-library/4.5",
  .libPaths()
))

library(ggplot2)
library(dplyr)
library(viridis)
# DIRECTORIES


mofa_out_dir       <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/new_sc_mofa_clean_14042026"
rgreat_dir         <- file.path(mofa_out_dir, "rgreat_results")
rgreat_dir_4view   <- file.path(mofa_out_dir, "rgreat_results_4view")
dotplot_dir        <- file.path(mofa_out_dir, "dotplots_top10_narrow_2")
dir.create(dotplot_dir, recursive = TRUE, showWarnings = FALSE)


# DOTPLOT FUNCTION
# TOP 10 TERMS - NARROW - LARGER DOTS - CLOSER TERMS


plot_great_dotplot <- function(tb_sig, label, mode, n_show = 10) {
  
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
    scale_colour_viridis_c(
      name      = "p.adjust",
      direction = -1,
      guide     = guide_colourbar(
        barwidth       = 0.4,
        barheight      = 3,
        direction      = "vertical",
        title.position = "top",
        title.hjust    = 0.5
      )
    ) +
    scale_size_continuous(name = "Region hits", range = c(2, 5)) +
    labs(title = paste(label, "-", mode),
         x     = "Fold enrichment",
         y     = NULL) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.y      = element_text(size = 8, lineheight = 0.8),
      axis.text.x      = element_text(size = 8),
      axis.title.x     = element_text(size = 9),
      plot.title       = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7),
      legend.position  = "right",
      plot.margin      = margin(3, 3, 3, 3)
    )
  
  ggsave(
    filename = file.path(dotplot_dir,
                         paste0("dotplot_rGREAT_", label, "_", mode, ".pdf")),
    plot   = p,
    width  = 6,
    height = 3.5,
    units  = "in",
    device = "pdf",
    bg     = "white"
  )
  
  cat("Saved dotplot:", label, "-", mode, "\n")
  return(p)
}


# HELPER FUNCTIONS


plot_from_csv <- function(csv_path, label, mode) {
  if (!file.exists(csv_path)) {
    message("CSV not found - skipping: ", csv_path)
    return(NULL)
  }
  tb <- read.csv(csv_path)
  if (nrow(tb) == 0) {
    message("No significant terms - skipping: ", label, " - ", mode)
    return(NULL)
  }
  plot_great_dotplot(tb, label, mode)
}

plot_both <- function(label, base_dir) {
  plot_from_csv(
    file.path(base_dir, paste0("rGREAT_", label, "_positive.csv")),
    label, "positive"
  )
  plot_from_csv(
    file.path(base_dir, paste0("rGREAT_", label, "_negative.csv")),
    label, "negative"
  )
}


# MODEL 1: OPALIN


cat("\n=== Opalin model ===\n")
plot_both("opalin_factor1_bulk",   rgreat_dir)
plot_both("opalin_factor2_opalin", rgreat_dir)
plot_both("opalin_factor3_opalin", rgreat_dir)
plot_both("opalin_factor4_bulk",   rgreat_dir)
plot_both("opalin_factor4_opalin", rgreat_dir)


# MODEL 2: PLEKHG1


cat("\n=== Plekhg1 model ===\n")
plot_both("plekhg1_factor1_plekhg1", rgreat_dir)
plot_both("plekhg1_factor2_bulk",    rgreat_dir)
plot_both("plekhg1_factor3_bulk",    rgreat_dir)


# MODEL 3: OPC


cat("\n=== OPC model ===\n")
plot_both("opc_factor1_bulk", rgreat_dir)
plot_both("opc_factor2_opc",  rgreat_dir)
plot_both("opc_factor3_opc",  rgreat_dir)


# MODEL 4: 4-VIEW


cat("\n=== 4-view model ===\n")
plot_both("4view_factor1_opalin",  rgreat_dir_4view)
plot_both("4view_factor1_opc",     rgreat_dir_4view)
plot_both("4view_factor2_bulk",    rgreat_dir_4view)
plot_both("4view_factor3_opalin",  rgreat_dir_4view)
plot_both("4view_factor3_opc",     rgreat_dir_4view)

cat("\nAll dotplots saved to:", dotplot_dir, "\n")

