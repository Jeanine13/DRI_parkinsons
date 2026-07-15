# LOAD MODEL IF NOT IN MEMORY
mofa_opalin_vst <- readRDS(file.path(mofa_out_dir, "mofa_opalin_vst_object.rds"))

# EXTRACT FACTOR SCORES AND MERGE WITH METADATA
factor_scores_opalin <- as.data.frame(get_factors(mofa_opalin_vst, as.data.frame = TRUE))
factor_scores_opalin <- merge(factor_scores_opalin, sample_meta_full,
                              by.x = "sample", by.y = "sample")

# VIEW FACTOR 1 SCORES ORDERED HIGH TO LOW
factor1 <- factor_scores_opalin[factor_scores_opalin$factor == "Factor1",
                                c("sample", "value", "group")]
factor1 <- factor1[order(factor1$value, decreasing = TRUE), ]
print(factor1)


factor1$condition <- ifelse(factor1$sample %in% 
                              c("PD1222", "PD1231", "PD726", "PD833", "PD936"),
                            "PD", "Control")

print(factor1[, c("sample", "value", "condition")])

#updated sample_meta_full to include proper 'conditon' column
sample_meta_full$condition <- ifelse(sample_meta_full$pd_binary == 1, "PD", "Control")



#UPDATED DUCNTION 


# CREATE FACTOR SCORES DIRECTORY
factor_scores_dir <- file.path(mofa_out_dir, "factor_scores")
dir.create(factor_scores_dir, recursive = TRUE, showWarnings = FALSE)

# FUNCTION TO EXTRACT AND SAVE FACTOR SCORES FOR A MODEL
save_factor_scores <- function(mofa_obj, model_name) {
  
  factor_scores <- as.data.frame(get_factors(mofa_obj, as.data.frame = TRUE))
  
  # ADD CONDITION MANUALLY
  factor_scores$condition <- ifelse(factor_scores$sample %in%
                                      c("PD1222", "PD1231", "PD726", "PD833", "PD936"),
                                    "PD", "Control")
  
  # ORDER BY FACTOR THEN BY VALUE
  factor_scores <- factor_scores[order(factor_scores$factor,
                                       factor_scores$value,
                                       decreasing = TRUE), ]
  
  # SAVE TO CSV
  write.csv(factor_scores,
            file.path(factor_scores_dir, paste0(model_name, "_factor_scores.csv")),
            row.names = FALSE)
  
  cat("Saved factor scores for", model_name, "\n")
  return(factor_scores)
}

# RUN FOR ALL FOUR MODELS
scores_opalin  <- save_factor_scores(mofa_opalin_vst,  "opalin")
scores_plekhg1 <- save_factor_scores(mofa_plekhg1_vst, "plekhg1")
scores_opc     <- save_factor_scores(mofa_opc_vst,     "opc")
scores_4view   <- save_factor_scores(mofa_4view_vst,   "4view")

cat("All factor scores saved to:", factor_scores_dir, "\n")





scores_opalin[scores_opalin$factor == "Factor2", 
              c("sample", "value", "condition")]












