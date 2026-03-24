library(MOFA2)

# select top 5000 HVFs per view
bulk_vars    <- apply(bulk_counts, 1, var)
bulk_top_idx <- head(order(bulk_vars, decreasing = TRUE), 5000)
bulk_hvf     <- bulk_counts[bulk_top_idx, ]

opalin_vars    <- apply(opalin_pseudo_norm, 1, var)
opalin_top_idx <- head(order(opalin_vars, decreasing = TRUE), 5000)
opalin_hvf     <- opalin_pseudo_norm[opalin_top_idx, ]

# create 2-view MOFA object (opalin vs bulk)
mofa_input_opalin <- list(
  bulk   = as.matrix(bulk_hvf),
  opalin = as.matrix(opalin_hvf)
)

mofa_opalin_cp <- create_mofa(mofa_input_opalin)
plot_data_overview(mofa_opalin_cp)

# set options
data_opts  <- get_default_data_options(mofa_opalin_cp)
model_opts <- get_default_model_options(mofa_opalin_cp)
model_opts$num_factors <- 9

train_opts <- get_default_training_options(mofa_opalin_cp)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

mofa_opalin_cp <- prepare_mofa(mofa_opalin_cp,
                               data_options     = data_opts,
                               model_options    = model_opts,
                               training_options = train_opts)

mofa_opalin_cp <- run_mofa(mofa_opalin_cp,
                           outfile      = "mofa_opalin_clusterpeaks.hdf5",
                           use_basilisk = TRUE)

sample_meta <- data.frame(
  sample    = c("PD1222", "PD1231", "PD726", "PD833", "PD936",
                "PDC167", "PDC184", "PDC197", "PDC87"),
  condition = c("PD", "PD", "PD", "PD", "PD",
                "Control", "Control", "Control", "Control")
)

samples_metadata(mofa_opalin_cp) <- sample_meta

plot_variance_explained(mofa_opalin_cp, max_r2 = 15)

plot_factor_cor(mofa_opalin_cp)


plot_factor(mofa_opalin_cp, factors = 1:4, color_by = "condition")

#Factor characterization



#table instead makes it easier to read

# get top weights as a table
opalin_weights_f1 <- get_weights(mofa_opalin_cp, 
                                 views = "opalin", 
                                 factors = 1, 
                                 as.data.frame = TRUE)

opalin_weights_f1 <- opalin_weights_f1[order(abs(opalin_weights_f1$value), decreasing = TRUE), ]
head(opalin_weights_f1, 20)

#match to gene names from my opalin_anno

opalin_anno_df <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/opalin_annotation.csv")

# check peak format in annotation
head(opalin_anno_df[, c("seqnames", "start", "end", "GENENAME")])


# create matching key in annotation
opalin_anno_df$peak_key <- paste(gsub("chr", "", opalin_anno_df$seqnames),
                                 opalin_anno_df$start,
                                 opalin_anno_df$end, sep = "-")

# get top 20 peaks
top_peaks_f1 <- head(opalin_weights_f1, 20)

# match to annotation
top_peaks_annotated <- merge(top_peaks_f1,
                             opalin_anno_df[, c("peak_key", "annotation", "GENENAME", "distanceToTSS")],
                             by.x = "feature",
                             by.y = "peak_key")

top_peaks_annotated[, c("feature", "value", "annotation", "GENENAME", "distanceToTSS")]

# take first peak from weights
top_peaks_f1$feature[1]

# check if it exists in annotation key
head(opalin_anno_df$peak_key)
grep("82730328", opalin_anno_df$peak_key)


# create matching key
opalin_anno_cp$peak_key <- paste(gsub("chr", "", opalin_anno_cp$seqnames),
                                 opalin_anno_cp$start,
                                 opalin_anno_cp$end, sep = "-")

# match top weighted peaks to annotation
top_peaks_annotated <- merge(top_peaks_f1,
                             opalin_anno_cp[, c("peak_key", "annotation", "GENENAME.x", "distanceToTSS")],
                             by.x = "feature",
                             by.y = "peak_key")

top_peaks_annotated[, c("feature", "value", "annotation", "GENENAME.x", "distanceToTSS")]
write.csv(top_peaks_annotated,
          "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/mofa_opalin_cp_factor1_top_peaks.csv",
          row.names = FALSE)

#view table
 clus_opalin_table <- read.csv("/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/mofa_opalin_cp_factor1_top_peaks.csv")

 
 plot_top_weights(mofa_opalin_cp,
                  view    = "opalin",
                  factor  = 1,
                  nfeatures = 20,
                  scale   = TRUE)
 
 
 
 #FACTOR 2 (OPALIN VS BULK)
 #TOO difficult to see in this plot
 plot_weights(mofa_opalin_cp,
              view      = "bulk",
              factor    = 2,
              nfeatures = 20)
 
 bulk_weights_f2 <- get_weights(mofa_opalin_cp,
                                views = "bulk",
                                factors = 2,
                                as.data.frame = TRUE)
 
 bulk_weights_f2 <- bulk_weights_f2[order(abs(bulk_weights_f2$value), decreasing = TRUE), ]
 
 #annotation from GO
 # annotate
 bulk_anno_df$peak_key <- paste(bulk_anno_df$seqnames, bulk_anno_df$start, bulk_anno_df$end, sep = "_")
 
 bulk_weights_f2_annotated <- merge(bulk_weights_f2,
                                    bulk_anno_df[, c("peak_key", "annotation", "GENENAME", "distanceToTSS")],
                                    by.x = "feature",
                                    by.y = "peak_key",
                                    all.x = TRUE)
 
 bulk_weights_f2_annotated <- bulk_weights_f2_annotated[order(abs(bulk_weights_f2_annotated$value), decreasing = TRUE), ]
 
 write.csv(head(bulk_weights_f2_annotated),
           "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/mofa_opalin_cp_factor2_bulk_top20_weights.csv",
           row.names = FALSE)
 

 #FACTOR 3 (BULK VS OPALIN)
 
 bulk_weights_f3 <- get_weights(mofa_opalin_cp,
                                views = "bulk",
                                factors = 3,
                                as.data.frame = TRUE)
 
 bulk_weights_f3 <- bulk_weights_f3[order(abs(bulk_weights_f3$value), decreasing = TRUE), ]
 top_peaks_f3 <- head(bulk_weights_f3, 20)
 
 top_peaks_f3_annotated <- merge(top_peaks_f3,
                                 bulk_anno_df[, c("peak_key", "annotation", "GENENAME", "distanceToTSS")],
                                 by.x = "feature",
                                 by.y = "peak_key")
 
 top_peaks_f3_annotated <- top_peaks_f3_annotated[!duplicated(top_peaks_f3_annotated$feature), ]
 
 top_peaks_f3_annotated[, c("feature", "value", "annotation", "GENENAME", "distanceToTSS")]
 
 write.csv(top_peaks_f3_annotated,
           "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/mofa_opalin_cp_factor3_bulk_top20_weights.csv",
           row.names = FALSE)
 

 #factor 4
 
 bulk_weights_f4 <- get_weights(mofa_opalin_cp,
                                views = "bulk",
                                factors = 4,
                                as.data.frame = TRUE)
 
 bulk_weights_f4 <- bulk_weights_f4[order(abs(bulk_weights_f4$value), decreasing = TRUE), ]
 top_peaks_f4 <- head(bulk_weights_f4, 20)
 
 top_peaks_f4_annotated <- merge(top_peaks_f4,
                                 bulk_anno_df[, c("peak_key", "annotation", "GENENAME", "distanceToTSS")],
                                 by.x = "feature",
                                 by.y = "peak_key")
 
 top_peaks_f4_annotated <- top_peaks_f4_annotated[!duplicated(top_peaks_f4_annotated$feature), ]
 
 top_peaks_f4_annotated[, c("feature", "value", "annotation", "GENENAME", "distanceToTSS")]
 
 write.csv(top_peaks_f4_annotated,
           "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/jd_analysis_sc/saved_objects/oligo_data_out/mofa_opalin_cp_factor4_bulk_top20_weights.csv",
           row.names = FALSE)
 
 #GENE SET ENRICHMENT ANALYSIS (FOLLOWS TUTORIAL)
 


 library(MOFAdata)
 
 # load reactome gene sets
 utils::data(reactomeGS)
 
 # check dimensions
 dim(reactomeGS)
 head(rownames(reactomeGS))
 head(colnames(reactomeGS))
 
 #we have cooridnate while ensembles uses IDs as columns. need to map peaks to ensemble ids using annotaions 
 
 # get ensembl IDs for opalin peaks
 opalin_anno_cp$peak_key2 <- paste(gsub("chr", "", opalin_anno_cp$seqnames),
                                   opalin_anno_cp$start,
                                   opalin_anno_cp$end, sep = "-")
 
 # subset to HVF peaks only
 hvf_peaks <- rownames(opalin_hvf)
 opalin_hvf_anno <- opalin_anno_cp[opalin_anno_cp$peak_key2 %in% hvf_peaks, ]
 
 # check ensembl column
 head(opalin_hvf_anno$ENSEMBL)
 
 # keep only peaks that have ensembl IDs in reactomeGS
 opalin_hvf_anno_filt <- opalin_hvf_anno[opalin_hvf_anno$ENSEMBL %in% colnames(reactomeGS), ]
 
 # create binary peak x pathway matrix
 # for each peak, find which pathways its gene belongs to
 peak_pathway_matrix <- matrix(0, 
                               nrow = nrow(reactomeGS),
                               ncol = nrow(opalin_hvf_anno_filt),
                               dimnames = list(rownames(reactomeGS), 
                                               opalin_hvf_anno_filt$peak_key2))
 
 for (i in seq_len(nrow(opalin_hvf_anno_filt))) {
   ensembl_id <- opalin_hvf_anno_filt$ENSEMBL[i]
   if (ensembl_id %in% colnames(reactomeGS)) {
     peak_pathway_matrix[, i] <- reactomeGS[, ensembl_id]
   }
 }
 
 # check dimensions
 dim(peak_pathway_matrix)
 
 # run GSEA on positive weights
 res_positive <- run_enrichment(mofa_opalin_cp,
                                feature.sets = peak_pathway_matrix,
                                view         = "opalin",
                                sign         = "positive")
 
 # run GSEA on negative weights
 res_negative <- run_enrichment(mofa_opalin_cp,
                                feature.sets = peak_pathway_matrix,
                                view         = "opalin",
                                sign         = "negative")
 
 
 # plot positive weights enrichment
 plot_enrichment(res_positive, factor = 1, max.pathways = 20)
 
 # plot negative weights enrichment
 plot_enrichment(res_negative, factor = 1, max.pathways = 20)
 
 
 #unsure of the above analysis
 
 
 
 
 #Repeat MOFA FOR BULK VS PLEKHG1
 # select top 5000 HVFs
 plekhg1_vars    <- apply(plekhg1_pseudo_norm, 1, var)
 plekhg1_top_idx <- head(order(plekhg1_vars, decreasing = TRUE), 5000)
 plekhg1_hvf     <- plekhg1_pseudo_norm[plekhg1_top_idx, ]
 
 # create 2-view MOFA object
 mofa_input_plekhg1 <- list(
   bulk    = as.matrix(bulk_hvf),
   plekhg1 = as.matrix(plekhg1_hvf)
 )
 
 mofa_plekhg1_cp <- create_mofa(mofa_input_plekhg1)
 plot_data_overview(mofa_plekhg1_cp)
 
 # set options
 data_opts  <- get_default_data_options(mofa_plekhg1_cp)
 model_opts <- get_default_model_options(mofa_plekhg1_cp)
 model_opts$num_factors <- 9
 
 train_opts <- get_default_training_options(mofa_plekhg1_cp)
 train_opts$convergence_mode <- "slow"
 train_opts$seed <- 42
 
 mofa_plekhg1_cp <- prepare_mofa(mofa_plekhg1_cp,
                                 data_options     = data_opts,
                                 model_options    = model_opts,
                                 training_options = train_opts)
 
 mofa_plekhg1_cp <- run_mofa(mofa_plekhg1_cp,
                             outfile      = "mofa_plekhg1_clusterpeaks.hdf5",
                             use_basilisk = TRUE)
 
 samples_metadata(mofa_plekhg1_cp) <- sample_meta
 
 plot_variance_explained(mofa_plekhg1_cp, max_r2 = 15)
 plot_factor_cor(mofa_plekhg1_cp)
 