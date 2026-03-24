#NEW MOFA SCRIPT (CLEANER VERSION)

# bulk
bulk_vars    <- apply(bulk_counts, 1, var)
bulk_top_idx <- head(order(bulk_vars, decreasing = TRUE), 5000)
bulk_hvf     <- bulk_counts[bulk_top_idx, ]

# opalin
opalin_vars    <- apply(opalin_vst, 1, var)
opalin_top_idx <- head(order(opalin_vars, decreasing = TRUE), 5000)
opalin_hvf     <- opalin_vst[opalin_top_idx, ]

# recreate MOFA object
mofa_input_opalin <- list(
  bulk   = as.matrix(bulk_hvf),
  opalin = as.matrix(opalin_hvf)
)

mofa_opalin <- create_mofa(mofa_input_opalin)

data_opts  <- get_default_data_options(mofa_opalin)
model_opts <- get_default_model_options(mofa_opalin)
model_opts$num_factors <- 3

train_opts <- get_default_training_options(mofa_opalin)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

mofa_opalin <- prepare_mofa(mofa_opalin,
                            data_options     = data_opts,
                            model_options    = model_opts,
                            training_options = train_opts)

mofa_opalin <- run_mofa(mofa_opalin,
                        outfile      = "mofa_opalin_model.hdf5",
                        use_basilisk = TRUE)


# sanity check
stopifnot(all(sort(sample_meta$sample) == sort(unlist(samples_names(mofa_opalin)))))

# add metadata
#samples_metadata(mofa_opalin) <- sample_meta

sample_meta <- data.frame(
  sample    = c("PD1222", "PD1231", "PD726", "PD833", "PD936",
                "PDC167", "PDC184", "PDC197", "PDC87"),
  condition = c("PD", "PD", "PD", "PD", "PD",
                "Control", "Control", "Control", "Control")
)



plot_factor_cor(mofa_opalin)
plot_variance_explained(mofa_opalin, max_r2 = 15)
#shows variance explained per view
head(get_variance_explained(mofa_opalin)$r2_total[[1]])
#plot variance per view 
plot_variance_explained(mofa_opalin, plot_total = T)[[2]]

# Variance explained for every factor in per view
head(get_variance_explained(mofa_opalin)$r2_per_factor[[1]])

#Exlporing factor 1 ( 64% bulk)


samples_metadata(mofa_opalin) <- sample_meta

library(GGally)
#plotting factors
plot_factor(mofa_opalin, factors = 1:3, color_by = "condition")

plot_weights(mofa_opalin, factor = 3, view = "bulk", nfeatures = 10)


#creating a 4 view mofa to compare all olgio subtypes with bulk
# select HVFs for all subtypes
plekhg1_vars    <- apply(plekhg1_vst, 1, var)
plekhg1_top_idx <- head(order(plekhg1_vars, decreasing = TRUE), 5000)
plekhg1_hvf     <- plekhg1_vst[plekhg1_top_idx, ]

opc_vars        <- apply(opc_vst, 1, var)
opc_top_idx     <- head(order(opc_vars, decreasing = TRUE), 5000)
opc_hvf         <- opc_vst[opc_top_idx, ]

# create 4-view MOFA object
mofa_input_all <- list(
  bulk    = as.matrix(bulk_hvf),
  opalin  = as.matrix(opalin_hvf),
  plekhg1 = as.matrix(plekhg1_hvf),
  opc     = as.matrix(opc_hvf)
)

mofa_all <- create_mofa(mofa_input_all)
plot_data_overview(mofa_all)

# set options
data_opts_all  <- get_default_data_options(mofa_all)
model_opts_all <- get_default_model_options(mofa_all)
model_opts_all$num_factors <- 3

train_opts_all <- get_default_training_options(mofa_all)
train_opts_all$convergence_mode <- "slow"
train_opts_all$seed <- 42

mofa_all <- prepare_mofa(mofa_all,
                         data_options     = data_opts_all,
                         model_options    = model_opts_all,
                         training_options = train_opts_all)

mofa_all <- run_mofa(mofa_all,
                     outfile      = "mofa_all_model.hdf5",
                     use_basilisk = TRUE)
# sanity check
stopifnot(all(sort(sample_meta$sample) == sort(unlist(samples_names(mofa_opalin)))))

# add metadata
#samples_metadata(mofa_opalin) <- sample_meta

sample_meta <- data.frame(
  sample    = c("PD1222", "PD1231", "PD726", "PD833", "PD936",
                "PDC167", "PDC184", "PDC197", "PDC87"),
  condition = c("PD", "PD", "PD", "PD", "PD",
                "Control", "Control", "Control", "Control")
)

plot_factor_cor(mofa_all)
plot_variance_explained(mofa_all, max_r2 = 15)
plot_factor(mofa_all, factors = 1:3, color_by = "condition")

#shows variance explained per view
head(get_variance_explained(mofa_all)$r2_total[[1]])
#plot variance per view 
plot_variance_explained(mofa_all, plot_total = T)[[2]]

# Variance explained for every factor in per view
head(get_variance_explained(mofa_all)$r2_per_factor[[1]])
