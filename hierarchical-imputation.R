# 01_data_cleaning.R
# Loads and cleans raw household survey data

library(dplyr)
library(forcats)

# Load the dataset
processed_data_updated <- read.csv("data/processed_data.csv")

# Remove duplicate rows
processed_data_cleaned <- processed_data_updated %>%
  distinct()

# Remove non-response entries
processed_data_cleaned <- processed_data_cleaned %>%
  filter(V702 != "Refuse de répondre")

# Convert character columns to factors
processed_data_cleaned <- processed_data_cleaned %>%
  mutate(across(where(is.character), as.factor))

# Randomly sample 50% of households
set.seed(123)
sampled_households <- sample(unique(processed_data_cleaned$id_menage), 
                             size = 0.5 * length(unique(processed_data_cleaned$id_menage)))

# Keep only individuals from sampled households with >1 individual
processed_data_cleaned <- processed_data_cleaned %>%
  filter(id_menage %in% sampled_households) %>%
  group_by(id_menage) %>%
  filter(n() > 1) %>%
  ungroup()

# 02_feature_selection.R
# Selects relevant variables and collapses rare factor levels

library(dplyr)
library(forcats)

numeric_vars <- c("V106_M", "V109", "V308_M", "V308_F", "V309_M", "V309_F", "V315_3")
categorical_vars <- c("V603", "V815_B", "V606", "V608_A", "V808_C", "V814", 
                      "V405", "V710", "V713_A", "V115_A", "V303", "V403")
specified_vars <- c(numeric_vars, categorical_vars)

# Function to collapse rare levels
collapse_rare_levels <- function(factor_var, threshold = 20) {
  levels_to_group <- names(which(table(factor_var) < threshold))
  if (length(levels_to_group) > 0) {
    factor_var <- fct_collapse(factor_var, Other = levels_to_group)
  }
  return(factor_var)
}

# Apply collapsing
processed_data_cleaned <- processed_data_cleaned %>%
  mutate(
    V106_M = collapse_rare_levels(V106_M),
    V305 = collapse_rare_levels(V305),
    V702 = collapse_rare_levels(V702),
    V606 = collapse_rare_levels(V606),
    V608_A = collapse_rare_levels(V608_A)
  )

# 03_mice_imputation.R
# Performs multilevel imputation using MICE

library(mice)
library(future)
library(future.apply)

# Prepare clustering variable
processed_data_cleaned$V107_cluster <- as.integer(processed_data_cleaned$V107)

# Setup MICE methods
imputation_methods <- make.method(processed_data_cleaned)
for (var in names(processed_data_cleaned)) {
  if (is.factor(processed_data_cleaned[[var]])) {
    imputation_methods[var] <- ifelse(nlevels(processed_data_cleaned[[var]]) > 10, "cart",
                                      ifelse(nlevels(processed_data_cleaned[[var]]) > 2, "polyreg", "logreg"))
  } else if (var %in% numeric_vars) {
    imputation_methods[var] <- "pmm"
  } else {
    imputation_methods[var] <- ""
  }
}
imputation_methods[c("Poids_Final", "id_menage", "Id_Ind")] <- ""

# Set up predictor matrix
predictor_matrix <- make.predictorMatrix(processed_data_cleaned)
predictor_matrix[, "V107_cluster"] <- -2
predictor_matrix[, c("id_menage", "Id_Ind", "Poids_Final")] <- 0
predictor_matrix["id_menage", ] <- 0
predictor_matrix["Id_Ind", ] <- 0

# Run MICE
plan(multisession)
imputed_data <- mice(processed_data_cleaned,
                     method = imputation_methods,
                     predictorMatrix = predictor_matrix,
                     m = 10,
                     maxit = 250,
                     seed = 123,
                     printFlag = TRUE,
                     future.seed = TRUE)

saveRDS(imputed_data, "outputs/MICE_Imputed_Object.rds")


# 03_mice_imputation.R
# Performs multilevel imputation using MICE

library(mice)
library(future)
library(future.apply)

# Prepare clustering variable
processed_data_cleaned$V107_cluster <- as.integer(processed_data_cleaned$V107)

# Setup MICE methods
imputation_methods <- make.method(processed_data_cleaned)
for (var in names(processed_data_cleaned)) {
  if (is.factor(processed_data_cleaned[[var]])) {
    imputation_methods[var] <- ifelse(nlevels(processed_data_cleaned[[var]]) > 10, "cart",
                                      ifelse(nlevels(processed_data_cleaned[[var]]) > 2, "polyreg", "logreg"))
  } else if (var %in% numeric_vars) {
    imputation_methods[var] <- "pmm"
  } else {
    imputation_methods[var] <- ""
  }
}
imputation_methods[c("Poids_Final", "id_menage", "Id_Ind")] <- ""

# Set up predictor matrix
predictor_matrix <- make.predictorMatrix(processed_data_cleaned)
predictor_matrix[, "V107_cluster"] <- -2
predictor_matrix[, c("id_menage", "Id_Ind", "Poids_Final")] <- 0
predictor_matrix["id_menage", ] <- 0
predictor_matrix["Id_Ind", ] <- 0

# Run MICE
plan(multisession)
imputed_data <- mice(processed_data_cleaned,
                     method = imputation_methods,
                     predictorMatrix = predictor_matrix,
                     m = 10,
                     maxit = 250,
                     seed = 123,
                     printFlag = TRUE,
                     future.seed = TRUE)

saveRDS(imputed_data, "outputs/MICE_Imputed_Object.rds")


# 04_jomo_imputation.R
# Performs hierarchical imputation using JOMO with MCMC diagnostics

library(jomo)
library(coda)
library(mcmcse)
library(haven)
library(dplyr)

# Load and prepare data
migrants_actuels <- read_sav("data/Tunisia-HIMS_migrants_actuels_2021.sav")
categorical_vars <- c("V101", "V108", "V107", "V109", "V114")
migrants_actuels_clean <- migrants_actuels %>%
  mutate(across(all_of(categorical_vars), as.factor))

Y.cat <- migrants_actuels_clean[, categorical_vars]
Y.numcat <- sapply(Y.cat, nlevels)
clus <- migrants_actuels_clean$V107

# MCMC diagnostics
imp_diagnostics <- jomo.MCMCchain(Y = Y.cat, clus = clus, nburn = 1000)
plot(imp_diagnostics$collectbeta[1, 1, 1:1000], type = "l")
acf(imp_diagnostics$collectbeta[1, 1, ])

# Thin chain
thinned_beta_5 <- imp_diagnostics$collectbeta[1, 1, seq(1, 1000, by = 5)]
write.csv(thinned_beta_5, "outputs/thinned_beta_5.csv", row.names = FALSE)

# Final imputation with tuned parameters
imputed_data_jomo <- jomo1rancat(Y.cat = Y.cat, Y.numcat = Y.numcat,
                                 clus = clus, nimp = 10, nburn = 1000, nbetween = 5000)

write.csv(imputed_data_jomo, "outputs/jomo_imputed_data.csv", row.names = FALSE)


# 05_evaluation_metrics.R
# Computes RMSE and Chi-square for imputed vs. observed data

library(e1071)
library(DescTools)
library(moments)

# RMSE for continuous vars
observed_data <- processed_data_cleaned[, numeric_vars]
imputed_sample <- complete(imputed_data, 1)
rmse_results <- sapply(numeric_vars, function(var) {
  obs <- observed_data[[var]]
  imp <- imputed_sample[[var]]
  sqrt(mean((obs[!is.na(obs)] - imp[!is.na(obs)])^2, na.rm = TRUE))
})
print(rmse_results)

# Chi-square for categorical vars
mice_imputed_data <- complete(imputed_data, 1)
chi_results <- lapply(categorical_vars, function(var) {
  obs_freq <- table(factor(observed_data[[var]], levels = levels(mice_imputed_data[[var]]))) + 1
  mice_freq <- table(factor(mice_imputed_data[[var]], levels = levels(mice_imputed_data[[var]]))) + 1
  chisq.test(obs_freq, p = prop.table(mice_freq), simulate.p.value = TRUE)
})


# 06_visualizations_diagnostics.R
# Creates trace plots and convergence diagnostics

library(mice)
library(coda)

# Trace plots for MICE
plot(imputed_data, trace = TRUE)

# PSRF for categorical variables
calculate_psrf_categorical <- function(var, imputed_data) {
  dummy_data <- lapply(1:imputed_data$m, function(i) model.matrix(~ . - 1, complete(imputed_data, i)[, var, drop = FALSE]))
  all_columns <- unique(unlist(lapply(dummy_data, colnames)))
  aligned <- lapply(dummy_data, function(df) {
    df <- as.data.frame(df)
    for (col in setdiff(all_columns, names(df))) df[[col]] <- 0
    df[, all_columns]
  })

  mcmc_chains <- lapply(seq_len(ncol(aligned[[1]])), function(j) {
    tryCatch(as.mcmc.list(lapply(aligned, function(mat) as.mcmc(mat[, j]))), error = function(e) NULL)
  })

  psrf_results <- lapply(mcmc_chains, function(chain) {
    tryCatch(gelman.diag(chain)$psrf[1, ], error = function(e) c(NA, NA))
  })

  data.frame(
    Category = names(dummy_data[[1]]),
    Point_Estimate = sapply(psrf_results, function(x) x[1]),
    Upper_CI = sapply(psrf_results, function(x) x[2]),
    Variable = var
  )
}

# Example
psrf_results_summary <- lapply(categorical_vars, function(v) calculate_psrf_categorical(v, imputed_data))
psrf_results_df <- do.call(rbind, psrf_results_summary)
print(psrf_results_df)


                               
