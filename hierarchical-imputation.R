install.packages("future.apply")
library(mice)
library(lattice)
library(VIM)
library(forcats)
library(dplyr)
library(car)
library(parallel)  # Parallel processing library

# Step 1: Read the processed CSV data
processed_data_updated <- read.csv("C:/Users/MARAM/Desktop/processed_data.csv")

# Step 2: Feature selection 
variables_to_keep <- c(
  "id_menage", "Id_Ind", "Poids_Final",  # Basic identifiers and weights
  "V100", "V101", "V104",                # Demographic variables (Gender, Age, Nationality)
  "V105_A", "V106_M", "V109", "V115_A",  # Additional demographic variables
  "V301", "V302", "V303", "V305", "V403", "V405",  # Employment-related variables
  "V702", "V710", "V713_A","V808_C", "V814", "V815_B", # Other important variables
  "V606", "V608_A",                      # Migration and employment-related variables
  "V107",                                # Country of migration (for clustering)
  "V315_3",                              # Reason for return
  "V308_M", "V308_F", "V309_M", "V309_F", # Numeric variables
  "V603"                                 # 
)

# Step 3: Extract and keep only the selected variables
processed_data_cleaned <- processed_data_updated[, variables_to_keep]

# Define variables based on types
numeric_vars <- c("V106_M", "V109", "V308_M", "V308_F", "V309_M", "V309_F", "V315_3")
categorical_vars <- c("V603", "V815_B", "V606", "V608_A", "V808_C", "V814", "V405", 
                      "V710", "V713_A", "V115_A", "V303", "V403")


# Combine all specified variables into one list
specified_vars <- c(numeric_vars, categorical_vars)

# Check for extra variables in processed_data_cleaned
extra_vars <- setdiff(names(processed_data_cleaned), specified_vars)

# If there are extra variables, create a new dataset with only the specified variables
if (length(extra_vars) > 0) {
  processed_data_cleaned_filtered <- processed_data_cleaned[, specified_vars]
  print("Extra variables found and removed:")
  print(extra_vars)
} else {
  processed_data_cleaned_filtered <- processed_data_cleaned
  print("No extra variables found.")
}

# Display summary of the new dataset to verify it only contains specified variables
summary(processed_data_cleaned_filtered)

summary(imputed_data)




# Step 4: Remove duplicate rows
processed_data_cleaned <- processed_data_cleaned[!duplicated(processed_data_cleaned), ]

# Step 5: Handle non-response categories and filter them out
processed_data_cleaned <- processed_data_cleaned %>%
  filter(V702 != "Refuse de répondre")

# Step 6: Convert character columns to factors
processed_data_cleaned <- processed_data_cleaned %>%
  mutate(across(where(is.character), as.factor))

# Step 7: Randomly sample a subset of households 
set.seed(123)  # Ensure reproducibility
sampled_households <- sample(unique(processed_data_cleaned$id_menage), 
                             size = 0.5 * length(unique(processed_data_cleaned$id_menage)))

# Step 8: Filter the dataset to keep only individuals from sampled households
processed_data_cleaned <- processed_data_cleaned %>%
  filter(id_menage %in% sampled_households)

# Step 9: Ensure household structure by keeping households with more than 1 individual
processed_data_cleaned <- processed_data_cleaned %>%
  group_by(id_menage) %>%
  filter(n() > 1) %>%
  ungroup()


# Step 11: Collapse rare levels for high-cardinality variables
collapse_rare_levels <- function(factor_var, threshold = 20) {  
  levels_to_group <- names(which(table(factor_var) < threshold))
  if (length(levels_to_group) > 0) {
    factor_var <- fct_collapse(factor_var, Other = levels_to_group)
  }
  return(factor_var)
}

# Apply collapsing to the high-cardinality variables
processed_data_cleaned <- processed_data_cleaned %>%
  mutate(
    V106_M = collapse_rare_levels(V106_M, threshold = 20),
    V305 = collapse_rare_levels(V305, threshold = 20),
    V702 = collapse_rare_levels(V702, threshold = 20),
    V606 = collapse_rare_levels(V606, threshold = 20),
    V608_A = collapse_rare_levels(V608_A, threshold = 20)
  )

# Step 12: Convert '9999' in factor variables to NA (Convert factor to character first)
processed_data_cleaned$V109 <- as.factor(na_if(as.character(processed_data_cleaned$V109), "9999"))
processed_data_cleaned$V117 <- as.factor(na_if(as.character(processed_data_cleaned$V117), "9999"))
processed_data_cleaned$V115_A <- as.factor(na_if(as.character(processed_data_cleaned$V115_A), "9999"))

# Step 13: Prepare the dataset for MICE
processed_data_cleaned$V107_cluster <- as.integer(processed_data_cleaned$V107)  # Convert V107 to integer for clustering

# Step 14: Set up imputation methods based on variable type
imputation_methods <- make.method(processed_data_cleaned)

# Define specific methods for variables to be imputed
for (var in names(processed_data_cleaned)) {
  if (is.factor(processed_data_cleaned[[var]])) {
    if (nlevels(processed_data_cleaned[[var]]) > 2) {
      if (nlevels(processed_data_cleaned[[var]]) > 10) {
        imputation_methods[var] <- "cart"  # Use CART for high-cardinality categorical variables
      } else {
        imputation_methods[var] <- "polyreg"  # Use polyreg for categorical variables with >2 levels
      }
    } else {
      imputation_methods[var] <- "logreg"  # Use logreg for binary categorical variables
    }
  } else if (var %in% numeric_vars) {
    imputation_methods[var] <- "pmm"  # Use pmm for numeric variables
  } else {
    # Skip imputation for variables not in numeric_vars or categorical_vars
    imputation_methods[var] <- ""  # Set to empty to skip imputation
  }
}

# Verify the assigned methods
print(imputation_methods)

# Step 15: Reconfigure the predictor matrix
# Ensure 'Poids_Final', 'id_menage', 'Id_Ind', and other non-imputed variables are excluded from imputation
imputation_methods[c("Poids_Final", "id_menage", "Id_Ind")] <- ""  # Exclude these variables
predictor_matrix <- make.predictorMatrix(processed_data_cleaned)
predictor_matrix[, "V107_cluster"] <- -2  # Use V107_cluster only as a predictor
predictor_matrix[, c("id_menage", "Id_Ind", "Poids_Final")] <- 0
predictor_matrix["id_menage", ] <- 0
predictor_matrix["Id_Ind", ] <- 0
predictor_matrix["Poids_Final", ] <- 0

# Step 16: Run MICE with specified methods and predictor matrix
library(future)
library(future.apply)
plan(multisession)  


# Run MICE with adjusted methods and predictor matrix
imputed_data <- mice(
  processed_data_cleaned, 
  method = imputation_methods, 
  predictorMatrix = predictor_matrix, 
  m = 10,         # Number of imputations 
  maxit = 250,     # Number of iterations
  seed = 123,    # Reproducibility
  printFlag = TRUE,
  future.seed = TRUE
)
plot(imputed_data,trace="all")


# Define the 19 variables
continuous_vars <- c("V106_M", "V109", "V308_M", "V308_F", "V309_M", "V309_F", "V315_3")
categorical_vars <- c("V603", "V815_B", "V606", "V608_A", "V808_C", "V814", 
                      "V405", "V710", "V713_A", "V115_A", "V303", "V403")

# Combine all specified variables into one list
selected_vars <- c(continuous_vars, categorical_vars)

# Extract only the 19 desired variables from processed_data_cleaned
observed_data <- processed_data_cleaned %>%
  select(all_of(selected_vars))

# Display the summary of the new observed_data
print("Summary of Observed Data:")
summary(observed_data)

summary(imputed_data)

# Calculate RMSE for each continuous variable
continuous_vars <- c("V106_M", "V109", "V308_M", "V308_F", "V309_M", "V309_F", "V315_3")

rmse_results <- sapply(continuous_vars, function(var) {
  observed <- observed_data[[var]]
  imputed <- complete(imputed_data, action = 1)[[var]]
  
  # Filter non-missing entries for comparison
  observed_non_missing <- observed[!is.na(observed)]
  imputed_non_missing <- imputed[!is.na(observed)]
  
  # Compute RMSE
  sqrt(mean((observed_non_missing - imputed_non_missing)^2, na.rm = TRUE))
})

print("RMSE for each continuous variable:")
print(rmse_results)

# Plot trace plot for a specific variable (V803)
plot(imputed_data, vars = "V803", trace = TRUE)

install.packages("e1071")
library(e1071)
library(dplyr)

# Convert continuous variables to numeric
processed_data_cleaned[continuous_vars] <- lapply(processed_data_cleaned[continuous_vars], function(x) as.numeric(as.character(x)))

# Check if all continuous variables are now numeric
str(processed_data_cleaned[continuous_vars])  # This will show if the columns are properly converted to numeric


# Calculate summary statistics using only observed data (non-missing values)
summary_stats <- data.frame(
  Variable = continuous_vars,
  Description = c(
    "Number of years since migration", 
    "Number of children in the household", 
    "Male household members", 
    "Female household members", 
    "Male children in the household", 
    "Female children in the household", 
    "Reason for return migration"
  ),
  Mean = sapply(continuous_vars, function(var) mean(processed_data_cleaned[[var]], na.rm = TRUE)),
  Median = sapply(continuous_vars, function(var) median(processed_data_cleaned[[var]], na.rm = TRUE)),
  SD = sapply(continuous_vars, function(var) sd(processed_data_cleaned[[var]], na.rm = TRUE)),
  Min = sapply(continuous_vars, function(var) min(processed_data_cleaned[[var]], na.rm = TRUE)),
  Max = sapply(continuous_vars, function(var) max(processed_data_cleaned[[var]], na.rm = TRUE)),
  Missing_Rate = sapply(continuous_vars, function(var) mean(is.na(processed_data_cleaned[[var]])) * 100),
  Skewness = sapply(continuous_vars, function(var) skewness(processed_data_cleaned[[var]], na.rm = TRUE)),
  Kurtosis = sapply(continuous_vars, function(var) kurtosis(processed_data_cleaned[[var]], na.rm = TRUE))
)

# Display the table
print(summary_stats)
summary(imputed_data)
summary(observed_data)






# observed vs imputed data
# Add pseudo-counts to avoid zero frequencies
all_levels <- union(levels(observed_data[[var]]), levels(mice_imputed_data[[var]]))

# Align frequencies with all levels
obs_freq <- table(factor(observed_data[[var]], levels = all_levels)) + 1
mice_freq <- table(factor(mice_imputed_data[[var]], levels = all_levels)) + 1
# Install DescTools if not already installed
if (!requireNamespace("DescTools", quietly = TRUE)) {
  install.packages("DescTools")
}

library(DescTools)

# Perform G-Test
g_test_result <- GTest(obs_freq, p = prop.table(mice_freq))

# Print G-Test results
print(g_test_result)

# Set a threshold for minimum frequencies
threshold <- 5

# Filter frequencies by threshold
sufficient_levels <- names(obs_freq[obs_freq >= threshold])

# Update frequencies to include only sufficient levels
obs_freq <- obs_freq[sufficient_levels]
mice_freq <- mice_freq[sufficient_levels]

# Perform the Chi-Square Test again
test_result <- chisq.test(obs_freq, p = prop.table(mice_freq), simulate.p.value = TRUE)
print(test_result)

# Loop through all categorical variables
mice_chi_square_results <- lapply(categorical_vars, function(var) {
  # Align levels
  all_levels <- union(levels(observed_data[[var]]), levels(mice_imputed_data[[var]]))
  obs_freq <- table(factor(observed_data[[var]], levels = all_levels)) + 1
  mice_freq <- table(factor(mice_imputed_data[[var]], levels = all_levels)) + 1
  
  # Apply threshold
  sufficient_levels <- names(obs_freq[obs_freq >= threshold])
  obs_freq <- obs_freq[sufficient_levels]
  mice_freq <- mice_freq[sufficient_levels]
  
  # Perform Chi-Square Test if valid
  if (length(obs_freq) > 1) {
    test_result <- chisq.test(obs_freq, p = prop.table(mice_freq), simulate.p.value = TRUE)
    list(
      Variable = var,
      Chi_Square_Statistic = test_result$statistic,
      P_Value = test_result$p.value
    )
  } else {
    # Return NA for insufficient data
    list(
      Variable = var,
      Chi_Square_Statistic = NA,
      P_Value = NA
    )
  }
})

# Compile results into a data frame
mice_results_summary <- do.call(rbind, lapply(mice_chi_square_results, function(x) {
  data.frame(
    Variable = x$Variable,
    Chi_Square_Statistic = x$Chi_Square_Statistic,
    P_Value = x$P_Value
  )
}))

# Print the results
print(mice_results_summary)


install.packages("moments")
library(moments)
# Check the structure of the imputed data
str(mice_imputed_data[, continuous_vars])

# Convert any non-numeric continuous variables to numeric
mice_imputed_data[continuous_vars] <- lapply(
  mice_imputed_data[continuous_vars],
  function(x) as.numeric(as.character(x))
)

# Recheck the structure
str(mice_imputed_data[, continuous_vars])



compute_mice_summary_stats <- function(data) {
  stats <- data.frame(
    Variable = names(data),
    Min = sapply(data, min, na.rm = TRUE),
    `1st Qu.` = sapply(data, function(x) quantile(x, 0.25, na.rm = TRUE)),
    Median = sapply(data, median, na.rm = TRUE),
    Mean = sapply(data, mean, na.rm = TRUE),
    `3rd Qu.` = sapply(data, function(x) quantile(x, 0.75, na.rm = TRUE)),
    Max = sapply(data, max, na.rm = TRUE),
    Skewness = sapply(data, skewness, na.rm = TRUE),
    Kurtosis = sapply(data, kurtosis, na.rm = TRUE)
  )
  return(stats)
}
# Compute enhanced summary statistics for continuous variables
mice_summary_stats <- compute_mice_summary_stats(mice_imputed_data[, continuous_vars])

# View the summary
print(mice_summary_stats)




listwise_data <- processed_data_cleaned[complete.cases(processed_data_cleaned[, continuous_vars]), ]
summary(listwise_data)


# Listwise deletion dataset (complete cases only)
listwise_data_cat <- processed_data_cleaned[complete.cases(processed_data_cleaned[, categorical_vars]), ]

# Generate frequencies for each categorical variable
listwise_frequencies <- lapply(categorical_vars, function(var) {
  table(listwise_data_cat[[var]])
})

print(listwise_frequencies)






# Check structure of imputed_data and processed_data_cleaned
str(imputed_data)  # Should show a data frame with columns and rows
str(processed_data_cleaned)  # Should show the original data structure

# View column names
print(colnames(imputed_data))
print(colnames(processed_data_cleaned))

# Extract the first imputed dataset
imputed_data <- as.data.frame(complete(imputed_data, action = 1))

# Verify dimensions
print(dim(imputed_data))
print(head(imputed_data))  # Check if it contains the expected data




# Save each imputed dataset to separate CSVs
for (i in 1:imputed_data$m) {
  imputed_dataset <- complete(imputed_data, action = i)
  write.csv(imputed_dataset, paste0("C:/Users/MARAM/Desktop/MICE_Imputed_Data_", i, ".csv"), row.names = FALSE)
}

# Save the first imputed dataset separately (optional)
imputed_dataset_first <- complete(imputed_data, action = 1)
write.csv(imputed_dataset_first, "C:/Users/MARAM/Desktop/MICE_Imputed_Data_First.csv", row.names = FALSE)

# Save the entire MICE object for future use
saveRDS(imputed_data, file = "C:/Users/MARAM/Desktop/MICE_Imputed_Object.rds")

# Test loading the saved MICE object to ensure it works
loaded_imputed_data <- readRDS("C:/Users/MARAM/Desktop/MICE_Imputed_Object.rds")
print("MICE object successfully loaded!")
str(observed_data)




library(mice)
library(coda)

# Function to calculate PSRF for one-hot encoded categorical variables
calculate_psrf_categorical <- function(var, imputed_data) {
  # Create dummy variables for each category across imputations
  dummy_data <- lapply(1:imputed_data$m, function(i) {
    imp_data <- complete(imputed_data, action = i)
    model.matrix(~ . - 1, data = imp_data[, var, drop = FALSE])  # One-hot encoding
  })
  
  # Align levels across imputations
  all_columns <- unique(unlist(lapply(dummy_data, colnames)))
  aligned_dummy_data <- lapply(dummy_data, function(df) {
    df <- as.data.frame(df)
    missing_cols <- setdiff(all_columns, colnames(df))
    for (col in missing_cols) df[[col]] <- 0  # Add missing columns with zeros
    df <- df[, all_columns, drop = FALSE]    # Ensure column order
    as.matrix(df)                            # Convert back to matrix
  })
  
  # Debugging step: Check the structure of aligned dummy data
  print(paste("Aligned dummy data for:", var))
  print(dim(aligned_dummy_data[[1]]))
  
  # Create MCMC chains for each binary variable
  mcmc_chains <- lapply(1:ncol(aligned_dummy_data[[1]]), function(j) {
    tryCatch({
      as.mcmc.list(lapply(aligned_dummy_data, function(df) as.mcmc(df[, j, drop = TRUE])))
    }, error = function(e) {
      NULL  # Handle any failures gracefully
    })
  })
  
  # Filter out NULL chains
  mcmc_chains <- Filter(Negate(is.null), mcmc_chains)
  
  # Calculate PSRF for each binary variable
  psrf_results <- lapply(mcmc_chains, function(chain) {
    tryCatch({
      gelman.diag(chain, autoburnin = TRUE)$psrf[1, ]
    }, error = function(e) {
      c(NA, NA)  # Return NA if PSRF calculation fails
    })
  })
  
  # Combine results into a data frame
  if (length(psrf_results) > 0) {
    results <- do.call(rbind, psrf_results)
    colnames(results) <- c("Point_Estimate", "Upper_CI")
    results <- cbind(Category = colnames(aligned_dummy_data[[1]]), results)
    results$Variable <- var
    return(results)
  } else {
    return(data.frame(Variable = var, Category = NA, Point_Estimate = NA, Upper_CI = NA))
  }
}

# List of categorical variables
categorical_vars <- c("V603", "V815_B", "V606", "V608_A", "V808_C", 
                      "V814", "V405", "V710", "V713_A", "V115_A", 
                      "V303", "V403")

# Apply the function to all categorical variables
psrf_categorical_results <- lapply(categorical_vars, function(var) calculate_psrf_categorical(var, imputed_data))

# Combine all results into a single data frame
psrf_categorical_summary <- do.call(rbind, psrf_categorical_results)

# Display the PSRF results for categorical variables
print("PSRF Results for Categorical Variables:")
print(psrf_categorical_summary)




