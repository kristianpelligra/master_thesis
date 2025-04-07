# Set working directory
setwd("C:/Users/User/Desktop/RMT")

# Loading libraries
library(e1071)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# Load and format data
data <- read_delim("heatmap_data_norm.csv", delim = ";") %>% 
  rename(LLV = "Left_Lateral", RLV = "Right_Lateral", TV = "Third", FV = "Fourth", ToV = "Total_Ventricle") %>% 
  filter(Subtype %in% c("AD", "SCI"))

# Get all protein names (columns 7:55)
protein_names <- colnames(data)[7:55]

# Protein groups correlating with key biomarkers
amyloid_associated = c("NPTX2", "NPTXR", "NPTX1", "CHL1", "NCAN", "PTPRN2", "TMEM132D",
                       "CCK", "PAM", "CLSTN1", "OMG", "APLP1", "NRCAM", "VGF", "CDH8",
                       "KIAA1549L")
tau_associated = c("ARPP21", "SNCB", "GAP43", "NRGN", "SPP1", "AQP4", "RPH3A", "SNCA",
                   "AMPH", "DDAH1", "PEBP1")

# Function to analyze combination of two proteins with ventricle volume
analyze_combination <- function(data, protein1, protein2, ventricle_col) {
  model_data <- data %>%
    select(Subtype, !!sym(protein1), !!sym(protein2), !!sym(ventricle_col), Age, Sex)
  
  model_data$Subtype <- factor(model_data$Subtype)
  
  # 5-fold cross validation
  set.seed(123)
  k_folds <- 5 
  preds <- numeric(nrow(model_data))
  actual <- numeric(nrow(model_data))
  folds <- createFolds(model_data$Subtype, k = k_folds)
  
  for (i in 1:k_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- folds[[i]]
    
    train_data <- model_data[train_indices, ]
    test_data <- model_data[test_indices, ]
    
    svm_model <- svm(
      Subtype ~ .,  
      data = train_data,
      kernel = "radial",
      probability = TRUE,
      cost = 1
    )
    
    pred_probs <- predict(svm_model, test_data, probability = TRUE)
    pred_attrs <- attr(pred_probs, "probabilities")
    
    preds[test_indices] <- pred_attrs[, "AD"]
    actual[test_indices] <- ifelse(test_data$Subtype == "AD", 1, 0)
  }
  
  roc_result <- roc(actual, preds)
  auc_value <- auc(roc_result)
  
  # Compute p-value to test if AUC is significantly different from 0.5
  # 3.92 is used because a 95% confidence interval corresponds to approximately ±1.96 standard deviations (z-score)
  # The total range of the 95% CI is therefore 1.96 * 2 = 3.92, allowing us to estimate the standard error
  
  ci_result <- ci.auc(roc_result)
  p_value <- 2 * (1 - pnorm(abs((auc_value - 0.5) / ((ci_result[3] - ci_result[1]) / 3.92))))
  
  return(list(
    protein1 = protein1,
    protein2 = protein2,
    ventricle = ventricle_col,
    auc = auc_value,
    p_value = p_value
  ))
}

# Ventricle columns
ventricle_cols <- c("LLV", "RLV", "TV", "FV", "ToV")

# Perform analysis for all combinations of proteins within and between groups
results <- list()
counter <- 1

# Amyloid-Amyloid combinations
for (i in seq_along(amyloid_associated)) {
  for (j in i:length(amyloid_associated)) {
    for (ventricle in ventricle_cols) {
      results[[counter]] <- analyze_combination(data, amyloid_associated[i], amyloid_associated[j], ventricle)
      counter <- counter + 1
    }
  }
}

# Tau-Tau combinations
for (i in seq_along(tau_associated)) {
  for (j in i:length(tau_associated)) {
    for (ventricle in ventricle_cols) {
      results[[counter]] <- analyze_combination(data, tau_associated[i], tau_associated[j], ventricle)
      counter <- counter + 1
    }
  }
}

# Amyloid-Tau combinations
for (amyloid_protein in amyloid_associated) {
  for (tau_protein in tau_associated) {
    for (ventricle in ventricle_cols) {
      results[[counter]] <- analyze_combination(data, amyloid_protein, tau_protein, ventricle)
      counter <- counter + 1
    }
  }
}

# Convert results to dataframe
results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    Protein1 = x$protein1,
    Protein2 = x$protein2,
    Ventricle = x$ventricle,
    AUC = x$auc,
    P_value = round(x$p_value, 10)
  )
}))

# Ensure p-values are numeric and handle NAs properly
results_df$P_value <- as.numeric(results_df$P_value)
results_df <- results_df %>% filter(!is.na(P_value))

# ---- Code for the first matrix (protein1 + protein2 + ventricle) ----

# Create matrices to store results
all_proteins <- unique(c(amyloid_associated, tau_associated))
result_matrix <- matrix(NA, nrow = length(all_proteins), ncol = length(all_proteins),
                        dimnames = list(all_proteins, all_proteins))
ventricle_matrix <- matrix(NA, nrow = length(all_proteins), ncol = length(all_proteins),
                           dimnames = list(all_proteins, all_proteins))

# Fill the matrix with highest AUC values and corresponding ventricles
for (i in 1:nrow(results_df)) {
  protein1 <- results_df$Protein1[i]
  protein2 <- results_df$Protein2[i]
  auc <- results_df$AUC[i]
  ventricle <- results_df$Ventricle[i]
  
  if (is.na(result_matrix[protein1, protein2]) || auc > result_matrix[protein1, protein2]) {
    result_matrix[protein1, protein2] <- auc
    result_matrix[protein2, protein1] <- auc  # Symmetric matrix
    ventricle_matrix[protein1, protein2] <- ventricle
    ventricle_matrix[protein2, protein1] <- ventricle  # Symmetric matrix
  }
}

# Format the AUC values with ventricle information for heatmap
result_matrix_with_ventricle1 <- matrix(
  paste(sprintf("%.3f", result_matrix), " (", ventricle_matrix, ")", sep=""),
  nrow = nrow(result_matrix),
  dimnames = dimnames(result_matrix)
)

# Function to analyze combination of two proteins
analyze_combination <- function(data, protein1, protein2) {
  model_data <- data %>%
    select(Subtype, !!sym(protein1), !!sym(protein2), Age, Sex)
  
  model_data$Subtype <- factor(model_data$Subtype)
  
  set.seed(123)
  k_folds <- 5
  preds <- numeric(nrow(model_data))
  actual <- numeric(nrow(model_data))
  folds <- createFolds(model_data$Subtype, k = k_folds)
  
  for (i in 1:k_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- folds[[i]]
    
    train_data <- model_data[train_indices, ]
    test_data <- model_data[test_indices, ]
    
    svm_model <- svm(
      Subtype ~ .,  
      data = train_data,
      kernel = "radial",
      probability = TRUE,
      cost = 1
    )
    
    pred_probs <- predict(svm_model, test_data, probability = TRUE)
    pred_attrs <- attr(pred_probs, "probabilities")
    
    preds[test_indices] <- pred_attrs[, "AD"]
    actual[test_indices] <- ifelse(test_data$Subtype == "AD", 1, 0)
  }
  
  roc_resultz <- roc(actual, preds)
  auc_value <- auc(roc_resultz)
  
  ci_resultz <- ci.auc(roc_resultz)
  p_value <- 2 * (1 - pnorm(abs((auc_value - 0.5) / ((ci_resultz[3] - ci_resultz[1]) / 3.92))))
  
  return(list(
    protein1 = protein1,
    protein2 = protein2,
    auc = auc_value,
    p_value = p_value,
    roc_obj = roc_resultz
  ))
}

# Perform analysis for all combinations of proteins within and between groups
resultz <- list()
counter <- 1

# Amyloid-Amyloid combinations
for (i in seq_along(amyloid_associated)) {
  for (j in i:length(amyloid_associated)) {
    protein1 <- amyloid_associated[i]
    protein2 <- amyloid_associated[j]
    cat("Analyzing:", protein1, "+", protein2, "(Amyloid-Amyloid)\n")
    resultz[[counter]] <- analyze_combination(data, protein1, protein2)
    counter <- counter + 1
  }
}

# Tau-Tau combinations
for (i in seq_along(tau_associated)) {
  for (j in i:length(tau_associated)) {
    protein1 <- tau_associated[i]
    protein2 <- tau_associated[j]
    cat("Analyzing:", protein1, "+", protein2, "(Tau-Tau)\n")
    resultz[[counter]] <- analyze_combination(data, protein1, protein2)
    counter <- counter + 1
  }
}

# Amyloid-Tau combinations
for (amyloid_protein in amyloid_associated) {
  for (tau_protein in tau_associated) {
    cat("Analyzing:", amyloid_protein, "+", tau_protein, "(Amyloid-Tau)\n")
    resultz[[counter]] <- analyze_combination(data, amyloid_protein, tau_protein)
    counter <- counter + 1
  }
}

# Convert results to dataframe
resultz_df <- do.call(rbind, lapply(resultz, function(x) {
  data.frame(
    Protein1 = x$protein1,
    Protein2 = x$protein2,
    AUC = x$auc,
    P_value = round(x$p_value, 10)
  )
}))

# Ensure p-values are numeric and handle NAs properly
resultz_df$P_value <- as.numeric(resultz_df$P_value)
resultz_df <- resultz_df %>% filter(!is.na(P_value))

# ---- Code for the second matrix (protein1 + protein2 only) ----

# Create a matrix for protein combinations without ventricle
result_matrix <- matrix(NA, nrow = length(all_proteins), ncol = length(all_proteins),
                        dimnames = list(all_proteins, all_proteins))

# Fill the matrix with AUC values from the protein-only analysis
for (i in 1:nrow(resultz_df)) {
  protein1 <- resultz_df$Protein1[i]
  protein2 <- resultz_df$Protein2[i]
  auc <- resultz_df$AUC[i]
  
  if (is.na(result_matrix[protein1, protein2]) || auc > result_matrix[protein1, protein2]) {
    result_matrix[protein1, protein2] <- auc
    result_matrix[protein2, protein1] <- auc  # Make the matrix symmetric
  }
}

# Format the AUC values for consistency with first matrix
formatted_result_matrix <- matrix(
  paste(sprintf("%.3f", result_matrix), sep=""),
  nrow = nrow(result_matrix),
  dimnames = dimnames(result_matrix)
)

# Ensure diagonal elements are set to NA, these will be colored in white in the resulting plots
diag(result_matrix_with_ventricle1) <- NA
diag(formatted_result_matrix) <- NA

# Assign to the variable expected by the heatmap code
result_matrix <- formatted_result_matrix

# Function to process matrix for heatmap
process_matrix <- function(matrix_data) {
  numeric_matrix <- apply(matrix_data, 2, function(x) as.numeric(sub("\\s*\\(.*\\)$", "", x)))
  rownames(numeric_matrix) <- rownames(matrix_data)
  colnames(numeric_matrix) <- colnames(matrix_data)
  
  # Set diagonal elements to NA (for white squares)
  # The aim is to produce Seaborn heatmaps
  diag(numeric_matrix) <- NA
  return(numeric_matrix)
}

# Process both matrices
numeric_matrix1 <- process_matrix(result_matrix_with_ventricle1)
numeric_matrix2 <- process_matrix(result_matrix)

# Define cluster origins
cluster_info <- data.frame(
  Protein = rownames(numeric_matrix1),
  Cluster = c(
    rep("Amyloid-associated", length(amyloid_associated)),  
    rep("Tau-associated", length(tau_associated))      
  ),
  stringsAsFactors = FALSE
)

# Reorder numeric_matrix rows and columns based on cluster information
sorted_proteins <- cluster_info$Protein[order(cluster_info$Cluster)]
numeric_matrix1 <- numeric_matrix1[sorted_proteins, sorted_proteins]
numeric_matrix2 <- numeric_matrix2[sorted_proteins, sorted_proteins]

# Reorder cluster origins accordingly
cluster_info <- cluster_info[order(cluster_info$Cluster), ]

# Define colors for each cluster
cluster_colors <- c("Amyloid-associated" = "orange", "Tau-associated" = "deeppink")

# Create row annotation
row_anno <- rowAnnotation(
  Cluster = cluster_info$Cluster,
  col = list(Cluster = cluster_colors),
  width = unit(1, "cm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  ),
  annotation_name_gp = gpar(fontsize = 14, fontface = "bold") # This changes the "Cluster" text size
)

# Define color palette
col_fun <- colorRamp2(c(0.4, 0.7, 1), c("blue", "darkblue", "black"))

# Create first heatmap, that includes ventricular volumes
ht1 <- Heatmap(numeric_matrix1,
               column_title = "AUC Heatmap - Combined Proteins + Ventricle",
               column_title_gp = gpar(fontsize = 16, fontface = "bold"),
               name = "AUC Value",
               col = col_fun,
               na_col = "white",
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 14),
               column_names_gp = gpar(fontsize = 14),
               width = unit(14, "cm"),
               height = unit(14, "cm"),
               heatmap_legend_param = list(
                 title = "AUC Value",
                 title_gp = gpar(fontsize = 14, fontface = "bold"),
                 labels_gp = gpar(fontsize = 14)
               ),
               left_annotation = row_anno
)

# Create second heatmap, that does not include ventricular volumes
ht2 <- Heatmap(numeric_matrix2,
               column_title = "AUC Heatmap - Combined Proteins",
               column_title_gp = gpar(fontsize = 16, fontface = "bold"),
               name = "AUC Value",
               col = col_fun,
               na_col = "white",
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 14),
               column_names_gp = gpar(fontsize = 14),
               width = unit(14, "cm"),
               height = unit(14, "cm"),
               heatmap_legend_param = list(
                 title = "AUC Value",
                 title_gp = gpar(fontsize = 14, fontface = "bold"),
                 labels_gp = gpar(fontsize = 14)
               )
)

pdf(
  file = "protein_heatmaps.pdf", 
  width = 16,       
  height = 11,      
  pointsize = 10, 
  paper = "special"
)

# Positioning the heatmaps side by side
draw(ht1 + ht2)
dev.off()
