# Comprehensive ROC Curve Analysis for Proteins

# Set working directory (adjust to your specific path)
setwd("C:/Users/User/Desktop/RMT")

# Load required libraries
library(e1071)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)

# Load and format data
data <- read_delim("heatmap_data_norm.csv", delim = ";") %>% 
  filter(Subtype %in% c("AD", "SCI")) %>% 
  select(2:3, 7:60, 61)

# Function to analyze one protein as a predictor
analyze_protein <- function(data, protein_col) {
  # Prepare model data
  model_data <- data %>%
    select(Subtype, !!sym(protein_col), Age, Sex)
  
  # Convert Subtype to factor
  model_data$Subtype <- factor(model_data$Subtype)
  
  # Set up cross-validation
  set.seed(123)
  k_folds <- 5
  preds <- numeric(nrow(model_data))
  actual <- numeric(nrow(model_data))
  folds <- createFolds(model_data$Subtype, k = k_folds)
  
  # Perform k-fold cross-validation
  for (i in 1:k_folds) {
    # Split data into training and testing sets
    train_indices <- unlist(folds[-i])
    test_indices <- folds[[i]]
    
    train_data <- model_data[train_indices, ]
    test_data <- model_data[test_indices, ]
    
    # Train SVM model
    svm_model <- svm(
      Subtype ~ .,  # Use all predictors
      data = train_data,
      kernel = "radial",
      probability = TRUE,
      cost = 1
    )
    
    # Predict probabilities
    pred_probs <- predict(svm_model, test_data, probability = TRUE)
    pred_attrs <- attr(pred_probs, "probabilities")
    
    # Store predictions and actual values
    preds[test_indices] <- pred_attrs[, "AD"]
    actual[test_indices] <- ifelse(test_data$Subtype == "AD", 1, 0)
  }
  
  # Calculate ROC and AUC
  roc_result <- roc(actual, preds)
  auc_value <- auc(roc_result)
  
  # Calculate confidence interval and p-value
  ci_result <- ci.auc(roc_result)
  p_value <- 2 * (1 - pnorm(abs((auc_value - 0.5) / ((ci_result[3] - ci_result[1]) / 3.92))))
  
  return(list(
    protein = protein_col,
    auc = auc_value,
    p_value = p_value,
    roc_obj = roc_result
  ))
}

# Prepare protein columns
protein_cols <- colnames(data)[3:51]

# Initialize results storage
results <- list()
counter <- 1

# Perform analysis for all proteins
for (protein in protein_cols) {
  cat("Analyzing:", protein, "\n")
  results[[counter]] <- analyze_protein(data, protein)
  counter <- counter + 1
}

# Convert results to dataframe
results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    protein = x$protein,
    auc = x$auc,
    p_value = round(x$p_value, 7)
  )
}))

# Sort results by AUC
results_df <- results_df %>% arrange(desc(auc))
print("Top 10 Proteins:")
print(head(results_df, 10))

# Ensure p-values are numeric and handle NA values
results_df$p_value <- as.numeric(results_df$p_value)
results_df <- results_df %>% filter(!is.na(p_value))

# Perform p-value adjustment using False Discovery Rate (FDR)
results_df$adjusted_p <- round(p.adjust(results_df$p_value, method = "fdr"), 12)

# Filter significant results
significant <- results_df %>% filter(adjusted_p < 0.05)
print("\nSignificant Proteins:")
print(significant)

# Plot ROC curves for top proteins
pdf("top_5_protein_roc_curves.pdf", width = 10, height = 8)
par(mfrow = c(1,1))

# Plot the best protein
best_protein_index <- which.max(sapply(results, function(x) x$auc))
plot(results[[best_protein_index]]$roc_obj, 
     main = "ROC Curves for Top Proteins",
     col = 1, lwd = 2)

# Add next top 4 proteins
top_indices <- order(sapply(results, function(x) x$auc), decreasing = TRUE)[2:5]
for (i in seq_along(top_indices)) {
  plot(results[[top_indices[i]]]$roc_obj, add = TRUE, col = i+1, lwd = 2)
}

# Create legend
legend_text <- c(
  with(results[[best_protein_index]], 
       paste0(protein, " (AUC = ", round(auc, 3), ")")),
  sapply(results[top_indices], function(x) {
    paste0(x$protein, " (AUC = ", round(x$auc, 3), ")")
  })
)
legend("bottomright", legend = legend_text, col = 1:5, lwd = 2)
dev.off()

cat("\n\nAnalysis complete. Results saved to disk.\n")
