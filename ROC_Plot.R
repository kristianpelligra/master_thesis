# ROC Curve Analysis for Specific Combinations
# Run Single_SVM, Single_SVM_Ventricle, Combined_SVM and Combined_SVM_Ventricle first

# Set working directory
setwd("C:/Users/User/Desktop/RMT")

# Load required libraries
library(e1071)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(readr)

# Load and format data
data <- read_delim("heatmap_data_norm.csv", delim = ";") %>% 
  rename(LLV = "Left_Lateral", RLV = "Right_Lateral", TV = "Third", FV = "Fourth", ToV = "Total_Ventricle") %>% 
  filter(Subtype %in% c("AD", "SCI"))

# Function to analyze specific combinations
analyze_combination <- function(data, columns) {
  # Prepare model data
  model_data <- data %>%
    select(Subtype, all_of(columns), Age, Sex)
  
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
    train_indices <- unlist(folds[-i])
    test_indices <- folds[[i]]
    
    train_data <- model_data[train_indices, ]
    test_data <- model_data[test_indices, ]
    
    # Train SVM model
    svm_model <- svm(
      Subtype ~ .,  
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
  
  return(list(
    columns = columns,
    auc = auc_value,
    roc_obj = roc_result
  ))
}

# Define combinations to analyze
combinations <- list(
  c("PTPRN2", "GAP43"),
  c("PTPRN2", "GAP43", "ToV"),
  c("GAP43", "ToV"),
  c("PTPRN2", "ToV"),
  c("GAP43"),
  c("PTPRN2")
)

# Custom color palette
colors <- c(
  "green",  # Deep Green
  "darkgreen",  # Bright Green
  "darkblue",  # Dark Blue
  "darkred",  # Brick Red
  "blue",  # Blue
  "red"   # Red
)

# Perform analysis for each combination
results <- lapply(combinations, function(cols) {
  analyze_combination(data, cols)
})

# Save ROC curves as a PNG
png("specific_combinations_roc_curves.png", width = 1200, height = 1000, res = 150)

# Increase text size for all elements
par(cex.main = 1.2,    # Title size
    cex.lab = 1.2,     # Axis label size
    cex.axis = 1.2,    # Axis tick mark label size
    cex = 1)       # General text size

# Plot the first combination
plot(results[[1]]$roc_obj, 
     main = "ROC Plot",
     col = colors[1], lwd = 4.5)

# Add other combinations
for (i in 2:length(results)) {
  plot(results[[i]]$roc_obj, add = TRUE, col = colors[i], lwd = 4.5)
}

# Create legend with increased text size
legend_text <- sapply(results, function(x) {
  paste0(paste(x$columns, collapse = " + "), " (AUC = ", round(x$auc, 3), ")")
})
legend("bottomright", 
       legend = legend_text, 
       col = colors[1:length(results)], 
       lwd = 4.5,
       cex = 1)  # Increase legend text size

# Close the PNG device
dev.off()

# Print AUC values for reference
for (i in seq_along(results)) {
  cat("Combination:", paste(results[[i]]$columns, collapse = " + "), 
      "- AUC:", round(results[[i]]$auc, 3), "\n")
}

