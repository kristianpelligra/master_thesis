# Two-Protein and Ventricle Volume ROC Curve Analysis

# Set working directory
setwd("C:/Users/User/Desktop/RMT")

# Load required libraries
library(e1071)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(reshape2)

# Load and format data
data <- read_delim("heatmap_data_norm.csv", delim = ";") %>% 
  rename(LLV = "Left_Lateral", RLV = "Right_Lateral", TV = "Third", FV = "Fourth", ToV = "Total_Ventricle") %>% 
  filter(Subtype %in% c("AD", "SCI"))

# Define protein groups
amyloid_associated <- c("NPTX2", "NPTXR", "NPTX1", "CHL1", "NCAN", "PTPRN2", "TMEM132D",
                        "CCK", "PAM", "CLSTN1", "OMG", "APLP1", "NRCAM", "VGF", "CDH8",
                        "KIAA1549L")
tau_associated <- c("ARPP21", "SNCB", "GAP43", "NRGN", "SPP1", "AQP4", "RPH3A", "SNCA",
                    "AMPH", "DDAH1", "PEBP1")

# Ventricle volumes to include
ventricle_volumes <- c("LLV", "RLV", "TV", "FV", "ToV")

# Enhanced analysis function incorporating ventricle volumes
analyze_combination <- function(data, amyloid_protein, tau_protein, ventricle_volume) {
  tryCatch({
    # Prepare model data
    model_data <- data %>%
      select(Subtype, !!sym(amyloid_protein), !!sym(tau_protein), 
             !!sym(ventricle_volume), Age, Sex)
    
    # Convert Subtype to factor
    model_data$Subtype <- factor(model_data$Subtype)
    
    # Ensure sufficient samples in each class
    if (length(unique(model_data$Subtype)) < 2 || 
        min(table(model_data$Subtype)) < 5) {
      warning(paste("Insufficient samples for", amyloid_protein, "+", 
                    tau_protein, "+", ventricle_volume))
      return(NULL)
    }
    
    # Cross-validation setup
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
      
      # Handle prediction probabilities
      preds[test_indices] <- if ("AD" %in% colnames(pred_attrs)) {
        pred_attrs[, "AD"]
      } else {
        0
      }
      actual[test_indices] <- ifelse(test_data$Subtype == "AD", 1, 0)
    }
    
    # ROC and statistical analysis
    roc_result <- roc(actual, preds)
    auc_value <- auc(roc_result)
    
    ci_result <- ci.auc(roc_result)
    p_value <- if (!is.na(ci_result[3]) && !is.na(ci_result[1])) {
      2 * (1 - pnorm(abs((auc_value - 0.5) / ((ci_result[3] - ci_result[1]) / 3.92))))
    } else {
      NA
    }
    
    return(list(
      amyloid_protein = amyloid_protein,
      tau_protein = tau_protein,
      ventricle_volume = ventricle_volume,
      auc = auc_value,
      p_value = p_value,
      roc_obj = roc_result
    ))
  }, 
  error = function(e) {
    warning(paste("Error in analysis for", amyloid_protein, "+", 
                  tau_protein, "+", ventricle_volume, ":", e$message))
    return(NULL)
  })
}

# Perform analysis for all combinations
results <- list()
for (amyloid_protein in amyloid_associated) {
  for (tau_protein in tau_associated) {
    for (ventricle_volume in ventricle_volumes) {
      result <- analyze_combination(data, amyloid_protein, tau_protein, ventricle_volume)
      if (!is.null(result)) {
        results <- append(results, list(result))
      }
    }
  }
}

# Convert results to dataframe
results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    AmyloidProtein = x$amyloid_protein,
    TauProtein = x$tau_protein,
    VentricleVolume = x$ventricle_volume,
    AUC = x$auc,
    P_value = round(x$p_value, 10)
  )
}))

# Data processing and statistical analysis
results_df <- results_df %>% filter(!is.na(P_value))
results_df$Adjusted_P <- round(p.adjust(results_df$P_value, method = "fdr"), 30)

# Sort and filter results
results_df <- results_df %>% arrange(desc(AUC))
significant_results <- results_df %>% filter(Adjusted_P < 0.05)

# Reporting
if (nrow(significant_results) > 0) {
  print("Significant Protein and Ventricle Combinations:")
  print(significant_results)
  
  # Visualize top ROC curves
  pdf("top_protein_ventricle_roc_curves.pdf", width = 10, height = 8)
  
  # Find indices of top combinations
  top_indices <- order(results_df$AUC, decreasing = TRUE)[1:min(5, nrow(results_df))]
  
  # Plot the best combination
  plot(results[[match(results_df$AUC[top_indices[1]], sapply(results, function(x) x$auc))]]$roc_obj, 
       main = "ROC Curves for Top Protein-Ventricle Combinations",
       col = 1, lwd = 2)
  
  # Add next top 4 combinations
  for (i in 2:length(top_indices)) {
    plot(results[[match(results_df$AUC[top_indices[i]], sapply(results, function(x) x$auc))]]$roc_obj, 
         add = TRUE, col = i, lwd = 2)
  }
  
  # Create legend
  legend_text <- sapply(top_indices, function(i) {
    paste0(results_df[i, "AmyloidProtein"], " + ", results_df[i, "TauProtein"], " + ", 
           results_df[i, "VentricleVolume"], " (AUC = ", round(results_df[i, "AUC"], 3), ")")
  })
  
  legend("bottomright", legend = legend_text, col = 1:length(top_indices), lwd = 2)
  
  dev.off()  # Close PDF only if it was created
} else {
  print("No significant results found (Adjusted P-value < 0.05).")
}

