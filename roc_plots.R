# Load required libraries
library(e1071)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# Set working directory
setwd("C:/Users/User/Desktop/RMT")

# Load and format data
data <- read_delim("heatmap_data_norm.csv", delim = ";") %>% 
  rename(
    LLV = "Left_Lateral",
    RLV = "Right_Lateral",
    TV  = "Third",
    FV  = "Fourth",
    ToV = "Total_Ventricle"
  ) %>%
  filter(Subtype %in% c("AD", "SCI")) %>%
  mutate(ToP = rowSums(across(7:55), na.rm = TRUE))

# Function to analyze SVM + ROC per combination
analyze_combination <- function(data, columns) {
  model_data <- data %>%
    dplyr::select(Subtype, all_of(columns), Age, Sex)
  
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
    
    svm_model <- svm(Subtype ~ ., data = train_data, kernel = "radial", probability = TRUE, cost = 1)
    pred_probs <- predict(svm_model, test_data, probability = TRUE)
    pred_attrs <- attr(pred_probs, "probabilities")
    
    preds[test_indices] <- pred_attrs[, "AD"]
    actual[test_indices] <- ifelse(test_data$Subtype == "AD", 1, 0)
  }
  
  roc_result <- roc(actual, preds)
  auc_value <- auc(roc_result)
  
  return(list(columns = columns, auc = auc_value, roc_obj = roc_result))
}

# Combinations
combinations_plot1 <- list(
  c("PTPRN2", "GAP43"),
  c("PTPRN2", "GAP43", "ToV"),
  c("ToV"),
  c("ToP")
)
combinations_plot2 <- list(
  c("GAP43", "ToV"),
  c("PTPRN2", "ToV"),
  c("GAP43"),
  c("PTPRN2")
)
colors1 <- c("darkred", "red", "orange", "gold")
colors2 <- c("navy", "blue", "deepskyblue", "turquoise")

# Perform SVM & ROC analysis
results_plot1 <- lapply(combinations_plot1, function(cols) analyze_combination(data, cols))
results_plot2 <- lapply(combinations_plot2, function(cols) analyze_combination(data, cols))

create_ggroc_plot <- function(results, colors, plot_title) {
  ggroc_data_list <- lapply(results, function(r) ggroc(r$roc_obj)$data)
  
  p <- ggplot() +
    mapply(function(df, col, idx) {
      geom_line(data = df, aes(x = 1 - specificity, y = sensitivity, color = factor(idx)), size = 1.5)
    }, df = ggroc_data_list, col = colors, idx = seq_along(results), SIMPLIFY = FALSE) +
    scale_color_manual(
      values = colors,
      labels = sapply(results, function(x) {
        paste0(paste(x$columns, collapse = " + "), " (AUC = ", round(x$auc, 3), ")")
      }),
      name = NULL
    ) +
    labs(
      title = plot_title,
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = c(0.65, 0.2),
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.5, "lines"),
      panel.grid = element_blank(),             # No gridlines
      axis.line = element_line(color = "black", size = 0.8),  # Add black axis lines
      axis.ticks = element_line(color = "black", size = 0.5)  # Optional: black tick marks
    )
  return(p)
}



# Build individual ROC ggplots
plot1 <- create_ggroc_plot(results_plot1, colors1, "Combined Proteins ± ToV")
plot2 <- create_ggroc_plot(results_plot2, colors2, "Single Proteins ± ToV")

# Combine with patchwork (A/B tags)
final_plot <- (plot1 | plot2) +
  plot_annotation(
    title = "ROC Curves: Protein Combinations and Single Features",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")),
    tag_levels = 'A'
  )

# Save final figure
ggsave("ROC_Patchwork_Tagged.png", plot = final_plot, width = 12, height = 6, dpi = 300)

# Print AUC summaries
cat("=== ROC Plot 1: Proteins vs. Proteins + ToV ===\n")
for (i in seq_along(results_plot1)) {
  cat("Combination:", paste(results_plot1[[i]]$columns, collapse = " + "), "- AUC:", round(results_plot1[[i]]$auc, 3), "\n")
}

cat("\n=== ROC Plot 2: Single Features ===\n")
for (i in seq_along(results_plot2)) {
  cat("Combination:", paste(results_plot2[[i]]$columns, collapse = " + "), "- AUC:", round(results_plot2[[i]]$auc, 3), "\n")
}