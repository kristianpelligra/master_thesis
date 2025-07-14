# Load required libraries
library(e1071)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(reshape2)

# Fix namespace conflict by explicitly using dplyr::select
select <- dplyr::select

data <- read_delim("heatmap_data_norm.csv", delim = ";")

# Print script starting message
cat("Starting AD vs. SCI classification analysis with ALL protein pairs...\n")

# Add error handling for data loading
tryCatch({
  # Load data
  cat("Loading data...\n")
  data <- read_delim("heatmap_data_norm.csv", delim = ";", show_col_types = FALSE) %>% 
    rename(LLV = "Left_Lateral", RLV = "Right_Lateral", TV = "Third", 
           FV = "Fourth", ToV = "Total_Ventricle") %>% 
    filter(Subtype %in% c("AD", "SCI"))
  
  cat("Data loaded successfully with", nrow(data), "rows and", ncol(data), "columns\n")
  
  # Check for required columns
  required_cols <- c("Subtype", "Age", "Sex")
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse=", ")))
  }
  
  # Make sure we have both AD and SCI cases
  subtype_counts <- table(data$Subtype)
  cat("Subtype distribution:", paste(names(subtype_counts), subtype_counts, sep="=", collapse=", "), "\n")
  
  if (length(subtype_counts) != 2 || !all(c("AD", "SCI") %in% names(subtype_counts))) {
    stop("Data must contain both AD and SCI subtypes")
  }
}, error = function(e) {
  stop(paste("Error loading data:", e$message))
})

# Get all protein names (columns 7:55)
protein_names <- colnames(data)[7:55]

# Protein groups
amyloid_associated = c("NPTX2", "NPTXR", "NPTX1", "CHL1", "NCAN", "PTPRN2", "TMEM132D",
                       "CCK", "PAM", "CLSTN1", "OMG", "APLP1", "NRCAM", "VGF", "CDH8",
                       "KIAA1549L")
tau_associated = c("ARPP21", "SNCB", "GAP43", "NRGN", "SPP1", "AQP4", "RPH3A", "SNCA",
                   "AMPH", "DDAH1", "PEBP1")

# Ensure all protein names exist in the data
amyloid_associated <- intersect(amyloid_associated, colnames(data))
tau_associated <- intersect(tau_associated, colnames(data))

if (length(amyloid_associated) == 0) {
  stop("No amyloid-associated proteins found in the data!")
}

if (length(tau_associated) == 0) {
  stop("No tau-associated proteins found in the data!")
}

# Define ventricle columns
ventricle_cols <- c("LLV", "RLV", "TV", "FV", "ToV")

# Make sure ventricle columns exist
ventricle_cols <- intersect(ventricle_cols, colnames(data))

if (length(ventricle_cols) == 0) {
  stop("No ventricle columns found in the data!")
}

# Function to create CV folds
create_folds <- function(data) {
  set.seed(123)
  k_folds <- 5
  folds <- createFolds(data$Subtype, k = k_folds)
  return(folds)
}

# Function to run SVM CV
run_svm_cv <- function(data, formula) {
  tryCatch({
    data$Subtype <- factor(data$Subtype)
    
    k_folds <- 5
    preds <- numeric(nrow(data))
    actual <- numeric(nrow(data))
    folds <- create_folds(data)
    
    for (i in 1:k_folds) {
      train_indices <- unlist(folds[-i])
      test_indices <- folds[[i]]
      
      train_data <- data[train_indices, ]
      test_data <- data[test_indices, ]
      
      svm_model <- svm(
        formula,  
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
    
    return(auc_value)
  }, error = function(e) {
    warning(paste("Error in SVM model:", e$message))
    return(NA)
  })
}

# Create lists to store AUC values
results_list <- list(
  "Single Protein" = numeric(),
  "Protein + Ventricle" = numeric(),
  "Combined Proteins" = numeric(),
  "Combined Proteins + Ventricle" = numeric()
)

# 1. Single Protein models
cat("Analyzing Single Protein models...\n")
for (protein in c(amyloid_associated, tau_associated)) {
  formula <- as.formula(paste("Subtype ~", protein, "+ Age + Sex"))
  model_data <- data %>% dplyr::select(Subtype, !!sym(protein), Age, Sex)
  auc <- run_svm_cv(model_data, formula)
  results_list[["Single Protein"]] <- c(results_list[["Single Protein"]], auc)
}

# 2. Protein + Ventricle models
cat("Analyzing Protein + Ventricle models...\n")
for (protein in c(amyloid_associated, tau_associated)) {
  for (ventricle in ventricle_cols) {
    formula <- as.formula(paste("Subtype ~", protein, "+", ventricle, "+ Age + Sex"))
    model_data <- data %>% dplyr::select(Subtype, !!sym(protein), !!sym(ventricle), Age, Sex)
    auc <- run_svm_cv(model_data, formula)
    results_list[["Protein + Ventricle"]] <- c(results_list[["Protein + Ventricle"]], auc)
  }
}

# 3. Combined Proteins models
cat("Analyzing Combined Proteins models...\n")
protein_pairs <- list()
counter <- 0

# Create all protein pairs
for (i in 1:(length(amyloid_associated)-1)) {
  for (j in (i+1):length(amyloid_associated)) {
    counter <- counter + 1
    protein_pairs[[counter]] <- c(amyloid_associated[i], amyloid_associated[j])
  }
}

for (i in 1:(length(tau_associated)-1)) {
  for (j in (i+1):length(tau_associated)) {
    counter <- counter + 1
    protein_pairs[[counter]] <- c(tau_associated[i], tau_associated[j])
  }
}

for (i in 1:length(amyloid_associated)) {
  for (j in 1:length(tau_associated)) {
    counter <- counter + 1
    protein_pairs[[counter]] <- c(amyloid_associated[i], tau_associated[j])
  }
}

# Run models for all pairs
for (pair in protein_pairs) {
  protein1 <- pair[1]
  protein2 <- pair[2]
  formula <- as.formula(paste("Subtype ~", protein1, "+", protein2, "+ Age + Sex"))
  model_data <- data %>% dplyr::select(Subtype, !!sym(protein1), !!sym(protein2), Age, Sex)
  auc <- run_svm_cv(model_data, formula)
  results_list[["Combined Proteins"]] <- c(results_list[["Combined Proteins"]], auc)
}

# 4. Combined Proteins + Ventricle models
cat("Analyzing Combined Proteins + Ventricle models...\n")
for (pair in protein_pairs) {
  for (ventricle in ventricle_cols) {
    protein1 <- pair[1]
    protein2 <- pair[2]
    formula <- as.formula(paste("Subtype ~", protein1, "+", protein2, "+", ventricle, "+ Age + Sex"))
    model_data <- data %>% dplyr::select(Subtype, !!sym(protein1), !!sym(protein2), !!sym(ventricle), Age, Sex)
    auc <- run_svm_cv(model_data, formula)
    results_list[["Combined Proteins + Ventricle"]] <- c(results_list[["Combined Proteins + Ventricle"]], auc)
  }
}

# Filter out NA values
for (model_type in names(results_list)) {
  results_list[[model_type]] <- results_list[[model_type]][!is.na(results_list[[model_type]])]
}

# Convert results to data frame
results_df <- bind_rows(lapply(names(results_list), function(model_type) {
  data.frame(
    ModelType = model_type,
    AUC = results_list[[model_type]]
  )
}))

# Calculate mean AUC for each model type - Now with 2 decimal places
mean_aucs <- results_df %>%
  group_by(ModelType) %>%
  summarize(MeanAUC = mean(AUC, na.rm = TRUE),
            Count = n()) %>%
  mutate(Label = sprintf("%s (Mean AUC: %.2f, n=%d)", ModelType, MeanAUC, Count)) %>%
  arrange(desc(MeanAUC))

# Reorder factor levels based on mean AUC
results_df$ModelType <- factor(results_df$ModelType, 
                               levels = mean_aucs$ModelType,
                               ordered = TRUE)

# Set color palette
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")

# Create the density plot with improved legend settings
p <- ggplot(results_df, aes(x = AUC, fill = ModelType)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = colors,
                    labels = mean_aucs$Label) +
  theme_minimal() +
  labs(
    title = "Density Plot of AUC Values for AD vs. SCI Classification Models",
    x = "AUC",
    y = "Density",
    fill = "Model Type"
  ) +
  theme_bw() +  # Apply a clean theme first
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "bottom",  # Force legend to bottom
    legend.box = "horizontal",
    legend.margin = margin(t = 10, r = 0, b = 0, l = 0),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  ) +
  scale_x_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1)) +
  # Very important: define the legend layout explicitly
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", 
                             title.hjust = 0.5))

# Save the PNG with white background
ggsave("density_plot.png", 
       plot = p, 
       width = 12, 
       height = 8, 
       dpi = 300, 
       bg = "white")

# Display the plot
print(p)