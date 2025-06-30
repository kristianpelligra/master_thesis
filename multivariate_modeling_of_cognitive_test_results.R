# Cognitive Test Modeling Script with Multiple Regression
# Using Specific Amyloid and Tau-Associated Proteins and Visualization

library(tidyverse)
library(broom)
library(boot)
library(ggplot2)
library(patchwork)

# Define Protein Groups 
amyloid_associated <- c(
  "NPTX2", "NPTXR", "NPTX1", "CHL1", "NCAN", "PTPRN2", "TMEM132D",
  "CCK", "PAM", "CLSTN1", "OMG", "APLP1", "NRCAM", "VGF", "CDH8", "KIAA1549L"
)

tau_associated <- c(
  "ARPP21", "SNCB", "GAP43", "NRGN", "SPP1", "AQP4", "RPH3A", "SNCA",
  "AMPH", "DDAH1", "PEBP1"
)

cognitive_tests <- c("RAVLT", "MMSE", "MoCA")

# Cross-Validation Function
cv_r2 <- function(formula, data, k = 5) {
  set.seed(123)
  data <- drop_na(data)
  
  cv_results <- boot(data, statistic = function(data, indices) {
    model <- lm(formula, data = data[indices, ])
    summary(model)$r.squared
  }, R = k)
  
  c(mean = mean(cv_results$t, na.rm = TRUE), sd = sd(cv_results$t, na.rm = TRUE))
}

# Model Statistics Function
compute_model_stats <- function(formula, data) {
  data <- drop_na(data)
  if (nrow(data) <= 10) {
    return(list(r2 = NA, adj_r2 = NA, f_stat = NA, p_value = NA,
                predictors = NA, sample_size = NA))
  }
  model <- lm(formula, data = data)
  summary_stats <- summary(model)
  
  list(
    r2 = summary_stats$r.squared,
    adj_r2 = summary_stats$adj.r.squared,
    f_stat = summary_stats$fstatistic[1],
    p_value = pf(summary_stats$fstatistic[1], summary_stats$fstatistic[2],
                 summary_stats$fstatistic[3], lower.tail = FALSE),
    predictors = length(coef(model)) - 1,
    sample_size = nrow(data)
  )
}

# Model Analysis Function
analyze_cognitive_models <- function(data, cognitive_test) {
  protein_cols <- c(amyloid_associated, tau_associated)
  ventricle_cols <- c("Left_Lateral", "Right_Lateral", "Third", "Fourth")
  
  results_df <- tibble(
    Model_Type = character(), Formula = character(),
    CV_R2_Mean = numeric(), CV_R2_SD = numeric(),
    R2 = numeric(), Adj_R2 = numeric(),
    F_Statistic = numeric(), P_Value = numeric(),
    Num_Predictors = numeric(), Sample_Size = numeric()
  )
  
  # Model Building Routines
  build_models <- function(formula, model_data, type) {
    if (nrow(model_data) > 10) {
      cv_result <- cv_r2(formula, model_data)
      model_stats <- compute_model_stats(formula, model_data)
      results_df <<- add_row(results_df,
                             Model_Type = type,
                             Formula = as.character(formula),
                             CV_R2_Mean = cv_result["mean"],
                             CV_R2_SD = cv_result["sd"],
                             R2 = model_stats$r2,
                             Adj_R2 = model_stats$adj_r2,
                             F_Statistic = model_stats$f_stat,
                             P_Value = model_stats$p_value,
                             Num_Predictors = model_stats$predictors,
                             Sample_Size = model_stats$sample_size
      )
    }
  }
  
  # Single Protein Models
  cat("Analyzing Protein Only models...\n")
  for (protein in protein_cols) {
    f <- as.formula(paste(cognitive_test, "~", protein, "+ Age + Sex"))
    df <- select(data, all_of(c(cognitive_test, protein, "Age", "Sex"))) %>% drop_na()
    build_models(f, df, "Protein Only")
  }
  
  # Protein + Ventricle Models
  cat("Analyzing Protein + Ventricle models...\n")
  for (protein in protein_cols) {
    for (ventricle in ventricle_cols) {
      f <- as.formula(paste(cognitive_test, "~", protein, "+", ventricle, "+ Age + Sex"))
      df <- select(data, all_of(c(cognitive_test, protein, ventricle, "Age", "Sex"))) %>% drop_na()
      build_models(f, df, "Protein + Ventricle")
    }
  }
  
  # Combined Protein Pairs
  cat("Analyzing Combined Proteins models...\n")
  make_pairs <- function(group) {
    combn(group, 2, simplify = FALSE)
  }
  protein_pairs <- c(make_pairs(amyloid_associated), make_pairs(tau_associated))
  protein_pairs <- c(protein_pairs, expand.grid(amyloid_associated, tau_associated, stringsAsFactors = FALSE) %>% 
                       split(1:nrow(.)))
  
  for (pair in protein_pairs) {
    f <- as.formula(paste(cognitive_test, "~", pair[[1]], "+", pair[[2]], "+ Age + Sex"))
    df <- select(data, all_of(c(cognitive_test, pair[[1]], pair[[2]], "Age", "Sex"))) %>% drop_na()
    build_models(f, df, "Combined Proteins")
  }
  
  # Combined Proteins + Ventricle
  cat("Analyzing Combined Proteins + Ventricle models...\n")
  for (pair in protein_pairs) {
    for (ventricle in ventricle_cols) {
      f <- as.formula(paste(cognitive_test, "~", pair[[1]], "+", pair[[2]], "+", ventricle, "+ Age + Sex"))
      df <- select(data, all_of(c(cognitive_test, pair[[1]], pair[[2]], ventricle, "Age", "Sex"))) %>% drop_na()
      build_models(f, df, "Combined Proteins + Ventricle")
    }
  }
  
  summary_results <- results_df %>%
    group_by(Model_Type) %>%
    summarize(
      Cognitive_Test = cognitive_test,
      Mean_CV_R2 = mean(CV_R2_Mean, na.rm = TRUE),
      Median_CV_R2 = median(CV_R2_Mean, na.rm = TRUE),
      SD_CV_R2 = sd(CV_R2_Mean, na.rm = TRUE),
      Mean_R2 = mean(R2, na.rm = TRUE),
      Mean_Adj_R2 = mean(Adj_R2, na.rm = TRUE),
      Mean_P_Value = mean(P_Value, na.rm = TRUE),
      Median_P_Value = median(P_Value, na.rm = TRUE),
      Mean_log10_P = mean(-log10(P_Value), na.rm = TRUE),
      Median_log10_P = median(-log10(P_Value), na.rm = TRUE),
      Num_Models = n(),
      Mean_Predictors = mean(Num_Predictors, na.rm = TRUE),
      Mean_Sample_Size = mean(Sample_Size, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(detailed_results = results_df, summary_results = summary_results)
}

# Main Execution
main <- function() {
  if (!file.exists("heatmap_data_norm_cogn.csv")) {
    stop("Error: Data file 'heatmap_data_norm_cogn.csv' not found.")
  }
  
  cat("Loading data file: heatmap_data_norm_cogn.csv\n")
  data <- read_delim("heatmap_data_norm_cogn.csv", delim = ";")
  
  all_results <- list()
  all_summaries <- list()
  
  for (test in cognitive_tests) {
    message(paste("Processing", test))
    results <- analyze_cognitive_models(data, test)
    all_results[[test]] <- results$detailed_results
    all_summaries[[test]] <- results$summary_results
  }
  
  final_detailed <- bind_rows(all_results) %>% mutate(Model_ID = row_number())
  write_csv(final_detailed, "cognitive_models_detailed.csv")
  
  final_summary <- bind_rows(all_summaries)
  write_csv(final_summary, "cognitive_models_summary.csv")
  
  png("cognitive_models_comparison.png", width = 12, height = 10, units = "in", res = 300)
  plot <- create_plots(all_summaries, final_detailed)
  print(plot)
  dev.off()
  
  cat("Plot saved as cognitive_models_comparison.png\n")
  
  top_models <- final_detailed %>%
    group_by(Model_Type) %>%
    arrange(desc(CV_R2_Mean)) %>%
    select(Model_Type, Formula, CV_R2_Mean, P_Value, Num_Predictors, Sample_Size)
  
  cat("Top performing models by cross-validated R²:\n")
  print(top_models)
  write_csv(top_models, "top_cognitive_models.csv")
  
  list(
    detailed_results = final_detailed,
    summary_results = final_summary,
    top_models = top_models
  )
}

# Run Script
results <- main()
print(results$summary_results)

cat("\nAmyloid-associated proteins used:\n")
print(amyloid_associated)

cat("\nTau-associated proteins used:\n")
print(tau_associated)
