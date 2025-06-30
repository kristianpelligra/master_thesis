# Correlation Plots for Cognitive Test Models
# Creates three aligned plots showing model predictions vs. actual values for 
# MMSE, MoCA, and RAVLT using two models:
# 1. PTPRN2 + NRGN
# 2. PTPRN2 + NRGN + Total_Ventricle

library(tidyverse)
library(cowplot)  # For plot alignment
library(ggpubr)   # For correlation statistics on plot

# Function to create correlation plot for a specific cognitive test and model
create_correlation_plot <- function(data, cog_test, model_formula, model_name, color) {
  # Remove NA values for the variables in the model
  # First extract all variables from the formula
  all_vars <- all.vars(as.formula(model_formula))
  
  # Create complete cases dataset for this specific model
  model_data <- data %>% 
    select(all_of(all_vars)) %>%
    drop_na()
  
  # Report sample size for this model
  cat(sprintf("  %s - %s: Using %d complete observations\n", 
              cog_test, model_name, nrow(model_data)))
  
  # Skip if too few observations
  if (nrow(model_data) < 10) {
    stop("Too few complete observations (less than 10)")
  }
  
  # Fit the model on complete cases only
  model <- lm(model_formula, data = model_data)
  
  # Get predicted values for these same complete cases
  predicted <- predict(model)
  
  # Create data frame for plotting - now both have same number of rows
  plot_data <- data.frame(
    Actual = model_data[[cog_test]],
    Predicted = predicted,
    Model = model_name
  )
  
  # Calculate model stats to display
  r_value <- cor(plot_data$Actual, plot_data$Predicted, use = "complete.obs")
  r_squared <- summary(model)$r.squared
  p_value <- cor.test(plot_data$Actual, plot_data$Predicted)$p.value
  n_obs <- nrow(plot_data)
  
  # Create line of equality (y=x)
  min_val <- min(c(plot_data$Actual, plot_data$Predicted), na.rm = TRUE)
  max_val <- max(c(plot_data$Actual, plot_data$Predicted), na.rm = TRUE)
  perfect_line <- data.frame(
    x = c(min_val, max_val),
    y = c(min_val, max_val)
  )
  
  # Create scatter plot with correlation line
  plot <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
    # Add perfect prediction line (y=x)
    geom_line(data = perfect_line, aes(x = x, y = y), 
              linetype = "dotted", color = "gray50", size = 0.7) +
    # Add regression line
    geom_smooth(method = "lm", color = "black", linetype = "dashed", 
                size = 0.5, fill = "lightgray") +
    # Add points
    geom_point(color = color, size = 3, alpha = 0.7) +
    # Add correlation stats
    annotate("text", x = min_val + (max_val-min_val)*0.05, 
             y = max_val - (max_val-min_val)*0.05,
             label = sprintf("R = %.2f\nR² = %.2f\np %s\nn = %d", 
                             r_value, r_squared, 
                             ifelse(p_value < 0.001, "< 0.001", 
                                    ifelse(p_value < 0.01, "< 0.01",
                                           ifelse(p_value < 0.05, "< 0.05", 
                                                  sprintf("= %.2f", p_value)))),
                             n_obs),
             hjust = 0, vjust = 1, size = 3.5) +
    # Labels and theme
    labs(
      title = cog_test,
      subtitle = model_name,
      x = paste("Actual", cog_test, "Score"),
      y = "Predicted Score"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  return(plot)
}

# Main function to generate all plots
generate_cognitive_correlation_plots <- function() {
  # Check if data file exists
  if (!file.exists("heatmap_data_norm_cogn.csv")) {
    stop("Error: Data file 'heatmap_data_norm_cogn.csv' not found.")
  }
  
  # Load data
  cat("Loading data...\n")
  data <- read_delim("heatmap_data_norm_cogn.csv", delim = ";")
  
  # Print summary of missing values for key variables
  cat("\nMissing value summary:\n")
  for (col in c("MMSE", "MoCA", "RAVLT", "PTPRN2", "NRGN", "Left_Lateral", "Right_Lateral", "Third", "Fourth", "Age", "Sex")) {
    if (col %in% names(data)) {
      missing <- sum(is.na(data[[col]]))
      total <- nrow(data)
      cat(sprintf("%s: %d missing out of %d (%.1f%%)\n", col, missing, total, missing/total*100))
    } else {
      cat(sprintf("%s: Column not found in data\n", col))
    }
  }
  
  # Add Total_Ventricle column if needed
  if (!"Total_Ventricle" %in% colnames(data)) {
    data <- data %>%
      mutate(Total_Ventricle = Left_Lateral + Right_Lateral + Third + Fourth)
    cat("\nCreated Total_Ventricle column\n")
  }
  
  # Cognitive tests to analyze
  cognitive_tests <- c("MMSE", "MoCA", "RAVLT")
  
  # Model specifications
  models <- list(
    model1 = list(
      formula_str = "~ PTPRN2 + NRGN + Age + Sex",
      name = "PTPRN2 + NRGN",
      color = "#1b9e77"  # Dark green
    ),
    model2 = list(
      formula_str = "~ PTPRN2 + NRGN + Total_Ventricle + Age + Sex",
      name = "PTPRN2 + NRGN + ToV",
      color = "#d95f02"  # Orange
    )
  )
  
  # Initialize list to store plots
  all_plots <- list()
  
  cat("\nCreating plots for each cognitive test and model:\n")
  
  # Create all plots
  for (test in cognitive_tests) {
    cat(sprintf("\nProcessing %s:\n", test))
    test_plots <- list()
    
    for (i in seq_along(models)) {
      model_info <- models[[i]]
      formula_str <- paste(test, model_info$formula_str)
      formula <- as.formula(formula_str)
      
      # Create plot with error handling
      tryCatch({
        plot <- create_correlation_plot(
          data = data,
          cog_test = test,
          model_formula = formula,
          model_name = model_info$name,
          color = model_info$color
        )
        
        test_plots[[i]] <- plot
        
      }, error = function(e) {
        # Print error message and create an empty plot as placeholder
        cat(sprintf("  Error creating plot for %s with model %s: %s\n", 
                    test, model_info$name, e$message))
        
        empty_plot <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = sprintf("Insufficient data for\n%s\n%s", test, model_info$name)) + 
          theme_void() +
          theme(plot.background = element_rect(fill = "white"))
        
        test_plots[[i]] <- empty_plot
      })
    }
    
    # Combine plots for this test
    combined <- plot_grid(
      test_plots[[1]] + theme(plot.subtitle = element_text(color = models$model1$color)),
      test_plots[[2]] + theme(plot.subtitle = element_text(color = models$model2$color)),
      ncol = 1,
      align = "v",
      labels = NULL
    )
    
    all_plots[[test]] <- combined
  }
  
  # Arrange all cognitive test plots horizontally
  cat("\nCombining all plots...\n")
  final_plot <- plot_grid(
    all_plots$MMSE, all_plots$MoCA, all_plots$RAVLT,
    ncol = 3,
    labels = c("A", "B", "C"),
    label_size = 18
  )
  
  # Add title to combined plot with improved visibility
  title_text <- "Cognitive Prediction by Top Protein Pair: Effect of Ventricular Volume"
  
  final_titled <- plot_grid(
    ggdraw() + 
      draw_label(
        title_text,
        fontface = "bold",
        size = 20,
        colour = "black"  # Ensure black text on white background
      ) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(15, 0, 5, 0)  # Add some padding
      ),
    final_plot,
    ncol = 1,
    rel_heights = c(0.07, 1)  # Slightly increase title space
  )
  
  # Save the plot
  output_file <- "cognitive_correlation_plots.png"
  ggsave(output_file, final_titled, width = 12.5, height = 8.33, dpi = 300)
  
  cat(sprintf("\nPlots saved as '%s'\n", output_file))
  
  return(final_titled)
}

# Execute the function to generate all plots
cat("Starting cognitive correlation plot generation\n")
correlation_plots <- generate_cognitive_correlation_plots()
print(correlation_plots)