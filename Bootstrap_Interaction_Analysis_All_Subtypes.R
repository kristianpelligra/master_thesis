# Load required libraries
library(tidyverse)
library(limma)
library(ggeffects)
library(ggplot2)
library(gridExtra)
library(rsample)  # For bootstrapping
library(broom)    # For extracting model summaries
library(vsn)

# Load and format data
d_df <- read_delim("heatmap_data.csv", delim = ";") %>%
  rename(
    Subtype = subtype, Sex = gender, Age = sampling_age, 
    A = a, T = t, 
    Left_Lateral = Left_Lateral_Ventricle, 
    Right_Lateral = Right_Lateral_Ventricle, 
    Third = Third_Ventricle, Fourth = Fourth_Ventricle, 
    Total_Ventricle = Total_Adjusted_CSF_Volume, Qalb = alb_ratio
  ) %>%
  mutate(across(c(Left_Lateral, Right_Lateral, Third, Fourth, Total_Ventricle, ICV), ~ . / 1000)) %>% # Normalize selected columns
  filter(Subtype %in% c("AD", "SCI", "MCI_A-T-", "MCI_A+T+")) %>%  # Filter for all four subtypes
  as.data.frame()

# Collecting protein columns
exp_data <- d_df[, 7:55]

# Quantile normalization
normalized_exp_data <- normalizeBetweenArrays(as.matrix(exp_data), method = "quantile")

# Variance stabilization normalization (VSN)
vsn_fit <- vsn2(as.matrix(normalized_exp_data))
vsn_corrected_df <- predict(vsn_fit, newdata = as.matrix(normalized_exp_data))
normData <- as.data.frame(vsn_corrected_df)  # Convert back to dataframe

# Replace processed data in the original dataframe
d_df <- cbind(d_df[, 1:6], normData, d_df[, 56:ncol(d_df)])

# Set up bootstrap resampling
set.seed(123)
bootstrap_samples <- bootstraps(d_df, times = 500)

# Convert Subtype and Sex to factors for modeling
d_df$Subtype <- as.factor(d_df$Subtype)
d_df$Sex <- as.factor(d_df$Sex)

# Function to create bootstrapped interaction analysis
create_bootstrapped_interaction_plot <- function(protein, ventricle) {
  # Bootstrapping function for model parameters
  bootstrap_model <- function(split) {
    # Extract bootstrap sample
    boot_sample <- analysis(split)
    
    # Full interaction model
    full_model <- lm(as.formula(paste(protein, "~", ventricle, "*Subtype + Sex + Age + ab42_40")), data = boot_sample)
    
    # Separate models for each subtype
    ad_data <- boot_sample[boot_sample$Subtype == "AD", ]
    sci_data <- boot_sample[boot_sample$Subtype == "SCI", ]
    mci_neg_data <- boot_sample[boot_sample$Subtype == "MCI_A-T-", ]
    mci_pos_data <- boot_sample[boot_sample$Subtype == "MCI_A+T+", ]
    
    ad_model <- lm(as.formula(paste(protein, "~", ventricle)), data = ad_data)
    sci_model <- lm(as.formula(paste(protein, "~", ventricle)), data = sci_data)
    mci_neg_model <- lm(as.formula(paste(protein, "~", ventricle)), data = mci_neg_data)
    mci_pos_model <- lm(as.formula(paste(protein, "~", ventricle)), data = mci_pos_data)
    
    # Return key statistics
    tibble(
      full_r_squared = summary(full_model)$r.squared,
      ad_slope = coef(ad_model)[2],
      ad_r_squared = summary(ad_model)$r.squared,
      sci_slope = coef(sci_model)[2],
      sci_r_squared = summary(sci_model)$r.squared,
      mci_neg_slope = coef(mci_neg_model)[2],
      mci_neg_r_squared = summary(mci_neg_model)$r.squared,
      mci_pos_slope = coef(mci_pos_model)[2],
      mci_pos_r_squared = summary(mci_pos_model)$r.squared
    )
  }
  
  # Apply bootstrapping
  bootstrap_results <- bootstrap_samples %>%
    mutate(model_stats = map(splits, bootstrap_model)) %>%
    unnest(model_stats)
  
  # Calculate bootstrap confidence intervals
  boot_summary <- bootstrap_results %>%
    summarise(
      full_r_squared = mean(full_r_squared, na.rm = TRUE),
      full_r_squared_ci_low = quantile(full_r_squared, 0.025, na.rm = TRUE),
      full_r_squared_ci_high = quantile(full_r_squared, 0.975, na.rm = TRUE),
      
      ad_slope = mean(ad_slope, na.rm = TRUE),
      ad_slope_ci_low = quantile(ad_slope, 0.025, na.rm = TRUE),
      ad_slope_ci_high = quantile(ad_slope, 0.975, na.rm = TRUE),
      ad_r_squared = mean(ad_r_squared, na.rm = TRUE),
      ad_r_squared_ci_low = quantile(ad_r_squared, 0.025, na.rm = TRUE),
      ad_r_squared_ci_high = quantile(ad_r_squared, 0.975, na.rm = TRUE),
      
      sci_slope = mean(sci_slope, na.rm = TRUE),
      sci_slope_ci_low = quantile(sci_slope, 0.025, na.rm = TRUE),
      sci_slope_ci_high = quantile(sci_slope, 0.975, na.rm = TRUE),
      sci_r_squared = mean(sci_r_squared, na.rm = TRUE),
      sci_r_squared_ci_low = quantile(sci_r_squared, 0.025, na.rm = TRUE),
      sci_r_squared_ci_high = quantile(sci_r_squared, 0.975, na.rm = TRUE),
      
      mci_neg_slope = mean(mci_neg_slope, na.rm = TRUE),
      mci_neg_slope_ci_low = quantile(mci_neg_slope, 0.025, na.rm = TRUE),
      mci_neg_slope_ci_high = quantile(mci_neg_slope, 0.975, na.rm = TRUE),
      mci_neg_r_squared = mean(mci_neg_r_squared, na.rm = TRUE),
      mci_neg_r_squared_ci_low = quantile(mci_neg_r_squared, 0.025, na.rm = TRUE),
      mci_neg_r_squared_ci_high = quantile(mci_neg_r_squared, 0.975, na.rm = TRUE),
      
      mci_pos_slope = mean(mci_pos_slope, na.rm = TRUE),
      mci_pos_slope_ci_low = quantile(mci_pos_slope, 0.025, na.rm = TRUE),
      mci_pos_slope_ci_high = quantile(mci_pos_slope, 0.975, na.rm = TRUE),
      mci_pos_r_squared = mean(mci_pos_r_squared, na.rm = TRUE),
      mci_pos_r_squared_ci_low = quantile(mci_pos_r_squared, 0.025, na.rm = TRUE),
      mci_pos_r_squared_ci_high = quantile(mci_pos_r_squared, 0.975, na.rm = TRUE)
    )
  
  # Generate predicted values for interaction effects
  pred <- ggpredict(lm(as.formula(paste(protein, "~", ventricle, "*Subtype")), data = d_df), 
                    terms = c(ventricle, "Subtype"))
  
  # Find the minimum x value for positioning
  min_x <- min(d_df[[ventricle]], na.rm = TRUE)
  
  # Create the plot
  p <- ggplot() +
    # Add actual data points
    geom_point(data = d_df, aes(x = .data[[ventricle]], y = .data[[protein]], color = Subtype), 
               alpha = 0.5, size = 1) +
    # Add model predictions
    geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1) + 
    geom_line(data = pred, aes(x = x, y = predicted, color = group), size = 0.8) +
    
    # Set colors for all four subtypes
    scale_color_manual(values = c("AD" = "blue", "SCI" = "green", "MCI_A-T-" = "red", "MCI_A+T+" = "orange")) +
    scale_fill_manual(values = c("AD" = "blue", "SCI" = "green", "MCI_A-T-" = "red", "MCI_A+T+" = "orange")) +
    
    labs(
      title = paste("Interaction Effect of Subtype &", gsub("_", " ", ventricle), "Ventricle Volume on", protein, "Levels"),
      subtitle = sprintf("Full model R² = %.3f", boot_summary$full_r_squared),
      x = paste(gsub("_", " ", ventricle), "Ventricle Volume (cm³)"),
      y = paste("Predicted", protein, "Levels"),
      color = "Subtype",
      fill = "Subtype"
    ) +
    
    # Add annotations for all four subtypes
    annotate("text", x = min_x, y = max(d_df[[protein]], na.rm = TRUE) * 0.98, hjust = 0, vjust = 0,
             label = sprintf("AD: Slope = %.3f, R² = %.3f", 
                             boot_summary$ad_slope, boot_summary$ad_r_squared),
             size = 3, color = "blue") +
    
    annotate("text", x = min_x, y = max(d_df[[protein]], na.rm = TRUE) * 0.95, hjust = 0, vjust = 0,
             label = sprintf("SCI: Slope = %.3f, R² = %.3f", 
                             boot_summary$sci_slope, boot_summary$sci_r_squared),
             size = 3, color = "green") +
    
    annotate("text", x = min_x, y = max(d_df[[protein]], na.rm = TRUE) * 0.92, hjust = 0, vjust = 0,
             label = sprintf("MCI_A-T-: Slope = %.3f, R² = %.3f", 
                             boot_summary$mci_neg_slope, boot_summary$mci_neg_r_squared),
             size = 3, color = "red") +
    
    annotate("text", x = min_x, y = max(d_df[[protein]], na.rm = TRUE) * 0.89, hjust = 0, vjust = 0,
             label = sprintf("MCI_A+T+: Slope = %.3f, R² = %.3f", 
                             boot_summary$mci_pos_slope, boot_summary$mci_pos_r_squared),
             size = 3, color = "orange") +
    
    theme_minimal() +
    theme(
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9), 
      plot.title = element_text(size = 10, face = "bold", color = "black"),
      plot.subtitle = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  # Combine results and plot
  list(
    plot = p, 
    results = data.frame(
      Protein = protein,
      Ventricle = ventricle,
      Full_Model_R_Squared = boot_summary$full_r_squared,
      Full_Model_R_Squared_CI_Low = boot_summary$full_r_squared_ci_low,
      Full_Model_R_Squared_CI_High = boot_summary$full_r_squared_ci_high,
      
      AD_Slope = boot_summary$ad_slope,
      AD_Slope_CI_Low = boot_summary$ad_slope_ci_low,
      AD_Slope_CI_High = boot_summary$ad_slope_ci_high,
      AD_R_Squared = boot_summary$ad_r_squared,
      AD_R_Squared_CI_Low = boot_summary$ad_r_squared_ci_low,
      AD_R_Squared_CI_High = boot_summary$ad_r_squared_ci_high,
      
      SCI_Slope = boot_summary$sci_slope,
      SCI_Slope_CI_Low = boot_summary$sci_slope_ci_low,
      SCI_Slope_CI_High = boot_summary$sci_slope_ci_high,
      SCI_R_Squared = boot_summary$sci_r_squared,
      SCI_R_Squared_CI_Low = boot_summary$sci_r_squared_ci_low,
      SCI_R_Squared_CI_High = boot_summary$sci_r_squared_ci_high,
      
      MCI_Neg_Slope = boot_summary$mci_neg_slope,
      MCI_Neg_Slope_CI_Low = boot_summary$mci_neg_slope_ci_low,
      MCI_Neg_Slope_CI_High = boot_summary$mci_neg_slope_ci_high,
      MCI_Neg_R_Squared = boot_summary$mci_neg_r_squared,
      MCI_Neg_R_Squared_CI_Low = boot_summary$mci_neg_r_squared_ci_low,
      MCI_Neg_R_Squared_CI_High = boot_summary$mci_neg_r_squared_ci_high,
      
      MCI_Pos_Slope = boot_summary$mci_pos_slope,
      MCI_Pos_Slope_CI_Low = boot_summary$mci_pos_slope_ci_low,
      MCI_Pos_Slope_CI_High = boot_summary$mci_pos_slope_ci_high,
      MCI_Pos_R_Squared = boot_summary$mci_pos_r_squared,
      MCI_Pos_R_Squared_CI_Low = boot_summary$mci_pos_r_squared_ci_low,
      MCI_Pos_R_Squared_CI_High = boot_summary$mci_pos_r_squared_ci_high
    )
  )
}

# Create PDF for plots
pdf("finalized_interaction_plots_all_subtypes_040425.pdf", width = 11, height = 8.5)

# Prepare for storing results
all_plots <- list()
all_results <- data.frame()

# Set protein and ventricle columns
protein_cols <- colnames(d_df[7:55])
ventricle_vars <- c("Left_Lateral", "Right_Lateral", "Third", "Fourth")

# Generate plots and results for each protein and ventricle
for (protein in protein_cols) {
  protein_plots <- list()
  for (ventricle in ventricle_vars) {
    # Create plot and get results
    plot_result <- create_bootstrapped_interaction_plot(protein, ventricle)
    
    # Store plot
    protein_plots[[ventricle]] <- plot_result$plot
    
    # Accumulate results
    all_results <- rbind(all_results, plot_result$results)
  }
  
  # Arrange plots for the current protein
  grid.arrange(grobs = protein_plots, ncol = 2)
}

# Close PDF device
dev.off()

# Print results
print(all_results)

all_results <- all_results %>% arrange(desc(Full_Model_R_Squared))

# Load required libraries
library(ggplot2)
library(gridExtra) # For arranging multiple plots
library(grid)      # For textGrob

# Statistical analysis to compute p-values
proteins_of_interest <- c("GAP43", "DDAH1", "AMPH", "AQP4")

# Modified to include all four subtypes
stat_summary <- data.frame()

for(protein in proteins_of_interest) {
  # Run ANOVA across all subtypes
  anova_result <- summary(aov(d_df[[protein]] ~ d_df$Subtype))
  p_value_anova <- anova_result[[1]]$`Pr(>F)`[1]
  
  # Generate pairwise stats
  subtype_stats <- data.frame(
    Protein = protein,
    AD_Mean = mean(d_df[d_df$Subtype == "AD", protein], na.rm = TRUE),
    AD_SD = sd(d_df[d_df$Subtype == "AD", protein], na.rm = TRUE),
    SCI_Mean = mean(d_df[d_df$Subtype == "SCI", protein], na.rm = TRUE),
    SCI_SD = sd(d_df[d_df$Subtype == "SCI", protein], na.rm = TRUE),
    MCI_Aneg_Tneg_Mean = mean(d_df[d_df$Subtype == "MCI_A-T-", protein], na.rm = TRUE),
    MCI_Aneg_Tneg_SD = sd(d_df[d_df$Subtype == "MCI_A-T-", protein], na.rm = TRUE),
    MCI_Apos_Tpos_Mean = mean(d_df[d_df$Subtype == "MCI_A+T+", protein], na.rm = TRUE),
    MCI_Apos_Tpos_SD = sd(d_df[d_df$Subtype == "MCI_A+T+", protein], na.rm = TRUE),
    ANOVA_P_Value = p_value_anova
  )
  
  # Add pairwise t-test p-values
  subtypes <- c("AD", "SCI", "MCI_A-T-", "MCI_A+T+")
  for(i in 1:(length(subtypes)-1)) {
    for(j in (i+1):length(subtypes)) {
      group1 <- subtypes[i]
      group2 <- subtypes[j]
      t_test_result <- t.test(
        d_df[d_df$Subtype == group1, protein],
        d_df[d_df$Subtype == group2, protein]
      )
      col_name <- paste0("P_", group1, "_vs_", group2)
      subtype_stats[[col_name]] <- t_test_result$p.value
    }
  }
  
  # Calculate fold changes between all group combinations
  for(i in 1:(length(subtypes)-1)) {
    for(j in (i+1):length(subtypes)) {
      group1 <- subtypes[i]
      group2 <- subtypes[j]
      mean1 <- mean(d_df[d_df$Subtype == group1, protein], na.rm = TRUE)
      mean2 <- mean(d_df[d_df$Subtype == group2, protein], na.rm = TRUE)
      fc_col_name <- paste0("FC_", group1, "_vs_", group2)
      subtype_stats[[fc_col_name]] <- mean1 / mean2
    }
  }
  
  stat_summary <- rbind(stat_summary, subtype_stats)
}

# Adjust p-values for multiple comparisons
p_value_columns <- grep("^P_", names(stat_summary))
all_p_values <- unlist(stat_summary[, p_value_columns])
all_p_adjusted <- p.adjust(all_p_values, method = "fdr")

# Put adjusted p-values back
p_adj_index <- 1
for(i in 1:nrow(stat_summary)) {
  for(j in p_value_columns) {
    col_name <- names(stat_summary)[j]
    adj_col_name <- paste0("FDR_", sub("^P_", "", col_name))
    stat_summary[i, adj_col_name] <- all_p_adjusted[p_adj_index]
    p_adj_index <- p_adj_index + 1
  }
}

# Save detailed statistics
write.csv(stat_summary, "protein_comparison_stats_all_subtypes.csv", row.names = FALSE)

# Updated function to generate boxplots with all subtypes in specified order
create_protein_boxplot <- function(data, protein_name, p_value_anova) {
  # Convert Subtype to factor with specified order
  data$Subtype <- factor(data$Subtype, levels = c("AD", "MCI_A+T+", "MCI_A-T-", "SCI"))
  
  ggplot(data, aes(x = Subtype, y = .data[[protein_name]], fill = Subtype)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2, aes(color = Subtype)) +
    scale_fill_manual(values = c("AD" = "blue", "MCI_A+T+" = "orange", 
                                 "MCI_A-T-" = "red", "SCI" = "green")) +
    scale_color_manual(values = c("AD" = "blue", "MCI_A+T+" = "orange", 
                                  "MCI_A-T-" = "red", "SCI" = "green")) +
    labs(
      title = paste("Distribution of", protein_name, "Levels by Subtype"),
      subtitle = paste("P-value:", signif(p_value_anova, digits = 3)),
      y = paste(protein_name, "Expression Level"),
      x = "Subtype"
    ) +
    theme_minimal() +
    theme(
      legend.position = "NULL",
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x labels for better readability
    ) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black")
}

# Before running your lapply function, convert Subtype to factor in the original dataset
d_df$Subtype <- factor(d_df$Subtype, levels = c("AD", "MCI_A+T+", "MCI_A-T-", "SCI"))

# Generate plots for each protein using the computed ANOVA p-values
protein_plots <- lapply(1:nrow(stat_summary), function(i) {
  create_protein_boxplot(d_df, stat_summary$Protein[i], stat_summary$ANOVA_P_Value[i])
})

# Arrange and save plots with updated title
combined_plot <- grid.arrange(
  grobs = protein_plots,
  ncol = 2,
  top = textGrob(
    "Protein Level Comparison Across AD, MCI_A+T+, MCI_A-T-, and SCI Groups",
    gp = gpar(fontsize = 18, fontface = "bold")
  )
)

# Save the plot
ggsave("protein_comparison_boxplots_all_subtypes.png", combined_plot, width = 14.5, height = 11)

read_delim("protein_comparison_stats_all_subtypes.csv", delim = ";") %>% 
  
  # Continue by creating a table visualizing all the pairwise t-test results. The table boxes should be colored based on the significance of the p-values. 
