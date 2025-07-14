# Load required libraries
library(tidyverse)
library(limma)
library(vsn)
library(ggpubr)
library(patchwork)

# Load and format data
d_df <- read_delim("heatmap_data.csv", delim = ";") %>%
  rename(
    Subgroup = subtype, Sex = gender, Age = sampling_age, 
    A = a, T = t, 
    Left_Lateral = Left_Lateral_Ventricle, 
    Right_Lateral = Right_Lateral_Ventricle, 
    Third = Third_Ventricle, Fourth = Fourth_Ventricle, 
    Total_Ventricle = Total_Adjusted_CSF_Volume, Qalb = alb_ratio
  ) %>%
  mutate(across(c(Left_Lateral, Right_Lateral, Third, Fourth, Total_Ventricle, ICV), ~ . / 1000)) %>%
  mutate(Subgroup = ifelse(Subgroup == "SCI", "SCD", Subgroup)) %>%
  filter(Subgroup %in% c("AD", "MCI_A+T+", "MCI_A-T-", "SCD")) %>% 
  as.data.frame()

# Protein normalization
protein_cols <- c("AQP4", "AMPH", "DDAH1", "SNCB", "GAP43", "ARPP21")
exp_data <- d_df[, protein_cols]

normalized_exp_data <- normalizeBetweenArrays(as.matrix(exp_data), method = "quantile")
vsn_fit <- vsn2(as.matrix(normalized_exp_data))
normData <- as.data.frame(predict(vsn_fit, newdata = as.matrix(normalized_exp_data)))

# Combine normalized proteins with metadata
metadata_cols <- setdiff(names(d_df), protein_cols)
d_df <- cbind(d_df[, metadata_cols], normData)

# Factor conversion
d_df$Subgroup <- factor(d_df$Subgroup, levels = c("AD", "MCI_A+T+", "MCI_A-T-", "SCD"))

# Define colors and comparisons
subgroup_colors <- c("AD" = "blue", "MCI_A+T+" = "orange", "MCI_A-T-" = "red", "SCD" = "green")
my_comparisons <- list(c("SCD", "AD"), c("SCD", "MCI_A+T+"), c("SCD", "MCI_A-T-"))

# Create violin plots function
create_protein_plot <- function(data, protein) {
  ggplot(data, aes_string(x = "Subgroup", y = protein, fill = "Subgroup")) +
    geom_violin(alpha = 0.7, trim = FALSE, color = "black", width = 0.8) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7, color = "black") +
    geom_jitter(width = 0.1, size = 2, alpha = 0.4) +
    scale_fill_manual(values = subgroup_colors, name = "Subgroup") +
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.01,
      size = 6
    ) +
    labs(
      title = paste(protein, "Across Subgroups"),
      x = NULL,
      y = paste("Log₂(", protein, "Normalized Intensity) [AU]")
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16),
      legend.position = "bottom",
      legend.key = element_rect(fill = "white", color = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 14)
    )
}

# Create all plots
protein_plots <- map(protein_cols, ~create_protein_plot(d_df, .x))
names(protein_plots) <- protein_cols

# Combine plots
final_plot <- wrap_plots(protein_plots, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "CSF Protein Comparison Across Clinical Subgroups",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold")),
    tag_levels = list(c("A", "B", "C", "D", "E", "F"))
  ) &
  theme(legend.position = "bottom")

# Save the final plot
ggsave("csf_protein_level_comparison_across_subgroups.png", final_plot, width = 12, height = 15.25)

# Print the plot
print(final_plot)