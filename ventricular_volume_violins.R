# Ventricular volumes for each subgroup are compared here

library(tidyverse)    
library(ggpubr)       
library(ggplot2)     
library(patchwork)    

# Load and Preprocess the Data

data <- read_delim("heatmap_data_norm.csv", delim = ";") %>%
  # Rename columns for consistency or clarity
  rename(
    Subgroup = "Subtype",
    `Total Ventricle` = Total_Ventricle
  ) %>%
  # Recode and relevel the 'Subgroup' factor
  mutate(
    Subgroup = recode(Subgroup, "SCI" = "SCD"),  # Rename SCI to SCD
    Subgroup = factor(Subgroup, levels = c("AD", "MCI_A+T+", "MCI_A-T-", "SCD"))  # Set order
  )

# Define Custom Colors per Subgroup

subgroup_colors <- c(
  "AD" = "blue", 
  "MCI_A+T+" = "orange", 
  "MCI_A-T-" = "red", 
  "SCD" = "green"
)

# Define Statistical Comparisons

my_comparisons <- list(
  c("SCD", "AD"),
  c("SCD", "MCI_A-T-"),
  c("SCD", "MCI_A+T+")
)

# Function to Create Violin + Boxplot + Jitter Plot

create_ventricle_plot <- function(data, ventricle_column, title) {
  ggplot(data, aes_string(x = "Subgroup", y = ventricle_column, fill = "Subgroup")) +
    geom_violin(alpha = 0.7, trim = FALSE, color = "black", width = 0.8) +  # Violin plot
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7, color = "black") +  # Boxplot inside violin
    geom_jitter(width = 0.1, size = 2, alpha = 0.4) +  # Jittered data points
    scale_fill_manual(values = subgroup_colors, name = "Subgroup") +  # Custom colors
    stat_compare_means(
      comparisons = my_comparisons,   # Pairwise Wilcoxon tests
      method = "wilcox.test",
      label = "p.signif",             # Show significance stars
      tip.length = 0.01,
      size = 6
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Volume (cm³)"
    ) +
    theme_minimal(base_size = 16) +  # Clean theme
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

# Generate Plots for Each Ventricle

plots <- list(
  create_ventricle_plot(data, "Left_Lateral", "Left Lateral Ventricle"),
  create_ventricle_plot(data, "Right_Lateral", "Right Lateral Ventricle"),
  create_ventricle_plot(data, "Third", "Third Ventricle"),
  create_ventricle_plot(data, "Fourth", "Fourth Ventricle"),
  create_ventricle_plot(data, "Total Ventricle", "Total Ventricle")
)

# Combine All Plots into a Single Figure

final_plot <- wrap_plots(plots, ncol = 2, guides = "collect", tag_level = "new") +
  plot_annotation(
    title = "Ventricle Volume Comparison Across Clinical Subgroups",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 26, face = "bold")
    ),
    tag_levels = list(c("A", "B", "C", "D", "E")),  # Labels for subplots
    tag_prefix = "",
    tag_suffix = ""
  ) &
  theme(legend.position = "bottom")


# Save Final Plot to File

ggsave("ventricle_volume_per_subgroup.png", plot = final_plot, width = 16, height = 18, dpi = 300)


# Calculate and Display Mean Ventricle Volumes by Subgroup

ventricle_means <- data %>%
  group_by(Subgroup) %>%
  summarise(
    Mean_Left_Lateral = mean(Left_Lateral, na.rm = TRUE),
    Mean_Right_Lateral = mean(Right_Lateral, na.rm = TRUE),
    Mean_Third = mean(Third, na.rm = TRUE),
    Mean_Fourth = mean(Fourth, na.rm = TRUE),
    Mean_Total = mean(`Total Ventricle`, na.rm = TRUE)
  )

# Print the calculated mean values
print(ventricle_means)