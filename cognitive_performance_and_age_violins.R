# Load Required Libraries
library(tidyverse)    # Data wrangling and visualization
library(ggsignif)     # Add significance annotations to ggplot
library(ggpubr)       # Statistical tests and plotting helpers
library(patchwork)    # Combine multiple ggplots

# Load and Preprocess Dataset
data <- read_delim("heatmap_data_norm_cogn.csv", delim = ";") %>%
  rename(Subgroup = "Subtype") %>%
  mutate(Subgroup = recode(Subgroup, "SCI" = "SCD"))

# Validate Required Columns
cognitive_columns <- c("MMSE", "MoCA", "KOD", "RAVLT", "RCF", "Age")
missing_columns <- setdiff(cognitive_columns, colnames(data))

if (length(missing_columns) > 0) {
  stop(paste("Missing columns in dataset:", paste(missing_columns, collapse = ", ")))
}

# Set Factor Levels and Define Colors
subgroup_colors <- c(
  "AD" = "blue",
  "MCI_A+T+" = "orange",
  "MCI_A-T-" = "red",
  "SCD" = "green"
)

data$Subgroup <- factor(data$Subgroup, levels = names(subgroup_colors))

# Define Function to Create Annotated Violin Plots
make_violin_plot <- function(data, score_col, title_prefix) {
  # Drop rows with missing values in the target column
  df <- data %>% drop_na(all_of(score_col))
  
  # Define plot title and y-axis bounds
  plot_title <- paste(title_prefix, "Performance")
  y_min <- min(df[[score_col]], na.rm = TRUE) * 0.95
  y_max <- max(df[[score_col]], na.rm = TRUE) * 1.20
  
  # Define subgroup comparisons (all vs. SCD)
  comparison_pairs <- list(
    c("AD", "SCD"),
    c("MCI_A+T+", "SCD"),
    c("MCI_A-T-", "SCD")
  )
  
  # Perform Wilcoxon tests
  p_vals <- sapply(comparison_pairs, function(groups) {
    if (all(groups %in% df$Subgroup)) {
      tryCatch(
        wilcox.test(
          df[[score_col]][df$Subgroup == groups[1]],
          df[[score_col]][df$Subgroup == groups[2]]
        )$p.value,
        error = function(e) NA_real_
      )
    } else {
      NA_real_
    }
  })
  
  # Convert p-values to significance symbols
  symbols <- case_when(
    p_vals <= 0.0001 ~ "****",
    p_vals <= 0.001  ~ "***",
    p_vals <= 0.01   ~ "**",
    p_vals <= 0.05   ~ "*",
    TRUE             ~ "ns"
  )
  
  # Set y-axis annotation positions
  y_positions <- seq(y_max * 0.95, y_max * 0.80, length.out = 3)
  
  # Generate base plot
  p <- ggplot(df, aes(x = Subgroup, y = .data[[score_col]], fill = Subgroup)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black", width = 0.8) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7, color = "black") +
    geom_jitter(width = 0.1, size = 2, alpha = 0.4) +
    scale_fill_manual(values = subgroup_colors, name = "Subgroup") +
    labs(
      title = plot_title,
      y = paste(title_prefix, "Score"),
      x = NULL
    ) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(size = 16, face = "bold"),
      axis.title.x = element_blank(),
      legend.position = "bottom"
    )
  
  # Add significance annotations
  for (i in seq_along(comparison_pairs)) {
    if (!is.na(p_vals[i])) {
      p <- p + geom_signif(
        comparisons = list(comparison_pairs[[i]]),
        annotations = symbols[i],
        y_position = y_positions[i],
        tip_length = 0.01,
        textsize = 6
      )
    }
  }
  
  return(p)
}

# Generate Plots for All Cognitive Measures and Age
plot_list <- list(
  make_violin_plot(data, "MMSE", "MMSE"),
  make_violin_plot(data, "MoCA", "MoCA"),
  make_violin_plot(data, "KOD", "KOD"),
  make_violin_plot(data, "RAVLT", "RAVLT"),
  make_violin_plot(data, "RCF", "RCF"),
  make_violin_plot(data, "Age", "Age")
)

# Combine All Plots into a Final Figure
final_figure <- wrap_plots(plot_list, ncol = 2, guides = "collect", tag_level = "new") +
  plot_annotation(
    title = "Cognitive Performance and Age Distribution Across Clinical Subgroups",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 26, face = "bold")
    ),
    tag_levels = list(c("A", "B", "C", "D", "E", "F")),
    tag_prefix = "",
    tag_suffix = ""
  ) &
  theme(legend.position = "bottom")

# Save the Combined Plot as an Image
ggsave(
  filename = "cognitive_age_violins.png",
  plot = final_figure,
  width = 16,
  height = 18,
  dpi = 300
)
