# --- Setup and Libraries ---
setwd("C:/Users/User/Desktop/RMT")

library(tidyverse)
library(limma)
library(ggeffects)
library(ggplot2)
library(rsample)
library(broom)
library(vsn)
library(patchwork)

# --- Load and Format Data ---
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
  filter(Subgroup %in% c("AD", "SCD")) %>%
  as.data.frame()

# --- Factor Conversion ---
d_df$Subgroup <- factor(d_df$Subgroup, levels = c("AD", "SCD"))
d_df$Sex <- as.factor(d_df$Sex)

# --- Define Variables ---
ventricle_vars <- c("Left_Lateral", "Right_Lateral", "Third", "Fourth", "Total_Ventricle")
meta_cols <- c("Subgroup", "Sex", "Age", "A", "T", "ab42_40", "Qalb", "ICV", ventricle_vars, paste0(ventricle_vars, "_norm"), "Total_Ventricle")

# --- Normalize Ventricle Volumes ---
global_means <- d_df %>% summarize(across(all_of(ventricle_vars), ~ mean(., na.rm = TRUE)))
d_df <- d_df %>% mutate(across(all_of(ventricle_vars), ~ ((. - global_means[[cur_column()]]) / global_means[[cur_column()]]), .names = "{.col}_norm"))

# --- Identify and Normalize Protein Data ---
protein_cols <- setdiff(colnames(d_df), meta_cols)
protein_cols <- intersect(protein_cols, c("GAP43", "DDAH1", "AQP4", "AMPH", "TREM2", "NPTX1"))
exp_data <- d_df[, protein_cols, drop = FALSE]
normalized_exp_data <- normalizeBetweenArrays(as.matrix(exp_data), method = "quantile")
vsn_fit <- vsn2(as.matrix(normalized_exp_data))
normData <- as.data.frame(predict(vsn_fit, newdata = as.matrix(normalized_exp_data)))
colnames(normData) <- protein_cols
d_df <- bind_cols(d_df %>% select(all_of(meta_cols)), normData)

# --- Bootstrap Resampling Setup ---
set.seed(123)
bootstrap_samples <- bootstraps(d_df, times = 500)

# --- Function: Interaction Plot ---
create_bootstrapped_interaction_plot <- function(protein, ventricle) {
  ventricle_z <- paste0(ventricle, "_norm")
  ventricle_display_name <- ifelse(ventricle == "Total_Ventricle", "Total", ventricle)
  
  bootstrap_model <- function(split) {
    boot_sample <- analysis(split)
    full_formula <- as.formula(paste(protein, "~", ventricle_z, "*Subgroup + Sex + Age + ab42_40"))
    full_model <- tryCatch({ lm(full_formula, data = boot_sample) }, error = function(e) NULL)
    ad_data <- boot_sample[boot_sample$Subgroup == "AD", ]
    scd_data <- boot_sample[boot_sample$Subgroup == "SCD", ]
    ad_model <- tryCatch({ if(nrow(ad_data) > 1) lm(as.formula(paste(protein, "~", ventricle_z)), data = ad_data) else NULL }, error = function(e) NULL)
    scd_model <- tryCatch({ if(nrow(scd_data) > 1) lm(as.formula(paste(protein, "~", ventricle_z)), data = scd_data) else NULL }, error = function(e) NULL)
    
    tibble(
      Full_R_Squared = if(!is.null(full_model)) summary(full_model)$r.squared else NA,
      AD_Slope = if(!is.null(ad_model)) coef(ad_model)[2] else NA,
      AD_R_Squared = if(!is.null(ad_model)) summary(ad_model)$r.squared else NA,
      SCD_Slope = if(!is.null(scd_model)) coef(scd_model)[2] else NA,
      SCD_R_Squared = if(!is.null(scd_model)) summary(scd_model)$r.squared else NA
    )
  }
  
  bootstrap_results <- bootstrap_samples %>% mutate(model_stats = map(splits, bootstrap_model)) %>% unnest(model_stats)
  
  boot_summary <- bootstrap_results %>% summarise(
    Full_R_Squared = mean(Full_R_Squared, na.rm = TRUE),
    AD_Slope = mean(AD_Slope, na.rm = TRUE),
    AD_R_Squared = mean(AD_R_Squared, na.rm = TRUE),
    SCD_Slope = mean(SCD_Slope, na.rm = TRUE),
    SCD_R_Squared = mean(SCD_R_Squared, na.rm = TRUE)
  )
  
  y_range <- range(d_df[[protein]], na.rm = TRUE)
  y_span <- diff(y_range)
  annotation_y_ad <- y_range[2] - 0.05 * y_span
  annotation_y_scd <- y_range[2] - 0.12 * y_span
  min_x <- min(d_df[[ventricle_z]], na.rm = TRUE)
  
  full_model <- tryCatch({ lm(as.formula(paste(protein, "~", ventricle_z, "*Subgroup")), data = d_df) }, error = function(e) NULL)
  if (is.null(full_model)) return(NULL)
  
  pred <- tryCatch({ ggpredict(full_model, terms = c(ventricle_z, "Subgroup")) }, error = function(e) NULL)
  if (is.null(pred)) return(NULL)
  
  p <- ggplot() +
    geom_point(data = d_df, aes(x = .data[[ventricle_z]], y = .data[[protein]], color = Subgroup), alpha = 0.5, size = 1) +
    geom_ribbon(data = pred, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1) + 
    geom_line(data = pred, aes(x = x, y = predicted, color = group), size = 0.8) +
    scale_color_manual(values = c("AD" = "blue", "SCD" = "green")) +
    scale_fill_manual(values = c("AD" = "blue", "SCD" = "green")) +
    labs(
      title = paste("Subgroup ×", gsub("_", " ", ventricle_display_name), "Volume Effect on", protein),
      subtitle = sprintf("Full model R² = %.3f", boot_summary$Full_R_Squared),
      x = "LDN-Normalized Ventricle Volume [AU]",
      y = paste("Log2(",protein," Normalized Intensity) [AU]"),
      color = "Subgroup",
      fill = "Subgroup"
    ) +
    annotate("text", x = min_x, y = annotation_y_ad, hjust = 0, vjust = 0, 
             label = sprintf("AD: Slope = %.3f, R² = %.3f", boot_summary$AD_Slope, boot_summary$AD_R_Squared),
             size = 4, color = "blue", fontface = "bold") +
    annotate("text", x = min_x, y = annotation_y_scd, hjust = 0, vjust = 0,
             label = sprintf("SCD: Slope = %.3f, R² = %.3f", boot_summary$SCD_Slope, boot_summary$SCD_R_Squared),
             size = 4, color = "green", fontface = "bold") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text = element_text(size = 14),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  
  list(plot = p)
}

# --- Function: Generate PNG for AQP4 ---
generate_aqp4_png <- function() {
  message("Creating separate PNG for AQP4")
  png_filename <- "interaction_plots_AQP4.png"
  if (dev.cur() > 1) dev.off()
  
  tryCatch({
    png(png_filename, width = 14, height = 16, units = "in", res = 300)
    protein <- "AQP4"
    protein_plots <- list()
    
    for (ventricle in ventricle_vars) {
      plot_result <- create_bootstrapped_interaction_plot(protein, ventricle)
      if (!is.null(plot_result)) {
        protein_plots[[ventricle]] <- plot_result$plot
      } else {
        protein_plots[[ventricle]] <- ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = paste("Error for", protein, "and", ventricle)) +
          theme_void() + xlim(0, 1) + ylim(0, 1)
      }
    }
    
    first_row <- wrap_plots(plotlist = protein_plots[c("Left_Lateral", "Right_Lateral")], ncol = 2)
    second_row <- wrap_plots(plotlist = protein_plots[c("Third", "Fourth")], ncol = 2)
    third_row <- wrap_plots(plotlist = list(protein_plots[["Total_Ventricle"]], ggplot() + theme_void()), ncol = 2)
    
    combined_plot <- first_row / second_row / third_row +
      plot_layout(guides = "collect") +
      plot_annotation(tag_levels = "A") &
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "center"
      )
    
    final_page <- combined_plot +
      plot_annotation(
        title = "Subgroup × Ventricular Volume Effects on AQP4",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
    
    print(final_page)
    dev.off()
    message("Successfully created PNG for AQP4")
  }, error = function(e) {
    message("Error creating PNG for AQP4: ", e)
    if (dev.cur() > 1) dev.off()
  })
}

# --- Run Plot for AQP4 ---
generate_aqp4_png()
