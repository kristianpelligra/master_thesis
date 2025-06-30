# --- Setup and Libraries ---
setwd("C:/Users/User/Desktop/RMT")

library(tidyverse)
library(limma)
library(vsn)
library(PCAtools)
library(grid)
library(cowplot)
library(ggplotify)
library(gridGraphics)

# --- Load and Prepare Data ---
d_df <- read_delim("heatmap_data.csv", delim = ";") %>%
  rename(
    Subgroup = subtype,
    Sex = gender,
    Age = sampling_age,
    "NfL" = nfl,
    A = a, T = t, N = n,
    `LLV Volume` = Left_Lateral_Ventricle,
    `RLV Volume` = Right_Lateral_Ventricle,
    "TV Volume" = Third_Ventricle,
    "FV Volume" = Fourth_Ventricle,
    `Albumin Quotient` = alb_ratio,
    `ToV Volume` = Total_Adjusted_CSF_Volume,
    `LH Mean Thickness` = lh_mean_thickness,
    `RH Mean Thickness` = rh_mean_thickness,
    `ab42/40` = ab42_40,
    `p-tau` = ptau, `t-tau` = ttau
  ) %>%
  mutate(Subgroup = ifelse(Subgroup == "SCI", "SCD", Subgroup)) %>%
  as.data.frame()

rownm <- make.names(d_df[, 1], unique = TRUE)

# --- Metadata Setup ---
metadata <- d_df[, c(2:6, 56:70)] %>%
  select(
    Age, Sex, Subgroup, ICV, `LLV Volume`, `RLV Volume`, `TV Volume`, `FV Volume`,
    `LH Mean Thickness`, `RH Mean Thickness`, A, T, N, `ab42/40`, `t-tau`, `p-tau`, `NfL`,
    `Albumin Quotient`, `ToV Volume`, `Total CSF Protein`
  ) %>%
  as.data.frame()

metadata$Age <- as.numeric(metadata$Age)
metadata$Sex <- factor(metadata$Sex, levels = c("F", "M"))
metadata$Subgroup <- factor(metadata$Subgroup, levels = c("AD", "MCI_A+T+", "MCI_A-T-", "SCD"))
rownames(metadata) <- rownm

# --- Normalize Data ---
columns_to_divide <- c("LLV Volume", "ICV", "RLV Volume", "TV Volume", "FV Volume", "ToV Volume")
d_df[columns_to_divide] <- d_df[columns_to_divide] / 1000

exp_data <- d_df[, 7:55]
normalized_exp_data <- normalizeBetweenArrays(as.matrix(exp_data), method = "quantile")
vsn_fit <- vsn2(as.matrix(normalized_exp_data))
normData <- as.data.frame(predict(vsn_fit, newdata = as.matrix(normalized_exp_data)))
rownames(normData) <- rownm
normData_t <- t(normData)

# --- PCA ---
p <- pca(normData_t, metadata = metadata)

# --- Wilcoxon tests for PC1 and PC3 by Subgroup ---
pc_scores <- as.data.frame(p$rotated)
pc_scores$Subgroup <- metadata$Subgroup

ref_group <- "SCD"
subgroups <- levels(pc_scores$Subgroup)
test_groups <- setdiff(subgroups, ref_group)

pc_wilcox_results <- data.frame()

for (pc in c("PC1", "PC3")) {
  for (grp in test_groups) {
    values_grp <- pc_scores %>% filter(Subgroup == grp) %>% pull(pc)
    values_ref <- pc_scores %>% filter(Subgroup == ref_group) %>% pull(pc)
    test_result <- wilcox.test(values_grp, values_ref, alternative = "two.sided")
    pc_wilcox_results <- rbind(pc_wilcox_results, data.frame(
      PC = pc,
      Comparison = paste(grp, "vs", ref_group),
      P_Value = test_result$p.value
    ))
  }
}
pc_wilcox_results$FDR <- p.adjust(pc_wilcox_results$P_Value, method = "fdr")
print(pc_wilcox_results)

# --- Generate All Plots ---
pbiplot <- biplot(
  p, x = "PC1", y = "PC3",
  showLoadings = TRUE,
  lengthLoadingsArrowsFactor = 1.5,
  sizeLoadingsNames = 5.6,
  lab = NULL,
  colby = "Subgroup",
  colkey = c("AD" = "blue", "MCI_A+T+" = "orange", "MCI_A-T-" = "red", "SCD" = "green"),
  shape = "Sex", shapekey = c("M" = 15, "F" = 17),
  drawConnectors = FALSE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  legendLabSize = 12,
  legendIconSize = 12,
  legendPosition = "right",
  pointSize = 4.2,
  ellipse = TRUE,
  ellipseType = "t",
  ellipseLevel = 0.95,
  ellipseFill = TRUE,
  ellipseAlpha = 1/4,
  ellipseLineSize = 0,
  returnPlot = TRUE
)

vars_to_use <- c("Age", "Sex", "ICV", "LLV Volume", "RLV Volume", "TV Volume", "FV Volume",
                 "LH Mean Thickness", "RH Mean Thickness", "A", "T", "N", "ab42/40", 
                 "t-tau", "p-tau", "NfL", "Albumin Quotient", "ToV Volume", "Total CSF Protein")

peigencor_r2 <- eigencorplot(
  p, components = getComponents(p, 1:4), metavars = vars_to_use,
  main = NULL, plotRsquared = TRUE,
  col = c("white", "cornsilk1", "gold", "forestgreen", "darkgreen"),
  cexCorval = 1.2, posColKey = "right", cexLabColKey = 1.8,
  fontCorval = 2, posLab = "bottomleft", rotLabX = 45,
  scale = TRUE, corFUN = "pearson", corUSE = "pairwise.complete.obs",
  signifSymbols = c("****", "***", "**", "*", ""),
  signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), returnPlot = TRUE
)

p_clinical_corr <- eigencorplot(
  p, components = getComponents(p, 1:4), metavars = vars_to_use,
  main = NULL, plotRsquared = FALSE,
  cexMain = 1.5, cexCorval = 1.0, fontCorval = 2, posLab = "bottomleft",
  rotLabX = 45, scale = TRUE, corFUN = "pearson", corUSE = "pairwise.complete.obs",
  signifSymbols = c("****", "***", "**", "*", ""),
  signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), returnPlot = TRUE
)

# --- Combine Plots in Grid ---
final_plot_grid <- plot_grid(
  pbiplot, peigencor_r2, p_clinical_corr,
  labels = c("A: PCA Biplot", "B: PC–Variable Correlations (r²)", "C: PC–Variable Correlations (r)"),
  label_size = 18, label_fontface = "bold",
  ncol = 2, align = "hv"
)

title <- ggdraw() + draw_label(
  "Principal Component Analysis (PCA)",
  fontface = 'bold', size = 24, hjust = 0.5
)

final_combined <- plot_grid(
  title, final_plot_grid,
  ncol = 1, rel_heights = c(0.1, 1)
)

# --- Save Final Image ---
ggsave("combined_pca_analysis_plots.png", final_combined,
       width = 16, height = 10, dpi = 300, bg = "white")
