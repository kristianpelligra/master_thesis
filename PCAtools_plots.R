# Set working directory
setwd("C:/Users/User/Desktop/RMT")

# Load required libraries
library(tidyverse)
library(limma)
library(vsn)
library(PCAtools)
library(ggplot2)
library(cli)
library(grid)

# Load and format data
d_df <- read_delim("heatmap_data.csv", delim = ";") %>%
  rename(
    Subtype = subtype, 
    Sex = gender, 
    Age = sampling_age, 
    "NfL" = nfl,
    A = a, 
    T = t, 
    N = n, 
    `LLV Volume` = Left_Lateral_Ventricle, 
    `RLV Volume` = Right_Lateral_Ventricle, 
    "TV Volume" = Third_Ventricle, 
    "FV Volume" = Fourth_Ventricle, 
    `Albumin Quotient` = alb_ratio, 
    `Total Ventricle` = Total_Adjusted_CSF_Volume, 
    `LH Mean Thickness` = lh_mean_thickness, 
    `RH Mean Thickness` = rh_mean_thickness, 
    `ab42/40` = ab42_40
  ) %>% 
  as.data.frame()

# Create unique row names
rownm <- make.names(d_df[, 1], unique = TRUE)

# Separate metadata and expression data
metadata <- d_df[, c(2:6, 56:69)] %>% 
  as.data.frame() %>%
  dplyr::select(
    "Age", "Sex", "Subtype", "ICV", "LLV Volume", 
    "RLV Volume", "TV Volume", "FV Volume", 
    "LH Mean Thickness", "RH Mean Thickness", 
    "A", "T", "N", "ab42/40", "ttau", "ptau", "NfL", "Albumin Quotient"
  )

# Convert Age to numeric and Sex to factor
metadata$Age <- as.numeric(metadata$Age)
metadata$Sex <- factor(metadata$Sex, levels = c('F', 'M'))

rownames(metadata) <- rownm

# Normalize and filter expression data
columns_to_divide <- c("LLV Volume", "ICV", "RLV Volume", "TV Volume", "FV Volume", "Total Ventricle")
d_df[columns_to_divide] <- d_df[columns_to_divide] / 1000

# Collecting and normalizing protein columns
exp_data <- d_df[, 7:55]
normalized_exp_data <- normalizeBetweenArrays(as.matrix(exp_data), method = "quantile")
vsn_fit <- vsn2(as.matrix(normalized_exp_data))
normData <- as.data.frame(predict(vsn_fit, newdata = as.matrix(normalized_exp_data)))
rownames(normData) <- rownm
normData_t <- t(normData)

# Perform PCA
p <- pca(normData_t, metadata = metadata)

horn <- parallelPCA(normData_t)
horn_ncomp <- horn$n

elbow <- findElbowPoint(p$variance)
elbow_ncomp <- elbow

pscree <- screeplot(
  p, components = getComponents(p, 1:5),
  axisLabSize = 10, 
  title = NULL,
  returnPlot = FALSE  # Set to TRUE to modify the plot
) +
  geom_vline(xintercept = horn_ncomp, color = "black", linetype = "dotted", linewidth = 1) + 
  geom_vline(xintercept = elbow_ncomp, color = "black", linetype = "dotted", linewidth = 1) +
  annotate("text", x = horn_ncomp + 0.3, y = max(p$variance) * 0.9, label = "Horn", color = "black", size = 4) +
  annotate("text", x = elbow_ncomp + 0.3, y = max(p$variance) * 0.9, label = "Elbow", color = "black", size = 4)

# Pairs plot with tighter layout
ppairs <- pairsplot(
  p, 
  components = getComponents(p, 1:3),
  triangle = TRUE, 
  trianglelabSize = 8,
  hline = 0, 
  vline = 0,
  pointSize = 0.5, 
  gridlines.major = FALSE, 
  gridlines.minor = FALSE,
  colby = 'Subtype',
  colkey = c('AD' = 'blue', 'MCI_A-T-' = 'red', 'MCI_A+T+' = 'orange', 'SCI' = 'green'), 
  legendLabSize = 10,
  legendIconSize = 10,
  title = '', 
  plotaxes = FALSE,
  margingaps = unit(c(0.005, 0.005, 0.005, 0.005), 'cm'),
  returnPlot = TRUE)

# PCA Bi-plot with ellipses
pbiplot <- biplot(
  p,
  title = "PCA Biplot",
  x = 'PC1', y = 'PC3', 
  subtitle = "PC1 vs PC3",
  showLoadings = TRUE,
  lengthLoadingsArrowsFactor = 2,
  sizeLoadingsNames = 5,
  lab = NULL,
  colby = 'Subtype',
  colkey = c('AD' = 'blue', 'MCI_A-T-' = 'red', 'MCI_A+T+' = 'orange', 'SCI' = 'green'), 
  encircle = FALSE,  # Set to FALSE, we'll add ellipses manually
  hline = 0,
  vline = 0,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  pointSize = 3,
  legendLabSize = 10,
  legendIconSize = 10,
  shape = 'Sex',
  shapekey = c('M' = 15, 'F' = 17),
  drawConnectors = FALSE,
  titleLabSize = 20,
  subtitleLabSize = 16,
  axisLabSize = 14,
  legendTitleSize = 12,
  legendPosition = 'right',
  returnPlot = TRUE) + 
  theme(
    plot.subtitle = element_text(face = "bold", size = 16)  # Bold and adjust subtitle size
  )

# Extract PC data for ellipses
pc_data <- as.data.frame(p$rotated)
pc_data$Subtype <- metadata$Subtype[match(rownames(pc_data), rownames(metadata))]

# Add 95% confidence ellipses
pbiplot <- pbiplot + 
  stat_ellipse(data = pc_data, 
               aes(x = PC1, y = PC3, fill = Subtype, group = Subtype), 
               geom = "polygon", 
               alpha = 0.25, 
               level = 0.95, 
               show.legend = FALSE) +
  scale_fill_manual(values = c('AD' = 'blue', 'MCI_A-T-' = 'red', 'MCI_A+T+' = 'orange', 'SCI' = 'green'))

# Display the plot
pbiplot

# Loading plot with adjusted aesthetics
ploadings <- plotloadings(
  p, 
  rangeRetain = 0.01, 
  labSize = 4,
  subtitle = "Top 1% Variables",
  subtitleLabSize = 10,
  axisLabSize = 10,
  shape = 21, 
  shapeSizeRange = c(2, 4),
  col = c('blue', 'white', 'red3'),
  legendPosition = 'none',
  drawConnectors = FALSE,
  returnPlot = TRUE
)

# First plot with Pearson r² Clinical Correlates
pearson_eigencor <- eigencorplot(
  p,
  main = "Principal Component Pearson r² Clinical Correlates",
  components = getComponents(p, 1:4),
  metavars = c("Age", "Sex", "Subtype", "LLV Volume", 
               "RLV Volume", "TV Volume", "FV Volume", 
               "LH Mean Thickness", "RH Mean Thickness", 
               "ab42/40", "ttau", "ptau", "NfL", "Albumin Quotient"),
  col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
  cexCorval = 1.0,
  posColKey = "right",
  cexLabColKey = 1.5,
  fontCorval = 2,
  posLab = 'bottomleft', 
  rotLabX = 45,
  scale = TRUE,
  cexMain = 1.5,  # Set title size here
  plotRsquared = TRUE,
  corFUN = 'pearson',
  corUSE = 'pairwise.complete.obs',
  signifSymbols = c('****', '***', '**', '*', ''),
  signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
  returnPlot = FALSE
)

pearson_eigencor

# Second plot for PC1-4 Clinical Correlations
correlate_eigencor <- eigencorplot(
  p,
  components = getComponents(p, 1:4),
  metavars = c("Age", "Sex", "Subtype", "LLV Volume", 
               "RLV Volume", "TV Volume", "FV Volume", 
               "LH Mean Thickness", "RH Mean Thickness", 
               "ab42/40", "ttau", "ptau", "NfL", "Albumin Quotient"),
  col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
  cexCorval = 0.7,
  colCorval = 'white',
  fontCorval = 2,
  posLab = 'bottomleft',
  rotLabX = 45,
  posColKey = 'right',
  cexLabColKey = 1.5,
  scale = TRUE,
  main = 'PC1-4 Clinical Correlations',
  colFrame = 'white',
  plotRsquared = FALSE,
  cexMain = 1.5  # Ensure title size is the same here
)

correlate_eigencor


library(cowplot)
library(ggplotify)

top_row <- plot_grid(pscree, pbiplot,
                     ncol = 2,
                     labels = c('A  Scree Plot', 'B  PCA Bi-plot'),
                     label_fontfamily = 'serif',
                     label_fontface = 'bold',
                     label_size = 18,
                     align = 'h',
                     rel_widths = c(0.8, 0.8))

bottom_row <- plot_grid(as.grob(pearson_eigencor),
                        as.grob(correlate_eigencor),
                        nrow = 2,
                        align = 'h')

plot_grid(top_row, ncol = 1,
          rel_heights = c(1, 1))

ggsave("plot_grid_output.pdf", plot = plot_grid(bottom_row, ncol = 1, rel_heights = c(1, 1)), 
       width = 10, height = 10, device = "pdf")

