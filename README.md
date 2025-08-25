## Outline of the Statistical Analysis Pipeline

### 1. Ventricular Volume Comparison Across Clinical Subgroups

- **Script:** `ventricular_volume_violins.R`  
- **Input Data:** `heatmap_data_norm.csv`  

This script analyzes and visualizes ventricular volumes across clinical subgroups (**AD**, **MCI_A+T+**, **MCI_A-T-**, **SCD**). The dataset is preprocessed to harmonize subgroup labels and apply custom color schemes. Violin plots with embedded boxplots and jittered points are generated for each ventricle (left lateral, right lateral, third, and fourth) as well as for the total ventricular volume, providing a clear comparison of distribution patterns across groups.

Pairwise **Wilcoxon rank-sum tests** are conducted between selected subgroups, with significance levels displayed as **star annotations** above the plots. All plots are combined into a multi-panel figure with consistent styling and saved as a high-resolution PNG. In addition, the script computes and prints **mean ventricular volumes** for each subgroup.

### 2. Cognitive Performance and Age Comparison Across Clinical Subgroups

- **Script:** `cognitive_performance_and_age_violins.R`  
- **Input Data:** `heatmap_data_norm_cogn.csv`  

This script analyzes and visualizes cognitive performance and age across clinical subgroups (**AD**, **MCI_A+T+**, **MCI_A-T-**, **SCD**). The dataset is preprocessed to harmonize subgroup labels and verify the presence of all required cognitive measures (MMSE, MoCA, KOD, RAVLT, RCF) and age. For each measure, violin plots with embedded boxplots and jittered points are generated. Pairwise Wilcoxon rank-sum tests are performed between subgroups and SCD, with significance levels displayed as star annotations above the plots. All results are combined into a multi-panel figure and saved as a high-resolution PNG, providing a clear comparison of cognitive and demographic differences across groups.

### 3. Principal Component Analysis (PCA) of CSF Proteomics and Clinical Variables  

- **Script:** `principal_component_analysis.R`  
- **Input Data:** `heatmap_data.csv`  

This script performs **principal component analysis (PCA)** on normalized CSF proteomics data, integrating clinical and imaging metadata (ventricle volumes, cortical thickness, tau/amyloid biomarkers, NfL, etc.). Data are quantile-normalized and variance-stabilized before PCA is applied. Results include:  
- A **PCA biplot** colored by subgroup and sex, with 95% confidence ellipses.  
- **Wilcoxon rank-sum tests** comparing subgroup differences along PC1 and PC3.  
- **Correlation heatmaps** showing associations between PCs and clinical variables (both r and r²).  

All plots are combined into a labeled multi-panel figure and saved as a high-resolution PNG. Statistical test results are also printed with FDR-adjusted p-values.  

### 4. Subgroup × Ventricle Interaction Effects on CSF Proteins  

- **Script:** `aqp4_interaction_model_plots.R`  
- **Input Data:** `heatmap_data.csv`  

This script investigates how **ventricular volumes interact with subgroup status (AD vs. SCD)** to influence CSF protein expression, focusing on selected proteins (e.g., AQP4, GAP43, DDAH1). Ventricle volumes are normalized, protein intensities are variance-stabilized, and **bootstrap resampling** is applied to estimate slope and R² stability for subgroup-specific regression models.  

For each protein–ventricle combination, the script generates **interaction plots** with subgroup regression lines, confidence intervals, and annotated bootstrap estimates (slopes and R² for AD and SCD). A multi-panel figure is produced for **AQP4**, summarizing effects across all ventricles, and saved as a high-resolution PNG.  

### 6. CSF Protein Comparison Across Clinical Subgroups  

- **Script:** `csf_protein_level_comparison.R`  
- **Input Data:** `heatmap_data.csv`  

This script compares **CSF protein levels** across clinical subgroups (**AD**, **MCI_A+T+**, **MCI_A-T-**, **SCD**). Selected proteins (AQP4, AMPH, DDAH1, SNCB, GAP43, ARPP21) are quantile- and variance-stabilized before visualization. For each protein, violin plots with embedded boxplots and jittered points are generated.  

Pairwise **Wilcoxon rank-sum tests** are performed between subgroups and SCD, and significance is displayed as star annotations above the plots. All protein plots are arranged into a labeled multi-panel figure and saved as a high-resolution PNG.  

### 7. SVM Classification of AD vs. SCD Using Protein and Ventricle Features  
- **Script:** `auc_density.R`  
- **Input Data:** `heatmap_data_norm.csv`  
This script builds **support vector machine (SVM)** models to classify AD vs. SCD from amyloid-/tau-associated proteins, with or without ventricular volumes. It evaluates four model families with 5-fold CV: (1) single-protein, (2) protein + ventricle, (3) protein-pair, (4) protein-pair + ventricle. For each configuration, **AUC** is computed from ROC curves. A **density plot** visualizes AUC distributions with **mean AUC** labeled; the figure is saved as a high-resolution PNG.

### 9. ROC Curve Analysis of Selected Protein and Ventricle Combinations  
- **Script:** `roc_plots.R`  
- **Input Data:** `heatmap_data_norm.csv`  
This script generates **ROC curves** for predefined feature sets to classify AD vs. SCD using SVMs with 5-fold CV.  
**Plot A:** PTPRN2 + GAP43; PTPRN2 + GAP43 + ToV; ToV alone; global protein sum (ToP).  
**Plot B:** GAP43 + ToV; PTPRN2 + ToV; GAP43 alone; PTPRN2 alone.  
It overlays ROC curves with labeled **AUC** values, arranges them into a two-panel figure, saves a high-resolution PNG, and prints AUC summaries to the console.

### 10. Cognitive Test Prediction Models Using CSF Proteins and Ventricular Volumes  
- **Script:** `multivariate_modeling_of_cognitive_test_results.R`  
- **Input Data:** `heatmap_data_norm_cogn.csv`  
This script models **cognitive scores** (**RAVLT, MMSE, MoCA**) using amyloid-/tau-associated proteins with/without ventricular volumes. It evaluates four regression model families: (1) protein only, (2) protein + ventricle, (3) protein pairs, (4) protein pairs + ventricle. Each model reports **cross-validated R² (boot-based)**, adjusted R², F-statistic, p-values, predictors, and sample size. Outputs include detailed and summary **CSV** files, a **comparison plot** (PNG), and a table of **top-performing models**.

