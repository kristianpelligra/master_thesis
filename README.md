## Outline of the Statistical Analysis Pipeline

### 1. Ventricular Volume Comparison Across Clinical Subgroups

- **Script:** `ventricular_volume_violins.R`  
- **Input Data:** `heatmap_data_norm.csv`  

This script analyzes and visualizes ventricular volumes across clinical subgroups (**AD**, **MCI_A+T+**, **MCI_A-T-**, **SCD**). The dataset is preprocessed to harmonize subgroup labels and apply custom color schemes. Violin plots with embedded boxplots and jittered points are generated for each ventricle (left lateral, right lateral, third, and fourth) as well as for the total ventricular volume, providing a clear comparison of distribution patterns across groups.

Pairwise **Wilcoxon rank-sum tests** are conducted between selected subgroups, with significance levels displayed as **star annotations** above the plots. All plots are combined into a multi-panel figure with consistent styling and saved as a high-resolution PNG. In addition, the script computes and prints **mean ventricular volumes** for each subgroup.
