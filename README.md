Here is an outline of the statistical test pipeline

Ventricular Volume Comparison Across Clinical Subgroups


This R script analyzes and visualizes ventricular volumes across clinical subgroups (AD, MCI_A+T+, MCI_A-T-, SCD). It preprocesses the dataset, standardizes subgroup labels, and assigns custom colors. For each ventricle (left, right, third, fourth, and total), it generates violin plots with embedded boxplots and jittered data points. Pairwise Wilcoxon rank-sum tests are performed between selected subgroups, and significance levels are displayed as star annotations above the plots. All plots are combined into a multi-panel figure with consistent styling and saved as a high-resolution PNG. The script also calculates and prints mean ventricular volumes for each subgroup.
