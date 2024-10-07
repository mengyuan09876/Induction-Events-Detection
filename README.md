# Significant Induction Events Detection

This repository contains a script for detecting significant induction events between host MAGs (Metagenome-Assembled Genomes) and corresponding phages based on log-transformed abundance values. This tool was developed and used for the submitted paper, **"Plasticizers determine a deeper reshape of soil virome than microplastics."** It is particularly useful in microbial ecology studies for monitoring phage activity under various environmental conditions.

## Overview

The script follows these main steps:

1. **Log Transformation**  
   Abundance values are log-transformed to the form \(\log(\text{value} + 1)\) to handle zeros and small values effectively, preparing the data for reliable ratio calculation.

2. **Coverage Ratio Calculation**  
   For each host MAG and its corresponding phage (identified by matching identifiers), the log-transformed host abundance is subtracted from the log-transformed phage abundance to produce a log-transformed fold-change value. This step highlights changes in phage activity relative to the host.

3. **Threshold Application**  
   A threshold of `1` is applied to define significant induction events. A log-difference of \( \geq 1 \) indicates that the phage abundance is approximately 2.7 times higher than the host abundance. This criterion is used to distinguish true induction events from background variability, particularly under stress conditions like Diethyl Phthalate (DP) treatment. For more details, see the Supplementary Information (SI).

4. **Export Results**  
   Significant events that meet or exceed the threshold are exported to a CSV file for further analysis, allowing for easy integration with other bioinformatics pipelines or data visualization tools.

5. **Heatmap Generation**  
   The script generates a heatmap of the significant induction events, visually representing the log-transformed ratios of phage to host abundance. This helps to quickly identify patterns of induction across samples and phages. The heatmap is saved as a PNG file.

## Requirements

To use this script, you’ll need Python installed with the following packages:
- `numpy`
- `pandas`
- `seaborn`
- `matplotlib`

You can install the dependencies with:
```bash
pip install -r requirements.txt   
