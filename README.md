# Significant Induction Events Detection

This repository contains a script for detecting significant induction events between host MAGs (Metagenome-Assembled Genomes) and corresponding phages based on log-transformed abundance values. This tool was developed and used for the submitted paper, **"Plasticizers determine a deeper reshape of soil virome than microplastics."** It is particularly useful in microbial ecology studies for monitoring phage activity under various environmental conditions.

## Overview

The script follows these main steps:

1. **Log Transformation**  
   Abundance values are log-transformed to the form \(\log(\text{value} + 1)\) to handle zeros and small values effectively, preparing the data for reliable ratio calculation.

2. **Coverage Ratio Calculation**  
   For each host MAG and its corresponding phage (identified by matching identifiers), the log-transformed host abundance is subtracted from the log-transformed phage abundance to produce a log-transformed fold-change value. In this context, **phage abundances are provided as RPKM values**, while **host abundances are represented as relative abundances**. 
3. **Threshold Application**  
   A threshold of `1` is applied to define significant induction events. A log-difference of \( \geq 1 \) indicates that the phage abundance is approximately 2.7 times higher than the host abundance. This criterion is used to distinguish true induction events from background variability.

4. **Export Results**  
   Significant events that meet or exceed the threshold are exported to a CSV file for further analysis, allowing for easy integration with other bioinformatics pipelines or data visualization tools.

5. **Heatmap Generation**  
   The script generates a heatmap of the significant induction events, visually representing the log-transformed ratios of phage to host abundance. This helps to quickly identify patterns of induction across samples and phages. The heatmap is saved as a PNG file.

##Example Data
Below is an example of the input data format expected by the script:

| H_V          | BL_a      | BL_b      | BL_c      | PE_a      | PE_b      | PE_c      | PVC_a     | PVC_b     | PVC_c     | DP_a        | DP_b       | DP_c       |
|--------------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-------------|------------|------------|
| M_sp.upd101  | 0.00155   | 0.00523   | 0.00572   | 0.00163   | 0.00458   | 0.00264   | 0.00165   | 0.00187   | 0.00340   | 32.24074    | 20.84618   | 20.84618   |
| V_sp.upd101_1| 0         | 0         | 0         | 0         | 0         | 0         | 0         | 0         | 0         | 304.03882   | 557.16807  | 557.16807  |
| V_sp.upd101_2| 0         | 0         | 0         | 0         | 0         | 0         | 0         | 0         | 0         | 0           | 122.70457  | 122.70457  |
| M_sp.upd105  | 0.01631   | 0.00071   | 0.00257   | 0.00047   | 0.00415   | 0.00536   | 0.00445   | 0.00592   | 0.00056   | 0.01382     | 0.64066    | 0.64066    |
| V_sp.upd105_1| 0         | 0         | 0         | 0         | 0         | 0         | 0         | 0         | 0         | 0           | 138.5803   | 138.5803   |
| V_sp.upd105_2| 0         | 0         | 0         | 0         | 0         | 0         | 0         | 0         | 0         | 92.3078     | 127.48879  | 127.48879  |

- **MAG Rows (`M_`)**: Represent the relative abundance of the host MAG.
- **Phage Rows (`V_`)**: Represent the RPKM values for corresponding phages.
- **Conditions**: Each treatment (e.g., BL, PE, PVC, DP) has three replicates (`_a`, `_b`, `_c`).
  
## Requirements

To use this script, you’ll need Python installed with the following packages:
- `numpy`
- `pandas`
- `seaborn`
- `matplotlib`
  

You can install the dependencies with:
```bash
pip install -r requirements.txt
