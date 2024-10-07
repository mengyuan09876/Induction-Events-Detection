# detect_induction_events.py

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def load_data(file_path):
    """Loads the abundance data from a CSV file."""
    return pd.read_csv(file_path, index_col=0)

def log_transform(data):
    """Applies log(x + 1) transformation to the abundance data."""
    return np.log1p(data)

def calculate_coverage_ratios(log_data):
    """Calculates the log-transformed coverage ratio between phages and MAGs."""
    mag_columns = [col for col in log_data.columns if col.startswith('M_')]
    phage_columns = [col for col in log_data.columns if col.startswith('V_')]
    
    coverage_ratios = {}
    for mag_col in mag_columns:
        identifier = mag_col.split('_')[1]
        phage_col = f"V_{identifier}"
        
        if phage_col in phage_columns:
            coverage_ratios[phage_col] = log_data[phage_col] - log_data[mag_col]
    
    return pd.DataFrame(coverage_ratios)

def identify_induction_events(coverage_ratios, threshold=1):
    """Identifies significant induction events based on a threshold."""
    induction_events = coverage_ratios >= threshold
    return coverage_ratios[induction_events]

def save_results(data, output_path):
    """Saves the significant induction events to a CSV file."""
    data.to_csv(output_path, index=True)
    print(f"Results saved to {output_path}")

def generate_heatmap(data, output_path="heatmap.png"):
    """Generates a heatmap from the significant induction events and saves it as an image."""
    plt.figure(figsize=(10, 8))
    sns.heatmap(data, cmap="YlGnBu", annot=True, fmt=".2f", cbar_kws={'label': 'Log Ratio'})
    plt.title("Significant Induction Events Heatmap")
    plt.xlabel("Phages")
    plt.ylabel("Samples")
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Heatmap saved to {output_path}")
    plt.show()

def main(file_path, output_path, threshold=1, heatmap_path="heatmap.png"):
    """Main function to run the induction event detection process."""
    abundance_data = load_data(file_path)
    log_data = log_transform(abundance_data)
    coverage_ratios = calculate_coverage_ratios(log_data)
    significant_events = identify_induction_events(coverage_ratios, threshold)
    save_results(significant_events, output_path)
    
    # Generate heatmap if there are any significant events
    if not significant_events.empty:
        generate_heatmap(significant_events.fillna(0), heatmap_path)

if __name__ == "__main__":
    # Replace 'abundance_data.csv' and 'significant_induction_events.csv' with your actual file paths
    main('abundance_data.csv', 'significant_induction_events.csv', threshold=1)
