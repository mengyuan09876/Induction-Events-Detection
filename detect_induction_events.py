# detect_induction_events.py

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

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

def generate_heatmap(input_path, output_path="heatmap.png"):
    """Generates a heatmap from the significant induction events CSV file and saves it as an image."""
    # Load the significant induction events from the saved CSV file
    data = pd.read_csv(input_path, index_col=0)
    
    # Create the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(data.fillna(0), cmap="YlGnBu", annot=True, fmt=".2f", cbar_kws={'label': 'Log Ratio'})
    plt.title("Significant Induction Events Heatmap")
    plt.xlabel("Phages")
    plt.ylabel("Samples")
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Heatmap saved to {output_path}")
    plt.show()

def main(input_file, output_file, threshold=1, heatmap_file="heatmap.png"):
    """Main function to run the induction event detection process."""
    abundance_data = load_data(input_file)
    log_data = log_transform(abundance_data)
    coverage_ratios = calculate_coverage_ratios(log_data)
    significant_events = identify_induction_events(coverage_ratios, threshold)
    save_results(significant_events, output_file)
    
    # Generate heatmap based on the saved CSV file for significant induction events
    generate_heatmap(output_file, heatmap_file)

if __name__ == "__main__":
    # Set up argument parsing for input parameters
    parser = argparse.ArgumentParser(description="Detect and visualize induction events between MAGs and phages.")
    parser.add_argument('--input', type=str, required=True, help="Path to the input CSV file with abundance data.")
    parser.add_argument('--output', type=str, default='significant_induction_events.csv', help="Path to save the significant induction events CSV file.")
    parser.add_argument('--threshold', type=float, default=1, help="Threshold for detecting induction events.")
    parser.add_argument('--heatmap', type=str, default='heatmap.png', help="Path to save the heatmap image.")

    # Parse the arguments
    args = parser.parse_args()

    # Run the main function with parsed arguments
    main(args.input, args.output, args.threshold, args.heatmap)
