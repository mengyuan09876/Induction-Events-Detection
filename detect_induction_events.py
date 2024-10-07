# detect_induction_events.py

import numpy as np
import pandas as pd

def load_data(file_path):
    """
    Loads the abundance data from a CSV file.
    Args:
        file_path (str): The path to the CSV file containing abundance data.
    Returns:
        pd.DataFrame: The loaded abundance data.
    """
    return pd.read_csv(file_path, index_col=0)

def log_transform(data):
    """
    Applies log(x + 1) transformation to the abundance data.
    Args:
        data (pd.DataFrame): The abundance data.
    Returns:
        pd.DataFrame: Log-transformed abundance data.
    """
    return np.log1p(data)

def calculate_coverage_ratios(log_data):
    """
    Calculates the log-transformed coverage ratio between phages and MAGs.
    Args:
        log_data (pd.DataFrame): The log-transformed abundance data.
    Returns:
        pd.DataFrame: DataFrame with the coverage ratios.
    """
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
    """
    Identifies significant induction events based on a threshold.
    Args:
        coverage_ratios (pd.DataFrame): The coverage ratios.
        threshold (float): The threshold for significant induction.
    Returns:
        pd.DataFrame: Filtered DataFrame with significant induction events.
    """
    induction_events = coverage_ratios >= threshold
    return coverage_ratios[induction_events]

def save_results(data, output_path):
    """
    Saves the significant induction events to a CSV file.
    Args:
        data (pd.DataFrame): The DataFrame to save.
        output_path (str): The path for the output CSV file.
    """
    data.to_csv(output_path, index=True)
    print(f"Results saved to {output_path}")

def main(file_path, output_path, threshold=1):
    """
    Main function to run the induction event detection process.
    Args:
        file_path (str): The path to the CSV file containing abundance data.
        output_path (str): The path for the output CSV file.
        threshold (float): The threshold for significant induction.
    """
    abundance_data = load_data(file_path)
    log_data = log_transform(abundance_data)
    coverage_ratios = calculate_coverage_ratios(log_data)
    significant_events = identify_induction_events(coverage_ratios, threshold)
    save_results(significant_events, output_path)

if __name__ == "__main__":
    # Example usage
    # Replace 'abundance_data.csv' and 'significant_induction_events.csv' with your actual file paths
    main('abundance_data.csv', 'significant_induction_events.csv', threshold=1)
