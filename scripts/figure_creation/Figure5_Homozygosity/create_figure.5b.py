from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import ScalarFormatter
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
import sys
import os

broad_genomic_range = "chr8:105300001-105400000"
def plot_smoothed_genomic_data(filepath):
    """
    Reads tab-separated genomic data from a file, extracts the position, 
    applies a rolling mean to reduce volatility, and generates a smooth 
    line plot comparing case and control homozygosity percentages.
    
    Args:
        filepath (str): The path to the tab-separated input file.
    """
    
    # 1. File Reading
    if not os.path.exists(filepath):
        print(f"Error: Input file not found at '{filepath}'")
        return

    try:
        # Read the tab-separated file directly
        df = pd.read_csv(filepath, sep='\t')
    except Exception as e:
        print(f"Error reading data: {e}")
        return

    # 2. Data Transformation
    
    # Extract Position from the SNP string and convert percentage columns to numeric
    df.columns = df.columns.str.strip()
    df['Position'] = df['snp'].apply(lambda x: int(x.split('_')[1]))
    chr_num = df["snp"].iloc[0].split('_')[0].replace("chr","")

    case_col = 'case_percent_homozygosed'
    control_col = 'control_percent_homozygosed'
    
    df[case_col] = pd.to_numeric(df[case_col])
    df[control_col] = pd.to_numeric(df[control_col])
    
    # Drop any resulting NaNs and sort by Position
    df = df.dropna().sort_values(by='Position').reset_index(drop=True)

    if len(df) < 3:
        print("Error: Need at least 3 data points.")
        return
    
    # --- PRIMARY SMOOTHING: APPLY ROLLING MEAN (Increased Window Size) ---
    # Increased to 50 to provide a much cleaner, less stochastic curve.
    window_size = 50
    
    df['Smoothed_Case'] = df[case_col].rolling(window=window_size, min_periods=1, center=True).mean()
    df['Smoothed_Control'] = df[control_col].rolling(window=window_size, min_periods=1, center=True).mean()
    
    # 3. Setup Plotting Environment
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 6))

    # Get data arrays (using the original position and smoothed values)
    X = df['Position'].values

    # --- Plotting Cases (Red Straight Line Segments) ---
    # Plotting smoothed data directly with straight lines (no spline interpolation)
    ax.plot(X, df['Smoothed_Case'], color='#D93025', linewidth=3, label='Bladder Cases')
    
    # --- Plotting Controls (Gray Straight Line Segments) ---
    # Plotting smoothed data directly with straight lines (no spline interpolation)
    ax.plot(X, df['Smoothed_Control'], color='#60646B', linewidth=3, label='Controls')


    # 4. Finalize Chart Aesthetics
    ax.set_title(f'{broad_genomic_range.capitalize()}', fontsize=16)
    ax.set_xlabel(f"Chr{chr_num} Position", fontsize=12)
    ax.set_ylabel('Percent Homozygosed', fontsize=12)

    # Format Y-axis ticks as percentage
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0%}'))
    
    # Set X and Y limits
    ax.set_xlim(X.min() - 100, X.max() + 100)
    ax.set_ylim(0, 1) # STRICTLY set limits from 0 to 1 (0% to 100%)

    def comma_formatter(x, pos):
        return f"{int(x):,}"

    ax.xaxis.set_major_formatter(FuncFormatter(comma_formatter))
    ax.tick_params(axis='x', pad=8)   # increase pad (default is 4)

    # Highlight start and end coordinates
    highlight_start = 56750001
    highlight_end   = 56850000

    # Add shaded region
    ax.axvspan(
        highlight_start,
        highlight_end,
        color="lightblue",
        label="Bonferroni significant region",
        alpha=0.3   # adjust opacity (0 = transparent, 1 = solid)
    )

    ax.legend(loc='lower right')
    
    # 5. Save the figure
    output_filename = f"homozygosity_plot.{broad_genomic_range}.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot Saved Successfully to {output_filename}")


if __name__ == '__main__':
    
    if len(sys.argv) < 2:
        print("Usage: python plot_genotype.py <path_to_input_file.tsv>")
        sys.exit(1)

    input_filepath = sys.argv[1]
    plot_smoothed_genomic_data(input_filepath)
