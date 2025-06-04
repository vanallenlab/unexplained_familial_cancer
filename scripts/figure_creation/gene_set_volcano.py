import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import fdrcorrection
import argparse
import os

def plot_volano(data, ):
	df = pd.read_csv(data,sep='\t',index_col=False)

	# === Compute -log10(p-value) ===
	df['-log10_p'] = -np.log10(df['p_value'])

	# === Bonferroni threshold ===
	bonferroni_thresh = 0.05 / len(df)
	bonferroni_line = -np.log10(bonferroni_thresh)

	# === Nominal p-value threshold ===
	nominal_thresh = 0.05
	nominal_line = -np.log10(nominal_thresh)

	# === FDR threshold (5%) ===
	_, fdr_corrected = fdrcorrection(df['p_value'], alpha=0.05)
	df['fdr'] = fdr_corrected
	fdr_5_percent = max(df.loc[df['fdr'] <= 0.05, 'p_value'], default=1.0)
	fdr_line = -np.log10(fdr_5_percent)

	# === Volcano Plot ===
	plt.figure(figsize=(10, 7))
	sns.set(style='whitegrid')

	# Scatter points
	sns.scatterplot(data=df, x='beta', y='-log10_p', s=30, edgecolor=None, alpha=0.7)

	# Add lines
	plt.axhline(y=bonferroni_line, color='red', linestyle='--', label='Bonferroni (0.05)')
	plt.axhline(y=nominal_line, color='orange', linestyle='--', label='Nominal (0.05)')
	if fdr_5_percent < 1.0:
	    plt.axhline(y=fdr_line, color='blue', linestyle='--', label='FDR (5%)')

	# Labels and formatting
	plt.xlabel('Effect Size (Beta)', fontsize=12)
	plt.ylabel('-log10(p-value)', fontsize=12)
	plt.title('Volcano Plot', fontsize=14)
	plt.legend()
	plt.tight_layout()
	plt.save_fig(args.output_name)

def main()
    """
    Main block
    """
	parser = argparse.ArgumentParser(
			description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('-d','--data', required=True, description="A file with all p-values, betas, etc for each pathway")
	parser.add_argument('-o', '--output-name', required=False, default = "pathway_volcano.png", help = "Output file path")
    args = parser.parse_args()

    plot_volano(args.data,args.output_name)

if __name__ == "__main__":
	main()