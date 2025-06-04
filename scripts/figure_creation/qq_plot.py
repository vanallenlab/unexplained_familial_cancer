from qqman import qqman
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def return_pvals(file, max_maf, annotation, column_of_interest = "Pvalue"):
	df = pd.read_csv(file, sep='\t')
	df = df[(df['Group'] == annotation) & (df['max_MAF'] == max_maf)]
	p_vals = list(df[column_of_interest])
	return p_vals

def main():
	"""
    Main block
    """
	parser = argparse.ArgumentParser(
    	description=__doc__,
    	formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('-d', '--data', required=True, help = 'SAIGE Output')
	parser.add_argument('-o', '--output', required=False, default = 'saige_results.png', help='output png dir') 
	args = parser.parse_args()

	### FILES ###
	df = pd.read_csv(args.data, sep='\t', index_col=False)

    # Get dimensions for our analysis
	max_mafs = sorted(list(set(df['max_MAF'])))
	nrows = len(max_mafs)

	groups = sorted(list(set(df['Group'])))
	ncols = len(groups)

	# Base width and height per subplot (in inches)
	base_width = 4
	base_height = 3

	# Adjust figsize based on grid shape
	fig_width = base_width * ncols
	fig_height = base_height * nrows


    ### QQ PLOT ###
	figure, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize = (fig_width,fig_height))
	axes = np.atleast_2d(axes)
	fontsize = 10
	for i in range(nrows):
		for j in range (ncols):
			af = max_mafs[i]
			grouping = groups[j]
			p_vals = return_pvals(args.data, af, annotation=grouping)
			qqman.qqplot(p_vals, ax=axes[i, j], title=f"{grouping} <={af * 100}%")
			axes[i, j].set_title(f"{grouping} <={af * 100}%", fontsize=fontsize)


	figure.tight_layout()
	plt.savefig(args.output,format="png")
	plt.close()

if __name__ == "__main__":
	main()