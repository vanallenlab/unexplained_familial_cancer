!gsutil -m cp gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/Assess_RF_3/3ac28007-f29e-4da1-8a6a-d2d390537757/call-PlotHWEvsTPProb_snps/tp_prob_hwe_output_snp.tsv hwe_snps.tsv
!gsutil -m cp gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/Assess_RF_2/b2ea7d39-ca72-4b60-891a-1a6055560759/call-plot_sensitivity_snps/sensitivity_output_snp.tsv sensitivity_snps.tsv
!gsutil -m cp gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/Assess_RF/2fabefbf-4aec-4df1-8f19-8b54426d5f15/call-GraphTotalSnps/tp_prob_vs_median_SNPs.tsv total_snps.tsv
!gsutil -m cp gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/Assess_RF/2fabefbf-4aec-4df1-8f19-8b54426d5f15/call-GraphDeNovoSnps/tp_prob_vs_median_DeNovo_SNPs.tsv denovo_snps.tsv


import pandas as pd
import numpy as np

# Load the data from the TSV files
hwe_df = pd.read_csv('hwe_snps.tsv', sep='\t')
sensitivity_df = pd.read_csv('sensitivity_snps.tsv', sep='\t')
denovo_df = pd.read_csv('denovo_snps.tsv', sep='\t')
total_df = pd.read_csv('total_snps.tsv', sep='\t')

# Merge all dataframes on 'tp_prob'
merged_df = pd.merge(hwe_df, sensitivity_df, on='tp_prob')
merged_df = pd.merge(merged_df, denovo_df, on='tp_prob')
merged_df = pd.merge(merged_df, total_df, on='tp_prob')

# Normalize the 'total_snps' column by dividing by the max value
merged_df['median_SNPs_norm'] = merged_df['median_SNPs_variants'] / merged_df['median_SNPs_variants'].max()
merged_df['median_DeNovo_SNPs_norm'] = merged_df['median_DeNovo_SNPs_variants'] / merged_df['median_DeNovo_SNPs_variants'].max()
merged_df['median_p_hwe_1_percent_euro_norm'] = merged_df['median_p_hwe_1_percent_euro'] / merged_df['median_p_hwe_1_percent_euro'].max()

# Define columns to optimize
max_cols = ['median_SNPs_norm', 'sensitivity']  # Columns to maximize
min_cols = ['median_DeNovo_SNPs_variants', 'median_p_hwe_1_percent_euro']  # Columns to minimize

#merged_df.head()
pd.set_option('display.max_rows', None)
# Calculate the Euclidean distance for each row
merged_df['euclidean_distance'] = np.sqrt(
    (1-merged_df['median_SNPs_norm']) ** 2 + (1- merged_df['sensitivity']) ** 2 + (merged_df['median_DeNovo_SNPs_norm']) ** 2 + (merged_df['median_p_hwe_1_percent_euro_norm']) ** 2
)
merged_df['euclidean_distance_sensitivity'] = np.sqrt(
    (1- merged_df['sensitivity']) ** 2 + (1-merged_df['median_SNPs_norm']) ** 2
)
merged_df['euclidean_distance_specificity'] = np.sqrt(
    (merged_df['median_DeNovo_SNPs_norm']) ** 2 + (merged_df['median_p_hwe_1_percent_euro_norm']) ** 2
)