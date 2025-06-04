# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0

workflow STEP_7_VARIANT_HWE {
  input {
    String step_4_output_dir       # Directory to STEP_4 Output VCFs
 
    File rf_predictions	# File output from hail RF in step 5
    File euro_sample_list # A list of unrelated europeans in the study for HWE analysis
  }

  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call gather_vcfs {
    input:
      dir = step_4_output_dir
  }

  call sort_vcf_list {
    input:
      unsorted_vcf_list = gather_vcfs.vcf_list
  }

  call gather_tsvs{}

  #Array[String] rf_prediction_tsvs = read_lines(gather_tsvs.out)

  scatter (i in range(length(sort_vcf_list.vcf_arr)-0)){
    String index_str = i

    # Remove all variants where "*" is in the ALT allele 
    call Prepare_VCF {
      input:
        vcf = sort_vcf_list.vcf_arr[i]
    }
    call Extract_Variants_From_VCF {
      input:
        vcf = Prepare_VCF.out_vcf
    }
    call RunVcfTools {
      input:
        vcf_slice = Prepare_VCF.out_vcf,
        callset_name = index_str,
        euro_sample_list = euro_sample_list
    }
    call JoinHWEandVariants {
      input:
        variants = Extract_Variants_From_VCF.variants,
        hwe_file = RunVcfTools.hwe_file,
        af_file = RunVcfTools.af_file,
        hwe_file_euro = RunVcfTools.hwe_file_euro,
        af_file_euro = RunVcfTools.af_file_euro,
        rf_predictions = gather_tsvs.out[i]
        #rf_predictions = rf_predictions
    }
  }
  call concatenateFiles as concatenateFiles_snps{
    input:
      files = JoinHWEandVariants.merged_output_snps,
      callset_name = "ufc_wgs_snps"
  }
  call concatenateFiles as concatenateFiles_indels{
    input:
      files = JoinHWEandVariants.merged_output_indels,
      callset_name = "ufc_wgs_indels"
  }
  call PlotHWEvsTPProb as PlotHWEvsTPProb_snps {
    input:
      merged_file = concatenateFiles_snps.out,
      snp_or_indel = "snp"
  }
  call PlotHWEvsTPProb as PlotHWEvsTPProb_indels {
    input:
      merged_file = concatenateFiles_indels.out,
      snp_or_indel = "indel"
  }
}


task gather_tsvs {
  input {
  }
  command <<<
  python3 <<CODE
  prefix = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_4_RANDOM_FOREST/sharded_tsvs/rf_predictions."
  suffix = ".tsv.bgz"

  with open("rf_predictions.sharded.list", "w") as f:
    for i in range(3103):  # 0 to 3102 inclusive
      f.write(f"{prefix}{i}{suffix}\n")
  CODE
  >>>
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
  output {
    Array[String] out = read_lines("rf_predictions.sharded.list")
  }
}

task Prepare_VCF {
  input {
    File vcf
  }
  command <<<
  set -euxo pipefail
  bcftools view -e 'ALT="*"' -Oz -o nostar.vcf.gz ~{vcf} 
  >>>
  runtime{
    docker:"vanallenlab/g2c_pipeline:latest"
    preemptible:3
  }
  output{
    File out_vcf = "nostar.vcf.gz"
  }
}

task gather_vcfs {
  input {
    String dir
  }
  command <<<
  set -euxo pipefail
  gsutil ls ~{dir}/*.vcf.bgz > vcf.list
  >>>
  output {
    File vcf_list = "vcf.list"
  }
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task sort_vcf_list {
  input {
    File unsorted_vcf_list
  }
  command <<<
  set -eu -o pipefail
  python3 <<CODE
  import re

  # Input file containing the list of file paths
  file_path = '~{unsorted_vcf_list}'

  # Read and parse file paths
  with open(file_path, 'r') as f:
      paths = f.readlines()

  # Function to extract chromosome and shard number for sorting
  def extract_key(path):
      # Extract chromosome and shard numbers using regex
      match = re.search(r'chr([0-9XY]+)-finalrun\.([0-9]+)', path)
      if match:
          chrom_str, shard_str = match.groups()

          # Convert chromosome to integer, treating 'X' as 23 and 'Y' as 24
          chrom_num = 23 if chrom_str == 'X' else 24 if chrom_str == 'Y' else int(chrom_str)
          shard_num = int(shard_str)

          # Return a tuple with chromosome and shard for sorting
          return (chrom_num, shard_num)
      else:
          return (float('inf'), float('inf'))  # Unmatched lines go to the end
  # Sort paths using the extracted keys
  sorted_paths = sorted(paths, key=extract_key)

  # Write to vcf.sorted.list and vcf_idx.sorted.list
  with open('vcf.sorted.list', 'w') as vcf_file, open('vcf_idx.sorted.list', 'w') as vcf_idx_file:
      for path in sorted_paths:
          clean_path = path.strip()
          vcf_file.write(clean_path + '\n')
          vcf_idx_file.write(clean_path + '.tbi\n')
  CODE
  >>>
  output {
    Array[String] vcf_arr = read_lines("vcf.sorted.list")
    Array[String] vcf_idx_arr = read_lines("vcf_idx.sorted.list")
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task PlotHWEvsTPProb {
  input {
    File merged_file  # The merged file containing p_hwe and tp_prob
    String snp_or_indel
  }

  command <<<
  python3 <<CODE
  import pandas as pd
  import matplotlib.pyplot as plt
  import numpy as np

  # Define the desired column names
  headers = ['CHROM', 'POS', 'P_HWE', 'P_HWE_EURO', 'REF_FREQ', 'ALT_FREQ', 'REF_FREQ_EURO', 'ALT_FREQ_EURO', 'allele_type', 'rf_probability', 'fp_prob', 'tp_prob', 'negative_log_p_hwe', 'negative_log_p_hwe_euro']

  # Read in the data
  df = pd.read_csv("~{merged_file}", header=None, names=headers, sep='\t')

  # 1. Calculate the number of lines in the file
  num_lines = len(df)
  print(f"Number of lines in the file: {num_lines}")

  # 2. Calculate Bonferroni significance
  alpha = 0.05
  bonferroni_significance = alpha / num_lines
  print(f"Bonferroni significance threshold: {bonferroni_significance}")

  # Function to calculate percentage of points passing Bonferroni significance
  def calculate_percentage(df, bonferroni_significance, tp_prob_values, column='P_HWE'):
      percentages = []
      for tp in tp_prob_values:
          passing_points = df[df['tp_prob'] >= tp][column] <= bonferroni_significance
          percentage = (passing_points.sum() / len(df[df['tp_prob'] >= tp])) * 100 if len(df[df['tp_prob'] >= tp]) > 0 else 0
          percentages.append(percentage)
      return percentages

  # 3. Prepare to calculate percentages for different AF thresholds
  tp_prob_values = np.arange(0, 1.01, 0.01)

  # No AF threshold (use full data)
  percentages_no_threshold = calculate_percentage(df, bonferroni_significance, tp_prob_values)
  percentages_no_threshold_euro = calculate_percentage(df, bonferroni_significance, tp_prob_values, column='P_HWE_EURO')

  # AF ≥ 1%
  df_af_1_percent = df[(df['ALT_FREQ'] >= 0.01) & (df['REF_FREQ'] >= 0.01)]
  bonferroni_significance_1_percent = alpha / len(df_af_1_percent)
  percentages_1_percent = calculate_percentage(df_af_1_percent, bonferroni_significance_1_percent, tp_prob_values)
  df_af_1_percent_euro = df[(df['REF_FREQ_EURO'] >= 0.01) & (df['ALT_FREQ_EURO'] >= 0.01)]
  bonferroni_significance_1_percent_euro = alpha / len(df_af_1_percent_euro)
  percentages_1_percent_euro = calculate_percentage(df_af_1_percent_euro, bonferroni_significance_1_percent_euro, tp_prob_values, column='P_HWE_EURO')

  # AF ≥ 5%
  df_af_5_percent = df[(df['ALT_FREQ'] >= 0.05) & (df['REF_FREQ'] >= 0.05)]
  bonferroni_significance_5_percent = alpha / len(df_af_5_percent)
  percentages_5_percent = calculate_percentage(df_af_5_percent, bonferroni_significance_5_percent, tp_prob_values)
  df_af_5_percent_euro = df[(df['REF_FREQ_EURO'] >= 0.05) & (df['ALT_FREQ_EURO'] >= 0.05)]
  bonferroni_significance_5_percent_euro = alpha / len(df_af_5_percent_euro)
  percentages_5_percent_euro = calculate_percentage(df_af_5_percent_euro, bonferroni_significance_5_percent_euro, tp_prob_values, column='P_HWE_EURO')

  # 4. Plot the results with different colors for each AF threshold
  plt.figure(figsize=(10, 6))

  # Non-European plots
  plt.plot(tp_prob_values, percentages_no_threshold, linestyle='-', color='b', alpha=0.7, label='No AF Threshold')
  plt.plot(tp_prob_values, percentages_1_percent, linestyle='-', color='r', alpha=0.7, label='AF ≥ 1%')
  plt.plot(tp_prob_values, percentages_5_percent, linestyle='-', color='g', alpha=0.7, label='AF ≥ 5%')

  # European plots
  plt.plot(tp_prob_values, percentages_no_threshold_euro, linestyle='--', color='b', alpha=0.7, label='No AF Threshold (Euro)')
  plt.plot(tp_prob_values, percentages_1_percent_euro, linestyle='--', color='r', alpha=0.7, label='AF ≥ 1% (Euro)')
  plt.plot(tp_prob_values, percentages_5_percent_euro, linestyle='--', color='g', alpha=0.7, label='AF ≥ 5% (Euro)')

  # Add a horizontal line at y=10
  plt.axhline(y=10, color='black', linestyle='--', linewidth=1, label='y = 10')

  # Add labels, title, and legend
  plt.xlabel('TP_Prob')
  plt.ylabel('Percentage of Biallelic Variants Failing HWE at Bonferonni Significance')
  plt.title('Percentage of Biallelic ~{snp_or_indel} Failing HWE vs TP_Prob')

  # Add legend for clarity
  plt.legend()

  # Save the plot
  plt.savefig("p_hwe_vs_tp_prob_plot_~{snp_or_indel}.png")

  # Create a DataFrame from the lists
  df_tp_prob_hwe = pd.DataFrame({
      'tp_prob': tp_prob_values,
      'median_p_hwe': percentages_no_threshold,
      'median_p_hwe_1_percent': percentages_1_percent,
      'median_p_hwe_5_percent': percentages_5_percent,
      'median_p_hwe_euro': percentages_no_threshold_euro,
      'median_p_hwe_1_percent_euro': percentages_1_percent_euro,
      'median_p_hwe_5_percent_euro': percentages_5_percent_euro
  })

  # Save the DataFrame to a TSV file
  df_tp_prob_hwe.to_csv("tp_prob_hwe_output_~{snp_or_indel}.tsv", sep='\t', index=False)
  CODE
  >>>

  output {
    File plot = "p_hwe_vs_tp_prob_plot_~{snp_or_indel}.png"  # Output the plot
    File output_file = "tp_prob_hwe_output_~{snp_or_indel}.tsv" # output file
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
    memory: "50 GB"
    disks: "local-disk 100 HDD"
  }
}

task JoinHWEandVariants {
  input {
    File variants  # The variants file
    File hwe_file  # The HWE file
    File af_file   # allele-frequency file
    File hwe_file_euro  # The HWE file
    File af_file_euro   # allele-frequency file
    File rf_predictions
  }

  command <<<
  set -euxo pipefail
  zcat ~{rf_predictions} > rf_predictions.tsv
  awk -F'\t' '!($2 ~ /\*/) {print}' rf_predictions.tsv > rf_predictions.nostar.tsv
  if [ $(wc -l < rf_predictions.nostar.tsv) -eq 0 ]; then
    touch output_snps.tsv
    touch output_indels.tsv
    exit 0
  fi 
  python3 <<CODE
  import pandas as pd
  import re
  import numpy as np
  import gc

  print("Checkpoint 0")
  # Read in AF Files
  df_af = pd.read_csv("~{af_file}",sep='\t',names=['CHROM','POS','N_ALLELES','N_CHR','REF_FREQ','ALT_FREQ'], skiprows=1)
  df_af_euro = pd.read_csv("~{af_file_euro}",sep='\t',names=['CHROM','POS','N_ALLELES_EURO','N_CHR_EURO','REF_FREQ_EURO','ALT_FREQ_EURO'], skiprows=1)

  # Filter and Merge AF Files
  df_af = df_af[['CHROM','POS','REF_FREQ','ALT_FREQ']]
  df_af = df_af[df_af.duplicated(subset=['CHROM', 'POS'], keep=False) == False]

  df_af_euro = df_af_euro[['CHROM','POS','REF_FREQ_EURO','ALT_FREQ_EURO']]
  df_af_euro = df_af_euro[df_af_euro.duplicated(subset=['CHROM', 'POS'], keep=False) == False]

  # Read in HWE Files
  df_hwe = pd.read_csv("~{hwe_file}",sep='\t',names=["CHROM","POS","OBS(HOM1/HET/HOM2)","E(HOM1/HET/HOM2)","ChiSq_HWE","P_HWE","P_HET_DEFICIT","P_HET_EXCESS"], skiprows=1)
  df_hwe = df_hwe[['CHROM','POS','P_HWE']]
  df_hwe = df_hwe[df_hwe.duplicated(subset=['CHROM', 'POS'], keep=False) == False]

  df_hwe_euro = pd.read_csv("~{hwe_file_euro}",sep='\t',names=["CHROM","POS","OBS(HOM1/HET/HOM2)_EURO","E(HOM1/HET/HOM2)_EURO","ChiSq_HWE_EURO","P_HWE_EURO","P_HET_DEFICIT_EURO","P_HET_EXCESS_EURO"], skiprows=1)
  df_hwe_euro = df_hwe_euro[['CHROM','POS','P_HWE_EURO']]
  df_hwe_euro = df_hwe_euro[df_hwe_euro.duplicated(subset=['CHROM', 'POS'], keep=False) == False]

  merged_hwe = pd.merge(df_hwe,df_hwe_euro, on=['CHROM','POS'], how='left')
  merged_af = pd.merge(df_af, df_af_euro, on=['CHROM','POS'], how='left')
  merged_allele_features = pd.merge(merged_hwe,merged_af, on=['CHROM','POS'], how="outer")
  merged_allele_features['POS'] =  merged_allele_features['POS'].astype(int)
  merged_allele_features['CHROM'] =  merged_allele_features['CHROM'].astype(str)
  print("Checkpoint 1")

  ### Phase 2 ###

  # Fix the 'alleles' column by replacing the extra double quotes
  df_rf_predictions = pd.read_csv("rf_predictions.nostar.tsv",sep='\t')
  #df_rf_predictions['n_alt_alleles'] = pd.to_numeric(df_rf_predictions['n_alt_alleles'], errors='coerce')
  #df_rf_predictions = df_rf_predictions[df_rf_predictions['n_alt_alleles'] == 1]

  # Convert 'rf_probability' column to strings (this ensures all values are treated as strings)
  df_rf_predictions['rf_probability'] = df_rf_predictions['rf_probability'].astype(str)
  df_rf_predictions['locus'] = df_rf_predictions['locus'].astype(str)
  df_rf_predictions[['CHROM', 'POS']] = df_rf_predictions['locus'].str.split(':', expand=True)
  df_rf_predictions['POS'] = df_rf_predictions['POS'].astype(int)
  df_rf_predictions['CHROM'] = df_rf_predictions['CHROM'].astype(str)
  df_rf_predictions = df_rf_predictions[['CHROM','POS','allele_type','rf_probability']]
  counts = df_rf_predictions.groupby(['CHROM', 'POS']).CHROM.transform('count')
  df_rf_predictions = df_rf_predictions[counts == 1]
  #df_rf_predictions['alleles'] = df_rf_predictions['alleles'].astype(str)

  # Merge the DataFrames on 'locus' and 'alleles'
  # Adjust column names if necessary based on the actual names in your files
  df = pd.merge(merged_allele_features, df_rf_predictions, on=['CHROM', 'POS'], how='left')
  print("df allele features info")
  print(merged_allele_features.head())
  print(len(merged_allele_features))
  print("rf_predictions tsv head")
  print(df_rf_predictions.head())
  print(len(df_rf_predictions))

  del df_rf_predictions
  del merged_allele_features
  gc.collect()

  df['P_HWE'] = df['P_HWE'].astype(float)
  df['P_HWE_EURO'] = df['P_HWE_EURO'].astype(float)
  df['REF_FREQ'] = df['REF_FREQ'].str.split(':').str[1].astype(float)
  df['ALT_FREQ'] = df['ALT_FREQ'].str.split(':').str[1].astype(float)
  df['REF_FREQ_EURO'] = df['REF_FREQ_EURO'].str.split(':').str[1].astype(float)
  df['ALT_FREQ_EURO'] = df['ALT_FREQ_EURO'].str.split(':').str[1].astype(float)
  df['rf_probability'] = df['rf_probability'].astype(str)

  print("Successful Merge!")

  # Function to extract fp_prob and tp_prob from the string
  def extract_probs(prob_str):
      # Use regular expressions to extract the numeric values for 'False' and 'True'
      fp_match = re.search(r'"False","value":([0-9.]+)', prob_str)
      tp_match = re.search(r'"True","value":([0-9.]+)', prob_str)

      # Extract the values, default to None if not found
      fp_prob = float(fp_match.group(1)) if fp_match else None
      tp_prob = float(tp_match.group(1)) if tp_match else None

      return fp_prob, tp_prob

  # Apply the function to each row in the 'rf_probability' column
  df['fp_prob'], df['tp_prob'] = zip(*df['rf_probability'].apply(extract_probs))
  print("Successful creation of tp_prob and fp_prob columns")

  # Ensure that 'tp_prob' is a float (it should already be a float, but just to be safe)
  df['tp_prob'] = df['tp_prob'].astype(float)

  # make a new column
  df['negative_log_p_hwe'] = -np.log10(df['P_HWE'])
  df['negative_log_p_hwe_euro'] = -np.log10(df['P_HWE_EURO'])
  print("Successful creation of negative log column")

  df_snps = df[df['allele_type'] == "snp"]
  df_indels = df[df['allele_type'] != "snp"]

  # Save the merged DataFrame to a file, if needed
  df_snps.to_csv("output_snps.tsv", sep='\t', index=False, header=False)
  df_indels.to_csv("output_indels.tsv",sep='\t',index=False, header=False)
  print("saved to output.tsv!")
  CODE
  >>>

  output {
    File merged_output_snps = "output_snps.tsv"  # Output the merged file
    File merged_output_indels = "output_indels.tsv"  # Output the merged file
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
    #memory: "G"
    #disks: "local-disk 100 HDD"
  }
}

task Extract_Variants_From_VCF {
  input{
    File vcf
  }
  command <<<
   set -euxo pipefail
   bcftools query -f '%CHROM:%POS\t\["%REF","%ALT"\]\n' ~{vcf}  > variants.tsv
  >>>
  output {
    File variants = "variants.tsv"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
    #disks: "local-disk 32 HDD"
    #memory: "8GB"
  }
}

task RunVcfTools {
  input {
    File vcf_slice
    File euro_sample_list
    String callset_name
  }

  command <<<
    # Extract Allele Frequency
    vcftools --gzvcf ~{vcf_slice} --freq --out ~{callset_name}

    #Gather HWE Statistics
    vcftools --gzvcf ~{vcf_slice} --hardy --out ~{callset_name}

    # Extract Allele Frequency for unrelated europeans
    vcftools --gzvcf ~{vcf_slice} --keep ~{euro_sample_list} --freq --out "~{callset_name}_euro"

    #Gather HWE Statistics for unrelated europeans
    vcftools --gzvcf ~{vcf_slice} --keep ~{euro_sample_list} --hardy --out "~{callset_name}_euro"
  >>>

  runtime{
    disks: "local-disk 20 HDD"
    #memory: "16GB"
    docker: "biocontainers/vcftools:0.1.15"
    preemptible: 3
  }

  output{
    File af_file = "~{callset_name}.frq"
    File hwe_file = "~{callset_name}.hwe"
    File af_file_euro = "~{callset_name}_euro.frq"
    File hwe_file_euro = "~{callset_name}_euro.hwe"
  }
}


task concatenateFiles {
  input {
    Array[File] files
    String callset_name
  }
  command <<<
    cat ~{sep=" " files} >> ~{callset_name}.tsv
  >>>
  runtime{
    docker:"ubuntu:latest"
    disks: "local-disk 50 HDD"
    preemptible: 3
  }
  output{
    File out = "~{callset_name}.tsv"
  }
}

