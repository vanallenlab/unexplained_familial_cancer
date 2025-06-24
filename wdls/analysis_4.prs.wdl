# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_4_PRS {
  input {
    String step_8_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_8_FILTER_TO_TP_VARIANTS/sharded_vcfs"
    String PGS_ID
    String cancer_type
    String analysis_4_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/"
    String workspace_bucket = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
  }
  File sample_data = "gs://~{workspace_bucket}/UFC_REFERENCE_FILES/analysis/~{cancer_type}/~{cancer_type}.metadata"
  File prs_file = "gs://~{workspace_bucket}/UFC_REFERENCE_FILES/pgs_files/~{PGS_ID}_hmPOS_GRCh38.txt.gz" 
  Int negative_shards = 0

  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call Tasks.gather_vcfs as gather_vcfs{
    input:
      dir = step_8_output_dir
  }

  call Tasks.sort_vcf_list {
    input:
      unsorted_vcf_list = gather_vcfs.vcf_list
  }

  call T1_prepare_scores_file {
    input:
      prs_file = prs_file
  }

  scatter (i in range(length(sort_vcf_list.vcf_arr) - negative_shards)){
    call T2_calculate_scores {
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        cancer_type = cancer_type,
        score_file = T1_prepare_scores_file.out1,
        PGS_ID = PGS_ID,
        shard = i
    }
  }

  call T3_sum_scores {
    input:
      pgs_scores = T2_calculate_scores.out1,
      cancer_type = cancer_type,
      PGS_ID = PGS_ID
  }

  call T4_control_for_ancestry {
    input:
      raw_prs = T3_sum_scores.out1,
      sample_data = sample_data,
      cancer_type = cancer_type,
      PGS_ID = PGS_ID
  }
  call T5_get_summary_statistics {
    input:
      adjusted_prs = T4_control_for_ancestry.out1,
      sample_data = sample_data,
      cancer_type = cancer_type,
      PGS_ID = PGS_ID
  }

  call Tasks.copy_file_to_storage as copy0{
   input:
      text_file = T3_sum_scores.out1,
      output_dir = analysis_4_output_dir
  }

  call Tasks.copy_file_to_storage as copy1{
    input:
      text_file = T4_control_for_ancestry.out1,
      output_dir = analysis_4_output_dir
  }
  call Tasks.copy_file_to_storage as copy2{
    input:
      text_file = T5_get_summary_statistics.out1,
      output_dir = analysis_4_output_dir
  } 
}

task T1_prepare_scores_file {
  input{
    File prs_file
  }
  command <<<
  python3 <<CODE
  import pandas as pd

  usecols = ['hm_source','hm_chr', 'hm_pos', 'effect_allele', 'other_allele', 'effect_weight']

  df = pd.read_csv("~{prs_file}", comment = "#", index_col=False, sep='\t',dtype={'hm_chr': str, 'hm_pos': str},usecols=usecols)
  df = df[df['hm_source'] == "liftover"]

  # Create primary and flipped IDs
  df1 = df[['hm_chr', 'hm_pos', 'effect_allele', 'other_allele', 'effect_weight']].copy()
  df1['ID'] = "chr" + df1['hm_chr'] + "_" + df1['hm_pos'] + "_" + df1['effect_allele'] + "_" + df1['other_allele']

  df2 = df[['hm_chr', 'hm_pos', 'effect_allele', 'other_allele', 'effect_weight']].copy()
  df2['ID'] = "chr" + df2['hm_chr'] + "_" + df2['hm_pos'] + "_" + df2['other_allele'] + "_" + df2['effect_allele']

  # Combine both versions
  score_df = pd.concat([
      df1[['ID', 'effect_allele', 'effect_weight']],
      df2[['ID', 'effect_allele', 'effect_weight']]
  ])

  # Drop duplicates just in case
  score_df = score_df.drop_duplicates()

  # Save to file
  score_df.to_csv("score_file.txt", sep="\t", index=False, header=False)
  CODE
  >>>
  output {
    File out1 = "score_file.txt"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
    memory: "8GB"
  }
}

task T2_calculate_scores {
  input{
    File vcf
    File score_file
    String cancer_type
    String PGS_ID
    String shard
  }
  command <<<
  # Calculate score
  plink --vcf ~{vcf} --score ~{score_file} sum --double-id --out tmp > plink.log || true

  if [ -f tmp.profile ]; then
    # Reformat
    awk '{print $2 "\t" $6}' tmp.profile > "~{cancer_type}.~{PGS_ID}.shard_~{shard}.profile"
  else
    echo -e "IID\tSCORE" > "~{cancer_type}.~{PGS_ID}.shard_~{shard}.profile"
  fi

  >>>
  output {
    File out1 = "~{cancer_type}.~{PGS_ID}.shard_~{shard}.profile"
  }
  runtime {
    docker: "elixircloud/plink:1.9-20210614"
    preemptible: 3
  }
}

task T3_sum_scores {
  input {
    Array[File] pgs_scores
    String cancer_type
    String PGS_ID
  }
  command <<<
  python3 <<CODE
  # This will hold cumulative scores
  total_scores = {}

  # Replace with your actual list of file paths or use glob
  pgs_files = "~{sep=' ' pgs_scores}".split()

  for filename in pgs_files:
      with open(filename) as f:
          next(f)  # skip header
          for line in f:
              sample, score_str = line.strip().split('\t')
              score = float(score_str)
              total_scores[sample] = total_scores.get(sample, 0.0) + score

  # Write final output
  with open("~{cancer_type}.~{PGS_ID}.raw.pgs", "w") as out:
      out.write("sample\tPGS\n")
      for sample, score in sorted(total_scores.items()):
          out.write(f"{sample}\t{round(score, 10)}\n")
  CODE
  >>>
  output {
    File out1 = "~{cancer_type}.~{PGS_ID}.raw.pgs"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task T4_control_for_ancestry {
  input {
    File raw_prs
    File sample_data
    String cancer_type
    String PGS_ID 
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  import statsmodels.api as sm
  from scipy.stats import zscore

  # Load input files
  prs_df = pd.read_csv("~{raw_prs}", sep='\t', index_col=False)
  sample_data = pd.read_csv("~{sample_data}", sep='\t', index_col=False)

  # Convert merge keys to string
  prs_df['sample'] = prs_df['sample'].astype(str)
  sample_data['original_id'] = sample_data['original_id'].astype(str)

  # Merge on sample ID
  merged_df = prs_df.merge(sample_data, left_on = "sample", right_on="original_id")

  merged_df[['PC1','PC2','PC3','PC4','age']] = merged_df[['PC1','PC2','PC3','PC4','age']].apply(zscore)
  merged_df['male_sex'] = merged_df['inferred_sex'].apply(lambda x: 1 if x == 'male' else 0)

  covariates = ['PC1', 'PC2', 'PC3', 'PC4']
  if merged_df['male_sex'].nunique() > 1:
    covariates.append('male_sex')

  # Set up regression: regress PGS ~ PC1 + PC2 + PC3 + PC4 + male_sex
  X = merged_df[covariates]
  X = sm.add_constant(X)
  y = merged_df['PGS']

  model = sm.OLS(y, X).fit()
  merged_df['PGS'] = model.resid

  # Calculate mean and std using only controls
  control_pgs = merged_df.loc[merged_df['cancer'] == 'control', 'PGS']
  control_mean = control_pgs.mean()
  control_std = control_pgs.std()

  # Standardize PGS for all samples using control-derived mean and std
  merged_df['PGS'] = (merged_df['PGS'] - control_mean) / control_std

  # Write final output
  merged_df[['sample', 'PGS']].to_csv("~{cancer_type}.~{PGS_ID}.pgs", sep='\t', index=False)

  CODE
  >>>
  output {
    File out1 = "~{cancer_type}.~{PGS_ID}.pgs"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
    memory: "8GB"
  }
}


task T5_get_summary_statistics {
  input {
    File adjusted_prs
    File sample_data
    String cancer_type
    String PGS_ID
  }
  command <<<
  python3 <<CODE
  from scipy.stats import ttest_ind
  from sklearn.metrics import roc_auc_score
  import statsmodels.api as sm
  import pandas as pd
  import numpy as np
  from scipy.stats import zscore

  f = open("~{cancer_type}.~{PGS_ID}.stats","w")

  df = pd.read_csv("~{adjusted_prs}",sep='\t',index_col=False)
  df['pgs_score'] = df['PGS']
  sample_data = pd.read_csv("~{sample_data}",sep='\t',index_col=False)
  df['sample'] = df['sample'].astype(str)
  sample_data['original_id'] = sample_data['original_id'].astype(str)
  merged_df = df.merge(sample_data, left_on = "sample", right_on = "original_id", how = "left")

  # 0. Normalize in PCs
  merged_df[['PC1','PC2','PC3','PC4','age']] = merged_df[['PC1','PC2','PC3','PC4','age']].apply(zscore)

  # 1. T-test between case and control PGS scores
  cases = merged_df[merged_df['cancer'] != 'control']['pgs_score']
  controls = merged_df[merged_df['cancer'] == 'control']['pgs_score']
  t_stat, p_val = ttest_ind(cases, controls, equal_var=False)
  f.write(f"T-test p-value (case vs control): {p_val:.4e}\n")

  # 2. Mean and std of the two distributions
  f.write(f"Case mean: {cases.mean():.4f}, std: {cases.std():.4f}\n")
  f.write(f"Control mean: {controls.mean():.4f}, std: {controls.std():.4f}\n")

  # 3. AUC of PRS
  # Binary labels: 1 = case, 0 = control
  y_true = merged_df['cancer'].apply(lambda x: 0 if x == 'control' else 1)
  auc = roc_auc_score(y_true, merged_df['pgs_score'])
  f.write(f"AUC of PGS score: {auc:.4f}\n")

  # 4. Logistic regression of cancer status on PGS with covariates
  # Prepare design matrix
  merged_df['sex_binary'] = merged_df['inferred_sex'].map({'male': 1, 'm': 1, 'female': 0, 'f': 0})

  # Base covariates
  covariates = ['pgs_score', 'PC1', 'PC2', 'PC3', 'PC4', 'age']

  # Add 'sex_binary' if it varies
  if merged_df['sex_binary'].nunique() > 1:
      covariates.append('sex_binary')

  Xy = merged_df[covariates + ['cancer']].dropna()

  # Now safely split cancer out to y
  X = Xy[covariates]
  y = Xy['cancer'].apply(lambda x: 0 if x == 'control' else 1)

  # Add constant for intercept
  X = sm.add_constant(X)

  # Fit logistic regression model
  logit_model = sm.Logit(y, X)
  result = logit_model.fit()

  # Print summary of coefficients
  f.write(result.summary2().as_text())
  f.close()
  CODE
  >>>
  output {
    File out1 = "~{cancer_type}.~{PGS_ID}.stats"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
    memory: "8GB"
  }
}
