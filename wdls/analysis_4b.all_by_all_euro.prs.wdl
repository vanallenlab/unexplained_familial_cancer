# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_4B_PRS {
  input {
    File raw_prs_basename = "PGS000783.raw.pgs" 
    String PGS_ID = "PGS000783"
    #Array[String] cancer_types = ["basal_cell","bladder"]
    Array[String] cancer_types = ["basal_cell","bladder","breast","colorectal","gastrointestinal","hematologic","kidney","lung","melanoma","neuroendocrine","nervous","non-hodgkins","ovary","prostate","sarcoma","squamous_cell","thyroid","uterus","viral"]
    String analysis_4_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/"
    String workspace_bucket = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
  }
  File sample_data = "gs://~{workspace_bucket}/UFC_REFERENCE_FILES/analysis/"
  Int negative_shards = 0

  String raw_prs_file = analysis_4_output_dir + raw_prs_basename

  scatter( cancer_type in cancer_types){
    File specific_cohort_sample_data = sample_data + cancer_type + "/" + cancer_type + ".metadata"
    call T4_control_for_ancestry {
      input:
        raw_prs = raw_prs_file,
        sample_data = specific_cohort_sample_data,
        cancer_type = cancer_type,
        PGS_ID = PGS_ID
    }

    call Tasks.copy_file_to_storage as copy1{
      input:
        text_file = T4_control_for_ancestry.out1,
        output_dir = analysis_4_output_dir + "ADJUSTED_PRS/"
    }

    call T4_control_for_ancestry as T4_control_for_ancestry_euro {
      input:
        raw_prs = raw_prs_file,
        sample_data = specific_cohort_sample_data,
        cancer_type = cancer_type,
        PGS_ID = PGS_ID,
        euro = "True"
    }

    call Tasks.copy_file_to_storage as copy2{
      input:
        text_file = T4_control_for_ancestry_euro.out1,
        output_dir = analysis_4_output_dir + "ADJUSTED_PRS/"
    }

    call T5_get_summary_statistics {
      input:
        adjusted_prs = T4_control_for_ancestry.out1,
        sample_data = specific_cohort_sample_data,
        cancer_type = cancer_type,
        PGS_ID = PGS_ID
    }
    call T5_get_summary_statistics as T5_get_summary_statistics_euro{
      input:
        adjusted_prs = T4_control_for_ancestry_euro.out1,
        sample_data = specific_cohort_sample_data,
        cancer_type = cancer_type,
        PGS_ID = PGS_ID
    }
  }

  call Tasks.concatenateFiles {
    input:
      files = T5_get_summary_statistics.out1,
      output_name = PGS_ID
  }

  call Tasks.concatenateFiles as concatenateFiles_euro{
    input:
      files = T5_get_summary_statistics_euro.out1,
      output_name = PGS_ID + "_euro"
  }

  call Tasks.copy_file_to_storage as copy3{
    input:
      text_file = concatenateFiles.out1,
      output_dir = analysis_4_output_dir + "ALL_BY_ALL/"
  }

  call Tasks.copy_file_to_storage as copy4{
    input:
      text_file = concatenateFiles_euro.out1,
      output_dir = analysis_4_output_dir + "ALL_BY_ALL/"
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
  df = df[df['hm_source'].notna() & (df['hm_source'].str.strip() != "")]

  # Keep only rows where hm_chr is in 1-22, X, or Y
  valid_chroms = set([str(i) for i in range(1, 23)] + ['X', 'Y'])
  df = df[df['hm_chr'].isin(valid_chroms)]

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

task T4_control_for_ancestry {
  input {
    File raw_prs
    File sample_data
    String cancer_type
    String PGS_ID
    String euro = "False"
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  import statsmodels.api as sm
  from scipy.stats import zscore

  # Load input files
  prs_df = pd.read_csv("~{raw_prs}", sep='\t', index_col=False)
  sample_data = pd.read_csv("~{sample_data}", sep='\t', index_col=False)

  print("1")
  # Limit to just Europeans for Sensitivity Analysis
  euro = "~{euro}" == "True"
  if euro:
    sample_data = sample_data[sample_data['intake_qc_pop'] == "EUR"]

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
  print("2")
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

  # Choose filename
  suffix = ".euro" if euro else ""
  outname = f"~{cancer_type}.~{PGS_ID}{suffix}.pgs"

  # Write final output
  merged_df[['sample', 'PGS']].to_csv(outname, sep='\t', index=False)

  CODE
  >>>
  output {
    File out1 = glob("~{cancer_type}.~{PGS_ID}*.pgs")[0]
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

  f.write("\n\n####~{cancer_type}####\n\n".capitalize())
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
  merged_df['smoking_history'] = (
      merged_df['smoking_history']
      .replace("NA", 0)        # turn 'NA' into 0
      .astype(float)           # force numeric dtype
  )

  # Base covariates
  covariates = ['pgs_score', 'PC1', 'PC2', 'PC3', 'PC4', 'age', 'smoking_history']

  # Add 'sex_binary' if it varies
  if merged_df['sex_binary'].nunique() > 1:
      covariates.append('sex_binary')

  Xy = merged_df[covariates + ['cancer']].dropna()

  # Now safely split cancer out to y
  y = Xy['cancer'].apply(lambda x: 0 if x == 'control' else 1)

  # Check smoking history distribution by group
  if 'smoking_history' in covariates:
      tab = pd.crosstab(y, Xy['smoking_history'], dropna=False)

      # Case group (y=1) and control group (y=0)
      cases = Xy.loc[y == 1, 'smoking_history']
      controls = Xy.loc[y == 0, 'smoking_history']

      # Drop if all cases or all controls are the same (0 or 1 only)
      if cases.nunique(dropna=True) <= 1 or controls.nunique(dropna=True) <= 1:
          covariates.remove('smoking_history')

  X = Xy[covariates]

  # Add constant for intercept
  X = sm.add_constant(X)

  # Fit logistic regression model
  logit_model = sm.Logit(y, X)
  result = logit_model.fit()

  # Get summary table as a DataFrame
  summary_df = result.summary2().tables[1]

  # Format p-values in scientific notation with higher precision
  summary_df['P>|z|'] = summary_df['P>|z|'].apply(lambda p: f"{p:.5e}" if not np.isnan(p) else "nan")

  # Write formatted table to file
  f.write(summary_df.to_string())

  # Print summary of coefficients
  #f.write(result.summary2().as_text())
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
