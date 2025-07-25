# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_1B_PERFORM_HWAS {
  input {
    File step_12_data = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/genitourinary_system/genitourinary_system.metadata"
    String cancer_type = "genitourinary_system"
 
    String analysis_1a_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/1A_HWAS/"
    String output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/1B_HWAS/" 
  }

  Int negative_shards = 0

  call Tasks.list_files_from_directory {
    input:
      dir = analysis_1a_dir,
      prefix = "aou",
      suffix = ".tsv.gz" 
  }

  scatter(i in range(length(list_files_from_directory.out1)-0)) {
    call T1_split_file {
      input:
        file_to_split = list_files_from_directory.out1[i],
        max_lines_per_chunk = 5000
    }
    scatter(i in range(length(T1_split_file.out1))){
      call T2_RunLogisticRegression {
        input:
          raw_roh_hwas = T1_split_file.out1[i],
          cancer_type = cancer_type,
          sample_data = step_12_data
      }
    }
    call Tasks.concatenateFiles_noheader as concat1 {
      input:
        files = T2_RunLogisticRegression.out1
    }
  }
  call Tasks.concatenateGzippedFiles_noheader as concat2 {
    input:
      files = concat1.out1,
      callset_name = cancer_type
  }
  call Tasks.copy_file_to_storage {
    input:
      text_file = concat2.out1,
      output_dir = output_dir 
  }
}

task T1_split_file {
  input {
    File file_to_split
    Int max_lines_per_chunk = 1000
  }

  command <<<
  set -euxo pipefail

    # Make output directory
    mkdir chunks

    # Unzip
    zcat ~{file_to_split} > file_to_split.txt
    rm ~{file_to_split}

    # Extract the header (first line)
    head -n 1 file_to_split.txt > header.tmp

    # Skip the header and split the rest
    tail -n +2 file_to_split.txt | \
      split -l ~{max_lines_per_chunk} - chunks/part_

    # Prepend the header to each split file
    for f in chunks/part_*; do
      cat header.tmp "$f" > "$f.with_header"
      mv "$f.with_header" "$f"
    done

    # Clean up
    rm header.tmp
  >>>

  output {
    Array[File] out1 = glob("chunks/part_*")
  }

  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }
}

task T2_RunLogisticRegression {
  input {
    File raw_roh_hwas
    String cancer_type
    File sample_data
  }

  command <<<
  set -euxo pipefail
 
  # Filter to samples of interest
  python3 <<CODE
  import os
  import pandas as pd
  import numpy as np
  import statsmodels.api as sm
  from scipy.stats import zscore

  # Read in dataframes
  haplotype_df = pd.read_csv("~{raw_roh_hwas}", sep='\t', index_col=False)
  df = pd.read_csv("~{sample_data}", sep='\t', index_col=False)
  df['is_case'] = df["cancer"].apply(lambda x: 0 if x == "control" else 1)
  df['sex_binary'] = df['inferred_sex'].apply(lambda x: 0 if x == "female" else 1)
  sample_list = [str(s).strip() for s in df['original_id']]

  # Clean haplotype_df column names
  haplotype_df.columns = haplotype_df.columns.astype(str).str.strip()

  # Subset haplotype_df to only include columns in sample_list
  haplotype_df = haplotype_df.loc[:, ['Haplotype'] + sample_list]

  # Covariates
  covariates = [f"PC{i}" for i in range(1, 5)] + ["age"]
  df[covariates] = df[covariates].apply(zscore)

  # Add 'sex_binary' if it varies
  if df['sex_binary'].nunique() > 1:
      covariates.append('sex_binary')
 
  # Drop rows with missing covariate values
  df = df.dropna(subset=covariates)

  # Prepare results
  results = []

  # Iterate over haplotypes
  for _, row in haplotype_df.iterrows():
      hap_id = row['Haplotype']
      values = row[sample_list].astype(float).values

      # Skip invariant haplotypes
      if np.all(values == values[0]):
          continue

      # Z-score normalize ROH values
      hap_values = (values - values.mean()) / values.std()

      # Design matrix
      X = df[covariates].copy()
      X['roh_score'] = hap_values
      X = sm.add_constant(X)

      y = df['is_case']

      if X.isnull().values.any():
          #print(f"NaN detected for haplotype {hap_id}")
          #print(X.isnull().sum())
          continue

      if np.isinf(X.values).any():
          #print(f"Inf detected for haplotype {hap_id}")
          continue
        

      model = sm.Logit(y, X)
      try:
          # Fit logistic regression
          #print(f"Running {hap_id}")
          result = model.fit(disp=0)
          beta = result.params['roh_score']
          pval = result.pvalues['roh_score']
      except:
          #print(f"Fail {hap_id}")
          beta, pval = np.nan, np.nan

      results.append({
          'haplotype': hap_id,
          'beta_roh_score': beta,
          'pvalue_roh_score': pval
      })

  # Save results
  results_df = pd.DataFrame(results)
  results_df.to_csv("haplotype_roh_logistic_results.tsv", sep="\t", index=False)

  CODE

  >>>
  output {
    File out1 = "haplotype_roh_logistic_results.tsv"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

