# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_1C_LOGISTIC_REGRESSION {
  input {
    String cancer_type = "thyroid"
 
    String analysis_1a_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/1B_GENES_0kb/"
    String output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/1C_RESULTS_GENES_0kb/" 
  }
  File prs_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/ADJUSTED_PRS/thyroid.PGS000797.pgs"
  File ppv_missense_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_5_GSEA/thyroid.num_ppv_missense.tsv"
  File step_12_data = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".metadata"
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
          sample_data = step_12_data,
          prs_file = prs_file,
          ppv_missense_file = ppv_missense_file
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
      callset_name = cancer_type + "_roh"
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
    File prs_file
    File ppv_missense_file
  }

  command <<<
  set -euxo pipefail
  python3 <<CODE
  import pandas as pd
  import numpy as np
  import statsmodels.api as sm
  from scipy.stats import zscore

  # ------------------------------------------------------------
  # Load data
  # ------------------------------------------------------------
  haplotype_df = pd.read_csv("~{raw_roh_hwas}", sep="\t", index_col=False)
  df = pd.read_csv("~{sample_data}", sep="\t")

  # Ensure original_id is string and stripped
  df["original_id"] = df["original_id"].astype(str).str.strip()

  # Phenotype / metadata columns
  df["is_case"] = (df["cancer"] != "control").astype(int)
  df["sex_binary"] = (df["inferred_sex"] != "female").astype(int)

  # ------------------------------------------------------------
  # Add PRS and PPV data
  # ------------------------------------------------------------
  prs_df = pd.read_csv("~{prs_file}", sep="\t")
  prs_df["original_id"] = prs_df["sample"].astype(str).str.strip()

  ppv_df = pd.read_csv("~{ppv_missense_file}", sep="\t")
  ppv_df["original_id"] = ppv_df["original_id"].astype(str).str.strip()

  # Merge with metadata, keep original_id as a column
  df = df.merge(prs_df[['original_id','PGS']], on="original_id", how="left")
  df = df.merge(ppv_df[['original_id','num_ppv_missense']], on="original_id", how="left")

  # ------------------------------------------------------------
  # Prepare haplotype matrix
  # ------------------------------------------------------------
  hap_matrix = haplotype_df.drop(columns=['Genomic_Region'])
  hap_matrix = hap_matrix.set_index('GeneName').T  # patients as rows, genes as columns
  hap_matrix.index.name = 'original_id'
  hap_matrix = hap_matrix.reset_index()  # make original_id a column
  hap_matrix.columns.name = None

  # Convert all gene columns to numeric
  gene_cols = [c for c in hap_matrix.columns if c != 'original_id']
  hap_matrix[gene_cols] = hap_matrix[gene_cols].apply(pd.to_numeric, errors='coerce')

  # Strip whitespace from original_id
  hap_matrix['original_id'] = hap_matrix['original_id'].astype(str).str.strip()

  # ------------------------------------------------------------
  # Merge metadata with haplotypes
  # ------------------------------------------------------------
  merged_df = df.merge(hap_matrix, on="original_id", how="inner")

  # ------------------------------------------------------------
  # Covariates
  # ------------------------------------------------------------
  pc_covariates = [f"PC{i}" for i in range(1, 5)]
  merged_df[pc_covariates] = merged_df[pc_covariates].apply(zscore)

  covariates = pc_covariates + ["PGS"]
  if merged_df["sex_binary"].nunique() > 1:
      covariates.append("sex_binary")

  # Drop rows with missing covariates
  merged_df = merged_df.dropna(subset=covariates)
  # ------------------------------------------------------------
  # Logistic regression per haplotype (gene)
  # ------------------------------------------------------------
  results = []

  for gene in gene_cols:
      hap_values = merged_df[gene].astype(float)

      # Skip invariant haplotypes
      if hap_values.nunique() <= 1:
          continue

      # Build design matrix
      X = merged_df[covariates].copy()
      X["roh_score"] = hap_values.values
      X = sm.add_constant(X)
      y = merged_df["is_case"]

      # Skip problematic matrices
      if X.isnull().any().any() or np.isinf(X.values).any():
          continue

      try:
          model = sm.Logit(y, X)
          result = model.fit(disp=0)

          results.append({
              "haplotype": gene,
              "beta_roh_score": result.params["roh_score"],
              "se_beta_roh_score": result.bse["roh_score"],
              "pvalue_roh_score": result.pvalues["roh_score"],
              "mean_case_value": hap_values[y==1].mean(),
              "mean_control_value": hap_values[y==0].mean()
          })
      except Exception:
          print(gene)
          continue

  # ------------------------------------------------------------
  # Save results
  # ------------------------------------------------------------
  results_df = pd.DataFrame(results)
  results_df.to_csv("haplotype_roh_logistic_results.tsv", sep="\t", index=False)
  print(f"Saved {len(results)} haplotype logistic regression results.")

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

task T2_RunLogisticRegression_old {
  input {
    File raw_roh_hwas
    String cancer_type
    File sample_data
    File prs_file
    File ppv_missense_file
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

  # --- Load data ---
  haplotype_df = pd.read_csv("~{raw_roh_hwas}", sep='\t', index_col=False)

  df = pd.read_csv("~{sample_data}", sep='\t', index_col=False)
  df['original_id'] = df['original_id'].astype(str).str.strip()   # <-- REQUIRED

  # Sample ordering based on IDs, not dataframe index
  sample_order = [s for s in df['original_id'] if s in haplotype_df.columns]

  df['is_case'] = df["cancer"].apply(lambda x: 0 if x == "control" else 1)
  df['sex_binary'] = df['inferred_sex'].apply(lambda x: 0 if x == "female" else 1)

  prs_df = pd.read_csv("~{prs_file}", sep='\t', index_col=False)
  prs_df['original_id'] = prs_df['sample'].astype(str).str.strip()

  ppv_df = pd.read_csv("~{ppv_missense_file}", sep='\t', index_col=False)
  ppv_df['original_id'] = ppv_df['original_id'].astype(str).str.strip()

  # Merges now work because original_id is a clean column
  df = df.merge(prs_df, on="original_id")
  df = df.merge(ppv_df, on="original_id", how="left")

  # Only now set the index
  df = df.set_index("original_id")

  # Clean haplotype_df column names
  haplotype_df.columns = haplotype_df.columns.astype(str).str.strip()

  # Subset haplotype_df to only include columns in sample_list
  haplotype_df = haplotype_df.loc[:, ['Genomic_Region'] + sample_order]

  # Covariates
  covariates = [f"PC{i}" for i in range(1, 5)] #+ ["age"]
  covariates = covariates + ['PGS','num_ppv_missense']
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
      hap_id = row['Genomic_Region']
      hap_values = row[sample_order].astype(float)

      # Skip invariant haplotypes
      if np.all(hap_values == hap_values[0]):
          continue

      # Design matrix
      X = df.loc[sample_order, covariates].copy()
      X['roh_score'] = hap_values.values
      X = sm.add_constant(X)

      y = df.loc[sample_order, 'is_case']

      if X.isnull().values.any():
          continue

      if np.isinf(X.values).any():
          continue
        

      model = sm.Logit(y, X)
      try:
          # Fit logistic regression
          result = model.fit(disp=0)
          beta = result.params['roh_score']
          se_beta = result.bse['roh_score']
          pval = result.pvalues['roh_score']

          # Mean haplotype values for cases and controls
          mean_case = hap_values[df['is_case'] == 1].mean()
          mean_ctrl = hap_values[df['is_case'] == 0].mean()

          results.append({
              'haplotype': hap_id,
              'beta_roh_score': beta,
              'se_beta_roh_score': se_beta,
              'pvalue_roh_score': pval,
              'mean_case_value': mean_case,
              'mean_control_value': mean_ctrl
          })
      except:
          continue

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

