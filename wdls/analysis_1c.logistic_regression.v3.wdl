# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_1C_LOGISTIC_REGRESSION {
  input {
    Array[String] cancer_types = ["basal_cell","bladder","breast","cervix","colorectal","lung","uterus","ovary","kidney","squamous_cell","melanoma","prostate","brain","neuroendocrine","sarcoma","non-hodgkin","hematologic","thyroid"]
    String sw_size = 100000
    String roh_size = 0
    String output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/results/" 
  }
  File analysis_1b_output = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/1B_ALL_BY_ALL/analysis_1b_output.jan9.sw_" + sw_size + ".roh_" + roh_size + "kb.tsv.gz.tsv.gz"

  Int negative_shards = 0

  call T1_split_file {
    input:
      file_to_split = analysis_1b_output,
      max_lines_per_chunk = 5000
  }
  scatter(cancer_type in cancer_types){
    File step_12_data = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".metadata"
    scatter(i in range(length(T1_split_file.out1) - 0)){
      call T2_RunLogisticRegression {
        input:
          raw_roh_hwas = T1_split_file.out1[i],
          cancer_type = cancer_type,
          sample_data = step_12_data
      }
    }
    call concatenateFiles_noheader {
      input:
        files = T2_RunLogisticRegression.out1,
        callset_name = cancer_type + ".sw_" + sw_size + ".roh_" + roh_size,
        short_file_name = cancer_type + ".sw_" + sw_size + ".roh_" + roh_size + ".short"
    }
    call Tasks.copy_file_to_storage as copy1{
      input:
        text_file = concatenateFiles_noheader.out1,
        output_dir = output_dir 
    }
    call Tasks.copy_file_to_storage as copy2{
      input:
        text_file = concatenateFiles_noheader.out2,
        output_dir = output_dir
    }
  }
}

task concatenateFiles_noheader {
  input {
    Array[File] files
    String callset_name = "out"
    String short_file_name = "out"
  }

  command <<<
  # Put Header Down
  head -n 1 ~{files[0]} > ~{callset_name}.tsv
  head -n 1 ~{files[0]} > ~{short_file_name}.tsv

  # Add files
  for f in ~{sep=" " files}; do
    tail -n +2 "$f"
  done >> tmp.tsv
  sort -k4,4g tmp.tsv >>  ~{callset_name}.tsv
  gzip -c ~{callset_name}.tsv > ~{callset_name}.tsv.gz

  sort -k4,4g tmp.tsv | head -n 50 >>  ~{short_file_name}.tsv
  gzip -c ~{short_file_name}.tsv > ~{short_file_name}.tsv.gz
  >>>

  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }

  output {
    File out1 = "~{callset_name}.tsv.gz"
    File out2 = "~{short_file_name}.tsv.gz"
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
  df = df[df['intake_qc_pop'] == "EUR"]

  # Ensure original_id is string and stripped
  df["original_id"] = df["original_id"].astype(str).str.strip()

  # Phenotype / metadata columns
  df["is_case"] = (df["cancer"] != "control").astype(int)
  df["sex_binary"] = (df["inferred_sex"] != "female").astype(int)

  # ------------------------------------------------------------
  # Prepare haplotype matrix (DEDUPLICATE GENES PROPERLY)
  # ------------------------------------------------------------

  # Drop genomic region
  hap_df = haplotype_df.drop(columns=["Genomic_Region"])

  # Clean GeneName
  hap_df["GeneName"] = hap_df["GeneName"].astype(str).str.strip()

  # Convert all sample columns to numeric (everything except GeneName)
  sample_cols = hap_df.columns.difference(["GeneName"])
  hap_df[sample_cols] = hap_df[sample_cols].apply(pd.to_numeric, errors="coerce")

  # Collapse multiple windows per gene (mean ROH per gene)
  hap_df = (
      hap_df
      .groupby("GeneName", as_index=False)
      .mean()
  )

  # Transpose: samples as rows, genes as columns
  hap_matrix = hap_df.set_index("GeneName").T
  hap_matrix.index.name = "original_id"

  # Reset index so original_id is a column
  hap_matrix = hap_matrix.reset_index()
  hap_matrix.columns.name = None

  # Identify gene columns ONCE
  gene_cols = hap_matrix.columns.difference(["original_id"])

  # Ensure gene columns are numeric (safe, no fragmentation)
  hap_matrix[gene_cols] = hap_matrix[gene_cols].astype(float)

  # Strip whitespace from IDs
  hap_matrix["original_id"] = hap_matrix["original_id"].astype(str).str.strip()

  # ------------------------------------------------------------
  # Merge metadata with haplotypes
  # ------------------------------------------------------------
  merged_df = df.merge(hap_matrix, on="original_id", how="inner")

  # ------------------------------------------------------------
  # Covariates
  # ------------------------------------------------------------
  pc_covariates = [f"PC{i}" for i in range(1, 5)]
  merged_df[pc_covariates] = merged_df[pc_covariates].apply(zscore)

  covariates = pc_covariates
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
      X["roh_score"] = hap_values.values #(np.exp(5 * np.clip(hap_values.values, 0, 1)) - 1) / (np.exp(5) - 1)
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

