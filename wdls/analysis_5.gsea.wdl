# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_5_GSEA {
  input {
    String step_10_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_10_VISUALIZE_VEP"
    Array[String] cancer_types = ["basal_cell","bladder","bone","soft_tissue","uterus","endometrial","ovary","breast","colorectal","kidney","hematologic","lymphoma","nervous","prostate","squamous_cell","melanoma","thyroid","respiratory"]
    String analysis_5_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_5_GSEA/"
    String workspace_bucket = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
    String biological_pathway = "tsg_dominant"
  }
  
  scatter( cancer_type in cancer_types){
    File sample_data = "gs://" + workspace_bucket + "/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".metadata"
    #File gene_list = "gs://" + workspace_bucket + "/UFC_REFERENCE_FILES/gene_list/" + biological_pathway + ".gene_list"
    File gene_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/all_tsg.list"
    Int negative_shards = 0

    # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
    call Tasks.list_files_from_directory {
      input:
        dir = step_10_output_dir,
        suffix = ".tsv.gz" 
    }

    scatter (i in range(length(list_files_from_directory.out1) - negative_shards)){
      call T1_get_rows {
        input:
          gene_list = gene_list,
          variant_tsv = list_files_from_directory.out1[i]
      }
    }
    call Tasks.concatenateFiles_noheader as concat1a {
      input:
        files = T1_get_rows.out1,
        callset_name = biological_pathway
    }
    call Tasks.concatenateFiles_noheader as concat1b {
      input:
        files = T1_get_rows.out2,
        callset_name = biological_pathway
    }
    call T2_gsea as T2_gsea_revel050{
      input:
        variants_tsv = concat1a.out1,
        sample_metadata = sample_data,
        cancer_type = cancer_type,
        path_threshold = "REVEL_050",
        allele_frequency = "1%"
    }
    call T2_gsea as T2_gsea_revel075{
      input:
        variants_tsv = concat1b.out1,
        sample_metadata = sample_data,
        cancer_type = cancer_type,
        path_threshold = "REVEL_075",
        allele_frequency = "1%"
    }
  }
  call Tasks.concatenateFiles_noheader as concat2a{
    input:
      files = T2_gsea_revel050.out1,
      callset_name = biological_pathway + "_" + "001"
  }
  call Tasks.concatenateFiles_noheader as concat2b{
    input:
      files = T2_gsea_revel075.out1,
      callset_name = biological_pathway + "_" + "001"
  }
  call Tasks.concatenateFiles_noheader as concat3{
    input:
      files = [concat2a.out2,concat2b.out2],
      callset_name = biological_pathway + "_" + "001"
  }
}


task T1_get_rows {
  input {
    File variant_tsv
    File gene_list
  }

  command <<<
  set -euxo pipefail
  python3<<CODE
  import pandas as pd
  df = pd.read_csv("~{variant_tsv}",sep='\t',index_col=False)
  df = df.drop_duplicates()
  # Read gene list from file
  with open("~{gene_list}") as f:
      genes = [line.strip() for line in f if line.strip()]

  # Filter where gene_impact is in genes
  df_050 = df[df['gene_impact'].isin(
      [f"{g}_{suffix}" for g in genes for suffix in ["REVEL050", "REVEL075"]]
  )]
  df_050.to_csv("revel050.tsv",sep='\t',index=False)

  # Filter where gene_impact is in genes
  df_075 = df[df['gene_impact'].isin(
      [f"{g}_{suffix}" for g in genes for suffix in ["REVEL075"]]
  )]
  df_050.to_csv("revel075.tsv",sep='\t',index=False) 
  CODE
  >>>
  output{
    File out1 = "revel050.tsv"
    File out2 = "revel075.tsv"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    disk: "local-disk 10 HDD"
    memory: "2GB"
    preemptible: 3
    cpu: 2
  }
}

task T2_gsea {
  input {
    File variants_tsv
    File sample_metadata
    String path_threshold
    String cancer_type
    String allele_frequency
  }

  command <<<
  set -euo pipefail

  python3 <<CODE
  import pandas as pd
  from scipy.stats import zscore
  import statsmodels.api as sm
  import numpy as np

  # Load the filtered data
  df = pd.read_csv("~{variants_tsv}", sep="\t",index_col=False)

  # Drop potential duplicates before summing
  df = df.drop_duplicates(subset=["gene_impact"])
  df = df.replace(2, 0) 
  # -----------------
  # 2. Sum across rows to get ALL_HIGH
  # -----------------
  patient_cols = df.columns.drop("gene_impact")
  all_high = df[patient_cols].sum(axis=0)  # sums down each patient column

  # Create a DataFrame with patient IDs and counts
  all_high_df = all_high.reset_index()
  all_high_df.columns = ["original_id", "num_pathogenic_variants"]

  # -----------------
  # 3. Merge with metadata
  # -----------------
  meta = pd.read_csv("~{sample_metadata}", sep="\t",index_col=False)
  meta['original_id'] = meta['original_id'].astype(str)
  all_high_df['original_id'] = all_high_df['original_id'].astype(str)
  merged = pd.merge(meta, all_high_df, on="original_id", how="left")

  # Create binary 'is_case' column: 1 if cancer, 0 if control
  merged['is_case'] = merged['cancer'].apply(lambda x: 0 if x == 'control' else 1)

  # Create binary 'sex_binary': 0 = female, 1 = male (or anything else)
  merged['sex_binary'] = merged['inferred_sex'].apply(lambda x: 0 if x == 'female' else 1)
  merged_df['smoking_history'] = (
      merged_df['smoking_history']
      .replace("NA", 0)        # turn 'NA' into 0
      .astype(float)           # force numeric dtype
  )

  # Define covariates to zscore
  covariates = [f"PC{i}" for i in range(1, 5)] + ['age',"num_pathogenic_variants"]
  merged[covariates] = merged[covariates].apply(zscore)

  # Z-score these covariates (handle missing data safely)
  merged[covariates] = merged[covariates].apply(lambda col: zscore(col.fillna(col.mean())))

  # Add 'sex_binary' to covariates only if it varies
  if merged['sex_binary'].nunique() > 1:
      covariates.append('sex_binary')
  # Add 'sex_binary' to covariates only if it varies
  if merged['smoking_history'].nunique() > 1:
      covariates.append('smoking_history')

  # Covariates
  X = merged[covariates]

  # Add intercept
  X = sm.add_constant(X)
  y = merged["is_case"]

  # Fit logistic regression
  model = sm.Logit(y, X)
  result = model.fit()

  # Extract stats
  coef = result.params[1]              # beta for predictor (assuming 1 predictor after intercept)
  se = result.bse[1]                   # standard error
  ci_low, ci_high = result.conf_int().iloc[1]  # 95% CI
  pval = result.pvalues[1]             # p-value

  # Write one-line output
  with open("logistic_regression_summary.txt", "w") as f:
      f.write("cancer\tPathogenic Threshold\tAF\tp_value\tbeta\t95%_CI_lower\t95%_CI_upper\n")
      f.write(f"~{cancer_type}\t~{path_threshold}\t~{allele_frequency}\t{pval:.3e}\t{coef:.4f}\t{ci_low:.4f}\t{ci_high:.4f}\n")

  CODE
  >>>

  output {
    File out1 = "logistic_regression_summary.txt"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "4G"
    disk: "local-disks 10 HDD"
  }
}
