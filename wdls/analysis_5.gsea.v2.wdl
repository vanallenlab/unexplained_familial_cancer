# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_5_GSEA {
  input {
    String step_10_cpg_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_10_VISUALIZE_VEP/ufc_cpg.tsv"
    Array[String] cancer_types = ["basal_cell","bladder","breast","colorectal","hematologic","kidney","lung","melanoma","neuroendocrine","nervous","non-hodgkins","ovary","prostate","sarcoma","squamous_cell","thyroid","uterus","cervix"]
    #Array[String] cancer_types = ["basal_cell","bladder"]
    String analysis_5_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_5_GSEA/"
    String workspace_bucket = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
    String biological_pathway = "all_cpg"
  }
  
  scatter( cancer_type in cancer_types){
    call T0_get_path_genes {
      input:
        cosmic_tsv = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cosmic_ufc.tsv",
        cancer_type = cancer_type
    }
    File sample_data = "gs://" + workspace_bucket + "/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".metadata"
    File gene_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/all_cpg_genes.list"
    Int negative_shards = 0


    call T1_get_rows {
      input:
        gene_list = T0_get_path_genes.out1,
        variant_tsv = step_10_cpg_file
    }
    call T2_gsea as T2_gsea_TIER0_001{
      input:
        variants_tsv = T1_get_rows.out7,
        sample_metadata = sample_data,
        cancer_type = cancer_type,
        path_threshold = "Tier0_001",
        allele_frequency = "0.01"
    }
    call T2_gsea as T2_gsea_TIER0_0001{
      input:
        variants_tsv = T1_get_rows.out8,
        sample_metadata = sample_data,
        cancer_type = cancer_type,
        path_threshold = "Tier0_0001",
        allele_frequency = "0.001"
    }
    call T2_gsea as T2_gsea_TIER2_001{
      input:
        variants_tsv = T1_get_rows.out5,
        sample_metadata = sample_data,
        cancer_type = cancer_type,
        path_threshold = "Tier2_001",
        allele_frequency = "0.01"
    }
    call T2_gsea as T2_gsea_TIER2_0001{
      input:
        variants_tsv = T1_get_rows.out6,
        sample_metadata = sample_data,
        cancer_type = cancer_type,
        path_threshold = "Tier2_0001",
        allele_frequency = "0.001"
    }
    call T2_gsea as T2_gsea_TIER3_001{
      input:
        variants_tsv = T1_get_rows.out1,
        sample_metadata = sample_data,
        cancer_type = cancer_type,
        path_threshold = "Tier3_001",
        allele_frequency = "0.01"
    }
    call T2_gsea as T2_gsea_TIER3_0001{
      input:
        variants_tsv = T1_get_rows.out2,
        sample_metadata = sample_data,
        cancer_type = cancer_type,
        path_threshold = "Tier3_0001",
        allele_frequency = "0.001"
    }
    call T2_gsea as T2_gsea_TIER4_001{
      input:
        variants_tsv = T1_get_rows.out3,
        sample_metadata = sample_data,
        cancer_type = cancer_type,
        path_threshold = "Tier4_001",
        allele_frequency = "0.01"
    }
    call T2_gsea as T2_gsea_TIER4_0001{
      input:
        variants_tsv = T1_get_rows.out4,
        sample_metadata = sample_data,
        cancer_type = cancer_type,
        path_threshold = "Tier4_0001",
        allele_frequency = "0.001"
    }
  }
  call Tasks.concatenateFiles_noheader as concat_tier0_001{
    input:
      files = T2_gsea_TIER0_001.out1,
      callset_name = biological_pathway + "_" + "001"
  }
  call Tasks.concatenateFiles_noheader as concat_tier0_0001{
    input:
      files = T2_gsea_TIER0_0001.out1,
      callset_name = biological_pathway + "_" + "0001"
  }
  call Tasks.concatenateFiles_noheader as concat_tier2_001{
    input:
      files = T2_gsea_TIER2_001.out1,
      callset_name = biological_pathway + "_" + "001"
  }
  call Tasks.concatenateFiles_noheader as concat_tier2_0001{
    input:
      files = T2_gsea_TIER2_0001.out1,
      callset_name = biological_pathway + "_" + "0001"
  }
  call Tasks.concatenateFiles_noheader as concat_tier3_001{
    input:
      files = T2_gsea_TIER3_001.out1,
      callset_name = biological_pathway + "_" + "001"
  }
  call Tasks.concatenateFiles_noheader as concat_tier3_0001{
    input:
      files = T2_gsea_TIER3_0001.out1,
      callset_name = biological_pathway + "_" + "0001"
  }
  call Tasks.concatenateFiles_noheader as concat_tier4_001{
    input:
      files = T2_gsea_TIER4_001.out1,
      callset_name = biological_pathway + "_" + "001"
  }
  call Tasks.concatenateFiles_noheader as concat_tier4_0001{
    input:
      files = T2_gsea_TIER4_0001.out1,
      callset_name = biological_pathway + "_" + "0001"
  }

  call Tasks.concatenateFiles_noheader as concat_final{
    input:
      files = [concat_tier0_001.out2,concat_tier0_0001.out2,concat_tier2_001.out2,concat_tier2_0001.out2,concat_tier3_001.out2,concat_tier3_0001.out2,concat_tier4_001.out2,concat_tier4_0001.out2],
      callset_name = "cosmic_genes"
  }
  call Tasks.copy_file_to_storage{
    input:
      text_file = concat_final.out1,
      output_dir = analysis_5_output_dir
  }

}

task T0_get_path_genes{
  input {
    File cosmic_tsv
    String cancer_type
  }
  command <<<
  set -euxo pipefail
  grep -E '~{cancer_type}|all' ~{cosmic_tsv} | cut -f1 > gene.list
  >>>
  output{
    File out1 = "gene.list"
  }
  runtime {
    preemptible:3
    docker: "ubuntu:latest"
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
      [f"{g}_{suffix}" for g in genes for suffix in ["Tier3_001"]]
  )]
  df_050.to_csv("tier3_001.tsv",sep='\t',index=False)

  # Filter where gene_impact is in genes
  df_075 = df[df['gene_impact'].isin(
      [f"{g}_{suffix}" for g in genes for suffix in ["Tier4_001"]]
  )]
  df_075.to_csv("tier4_001.tsv",sep='\t',index=False)

  # Filter where gene_impact is in genes
  df_050 = df[df['gene_impact'].isin(
      [f"{g}_{suffix}" for g in genes for suffix in ["Tier3_0001"]]
  )]
  df_050.to_csv("tier3_0001.tsv",sep='\t',index=False)

  # Filter where gene_impact is in genes
  df_075 = df[df['gene_impact'].isin(
      [f"{g}_{suffix}" for g in genes for suffix in ["Tier4_0001"]]
  )]
  df_075.to_csv("tier4_0001.tsv",sep='\t',index=False)

  # Filter where gene_impact is in genes
  df_VUS = df[df['gene_impact'].isin(
      [f"{g}_{suffix}" for g in genes for suffix in ["Tier2_001"]]
  )]
  df_VUS.to_csv("tier2_001.tsv",sep='\t',index=False)

  # Filter where gene_impact is in genes
  df_VUS = df[df['gene_impact'].isin(
      [f"{g}_{suffix}" for g in genes for suffix in ["Tier2_0001"]]
  )]
  df_VUS.to_csv("tier2_0001.tsv",sep='\t',index=False) 
  # Filter where gene_impact is in genes
  df_VUS = df[df['gene_impact'].isin(
      [f"{g}_{suffix}" for g in genes for suffix in ["Tier0_001"]]
  )]
  df_VUS.to_csv("tier0_001.tsv",sep='\t',index=False)

  # Filter where gene_impact is in genes
  df_VUS = df[df['gene_impact'].isin(
      [f"{g}_{suffix}" for g in genes for suffix in ["Tier0_0001"]]
  )]
  df_VUS.to_csv("tier0_0001.tsv",sep='\t',index=False)
  CODE
  >>>
  output{
    File out1 = "tier3_001.tsv"
    File out2 = "tier3_0001.tsv"
    File out3 = "tier4_001.tsv"
    File out4 = "tier4_0001.tsv"
    File out5 = "tier2_001.tsv"
    File out6 = "tier2_0001.tsv"
    File out7 = "tier0_001.tsv"
    File out8 = "tier0_0001.tsv"
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
  import statsmodels.api as sm
  from scipy.stats import zscore
  import numpy as np
  import scipy.stats as stats

  # Load the filtered data
  df = pd.read_csv("~{variants_tsv}", sep="\t",index_col=False)

  # Drop potential duplicates before summing
  df = df.drop_duplicates(subset=["gene_impact"])
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

  # Define covariates to zscore
  covariates = [f"PC{i}" for i in range(1, 5)]
  merged[covariates] = merged[covariates].apply(zscore)
  merged['num_pathogenic_variants_bin'] = (merged['num_pathogenic_variants'] >= 1).astype(int)
  #merged = merged.rename(columns={
  #  'num_pathogenic_variants_bin': 'num_pathogenic_variants',
  #})
  covariates.append('num_pathogenic_variants')

  # Add 'sex_binary' to covariates only if it varies
  if merged['sex_binary'].nunique() > 1:
      covariates.append('sex_binary')

  # Covariates
  X = merged[covariates]

  # Add intercept
  X = sm.add_constant(X)
  y = merged["is_case"]

  # Case/control summaries
  cases_vals = merged.loc[merged["is_case"] == 1, "num_pathogenic_variants"]
  controls_vals = merged.loc[merged["is_case"] == 0, "num_pathogenic_variants"]

  def five_num(x):
      if x.empty:
          return "NA\tNA\tNA\tNA\tNA"
      q1, q2, q3 = np.percentile(x, [25, 50, 75])
      return f"{x.min():.2f}\t{q1:.2f}\t{q2:.2f}\t{q3:.2f}\t{x.max():.2f}"


  cases_summary = five_num(cases_vals)
  controls_summary = five_num(controls_vals)

  try:
      # Fit logistic regression
      model = sm.Logit(y, X)
      result = model.fit(disp=False)

      coef = result.params['num_pathogenic_variants']              # beta
      se = result.bse['num_pathogenic_variants']                   # standard error
      ci_low, ci_high = result.conf_int().loc['num_pathogenic_variants']  # 95% CI
      pval = result.pvalues['num_pathogenic_variants']             # p-value

      row = f"~{cancer_type}\t~{path_threshold}\t~{allele_frequency}\t{pval:.3e}\t{coef:.4f}\t{ci_low:.4f}\t{ci_high:.4f}\t{cases_summary}\t{controls_summary}\n"

  except Exception:
      # If it fails, return NA line
      row = f"~{cancer_type}\t~{path_threshold}\t~{allele_frequency}\tNA\tNA\tNA\tNA\t{cases_summary}\t{controls_summary}\n"

  # Write header + result
  with open("logistic_regression_summary.txt", "w") as f:
      f.write("cancer\tPathogenic_Threshold\tAF\tp_value\tbeta\t95%_CI_lower\t95%_CI_upper\tcases_min\tcases_Q1\tcases_Q2\tcases_Q3\tcases_max\tcontrols_min\tcontrols_Q1\tcontrols_Q2\tcontrols_Q3\tcontrols_max\n")
      f.write(row)

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
