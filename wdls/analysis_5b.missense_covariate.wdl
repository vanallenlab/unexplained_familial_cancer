# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_5B_MISSENSE_COVARIATE {
  input {
    File step_10_cpg_output = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/STEP_10_VISUALIZE_VEP/77bd97e7-cc41-422e-948a-ee33a709c73b/call-merge_variant_counts/shard-4/ufc.cpg.variant_counts.tsv"
    #File analysis_5_output
    File cosmic_tsv = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cosmic_ufc.tsv"
    #File sample_list
    #Array[String] cancer_types = ["breast"]
    Array[String] cancer_types = ["basal_cell","bladder","breast","colorectal","hematologic","kidney","lung","melanoma","neuroendocrine","nervous","non-hodgkins","ovary","prostate","sarcoma","squamous_cell","thyroid","uterus","cervix"]
  }
  
  call T1_filter_step_5 {
    input:
      analysis_5_output = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/ANALYSIS_5_GSEA/4b4224e2-7b5c-407f-a0e4-38c57ae3059d/call-concatenateFiles_noheader/all_cpg_001.tsv"
  }
  scatter (cancer_type in cancer_types) {
    File sample_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".list"
    String output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_5_GSEA/"
    call T2_normalize_silico_missense {
      input:
        step_10_cpg_output = step_10_cpg_output,
        cosmic_ufc = cosmic_tsv,
        cancer_type = cancer_type,
        sample_list = sample_list,
        t1_output = T1_filter_step_5.out1
    }
    call Tasks.copy_file_to_storage{
     input:
        text_file = T2_normalize_silico_missense.out1,
        output_dir = output_dir
    }
  }
}

task T1_filter_step_5 {
  input {
    File analysis_5_output
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  df = pd.read_csv("~{analysis_5_output}",sep='\t',index_col=False)

  # For each cancer, get the index of the row with the lowest p_value
  idx = df.groupby('cancer')['p_value'].idxmin()

  # Keep only those rows
  df_lowest_p = df.loc[idx].reset_index(drop=True)

  df_lowest_p.to_csv("analysis_5_filtered.tsv",sep='\t',index=False)
  CODE
  >>>
  output {
    File out1 = "analysis_5_filtered.tsv"
  }
  runtime {
    docker:"vanallenlab/pydata_stack"
    preemptible:3
  }
} 


task T2_normalize_silico_missense {
  input {
    File step_10_cpg_output
    File t1_output
    File cosmic_ufc
    File sample_list
    String cancer_type
  }
  command <<<
  python3 <<CODE
  import pandas as pd

  # --- Load input files ---
  df = pd.read_csv("~{step_10_cpg_output}", sep="\t",index_col=False)         # step 10 output, wide format
  analysis_df = pd.read_csv("~{t1_output}", sep="\t",index_col=False)
  analysis_df = analysis_df[analysis_df['cancer'].str.contains("~{cancer_type}", na=False)]

  # Pull tier + AF label used in the column identifiers
  tier = analysis_df['Pathogenic_Threshold'].iloc[0]   # e.g., "Tier3"

  # Load COSMIC lookup list
  cosmic_df = pd.read_csv("~{cosmic_ufc}", sep="\t", index_col=False,names=['gene', 'cancer_type'])
  cosmic_genes = cosmic_df[
      cosmic_df['cancer_type'].str.contains("~{cancer_type}|all", na=False)]['gene'].tolist()
  print("check1")
  # Construct the strings we need to look for in the "gene_impact" column
  match_strings = {f"{gene}_{tier}" for gene in cosmic_genes}

  # Filter df down to rows (genes) that match our PPV missense criteria
  df_filtered = df[df['gene_impact'].isin(match_strings)]

  # Sum variant counts per patient (each patient = a column, so sum rows)
  patient_counts = df_filtered.iloc[:, 1:].sum(axis=0)
  print("chek2")
  # Turn into final output format
  output_df = (
      patient_counts
      .reset_index()
      .rename(columns={"index": "original_id", 0: "num_ppv_missense"})
  )

  # --- Save output ---
  output_df.to_csv("~{cancer_type}.num_ppv_missense.tsv", sep="\t", index=False)

  CODE
  >>>
  output {
    File out1 = "~{cancer_type}.num_ppv_missense.tsv"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}
