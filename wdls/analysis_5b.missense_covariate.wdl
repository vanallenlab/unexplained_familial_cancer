# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_5B_MISSENSE_COVARIATE {
  input {
    File step_10_cpg_output
    File analysis_5_output
    File cosmic_ufc
    File sample_list
    Array[String] cancer_types = ["breast"]
  }
  
  call T1_filter_step_5 {
    input:
      analysis_5_output = analysis_5_output
  }
  scatter (cancer_type in cancer_types) {
    File sample_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".list"
    String output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/ANALYSIS_5/" + cancer_type + ".ppv_missense"
    call T2_normalize_silico_missense {
      input:
        step_10_cpg_output = step_10_cpg_output,
        cosmic_ufc = cosmic_ufc,
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
    out1 = "analysis_5_filtered.tsv"
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
  # --- Load data ---
  df = pd.read_csv("~{step_10_cpg_output}", sep="\t",index_col=False)
  analysis5_df = pd.read_csv("~{t1_output}", sep="\t",index_col=False)
  analysis5_df = analysis5_df[analysis5_df['cancer'] == "~{cancer_type}"]


  # --- Melt into long form and filter where value == 1 ---
  df_long = (
      df.stack()
        .reset_index()
        .rename(columns={"level_1": "patient", 0: "value"})
  )

  # Keep only where the variant count == 1
  df_filtered = df_long[df_long["value"] == 1][["patient", criteria_col]]

  # --- Save output ---
  df_filtered.to_csv("patient_criteria.tsv", sep="\t", index=False)

  cosmic_df = pd.read_csv("~{cosmic_ufc}",sep='\t',index_col=False,names=['gene','cancer_type'])
  cosmic_genes = cosmic_df[cosmic_df['cancer_type'].str.contains("~{cancer_type}|all", na=False)]['gene'].tolist()


  # Make sure that we filter patient_criteria.tsv to just rows that have gene_Tier~{tier}_~{AF}
  tier = analysis5_df['Pathogenic_Threshold'][0]
  AF = analysis5_df['AF'][0]

  # get samples
  with open("~{sample_list}") as f:
    all_samples = [line.strip() for line in f if line.strip()]

  #match_strings = {f"{gene}_Tier{tier}_{AF}" for gene in cosmic_genes}
  match_strings = {f"{gene}_{tier}_{AF}" for gene in cosmic_genes}

  # --- Step 3: Loop through patient_criteria.tsv ---
  # Column 0 = patient ID, column 1 = variant/gene info
  patients_with_hit = set()

  for _, row in df_filtered.iterrows():
      patient_id = row[0]
      gene_label = str(row[1])
      if gene_label in match_strings:
          patients_with_hit.add(patient_id)

  # --- Step 4: Build output DataFrame ---
  output_df = pd.DataFrame({
      'original_id': all_samples,
      'PPV_Missense': [1 if s in patients_with_hit else 0 for s in all_samples]
  })

  # --- Step 5: Write to TSV ---
  output_df.to_csv("PPV_Missense_flags.tsv", sep='\t', index=False)

  CODE
  >>>
  output {
    File out1 = "patient_criteria.tsv"
    File out2 = "PPV_Missense_flags.tsv"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}
