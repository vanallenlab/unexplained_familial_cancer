# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0

workflow ANALYSIS_6A_DO_IT_ALL {
  input {
    String cancer_type = "breast"
    String workspace_bucket = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"

    # SAIGE_RESULTS
    File analysis_3b_saige_results = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_3_SAIGE_GENE/results/breast.patient_report"

    # PRS metrics
    #Float analysis_4_p_cutoff = 0.05
    File prs_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/ADJUSTED_PRS/breast.PGS000783.pgs"

    # Damaging Missense
    File step_10_cpg_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_10_VISUALIZE_VEP/v2/ufc.cpg.variant_counts.tsv.gz"
    File cosmic_tsv = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cosmic_ufc.tsv"
    #File analysis5_results 
    #Array[String] tier_of_interest
    #Array[String] analysis_5_afs
    #File step10_results

    # ANALYSIS 1
    #File? homozygosity_tsv
    #Array[File]? analysis_1_genes

    # SV ANALYSIS?
    # File sv_tsv  # optional - commented out
    
  }

  # Metadata for cancer type
  File metadata_tsv = "gs://" + workspace_bucket + "/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".metadata"
  File sample_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".list"

  call T1_normalize_saige {
    input: 
      saige_results_tsv=analysis_3b_saige_results,
      sample_list = sample_list
  }

  call T2_normalize_silico_missense {
    input:
      step_10_cpg_output = step_10_cpg_file, 
      cosmic_ufc = cosmic_tsv,
      tier = "REVEL075",
      AF = "001",
      sample_list = sample_list,
      cancer_type = "breast" 
  }

  call aggregate_tsvs {
    input:
      saige_tsv=T1_normalize_saige.out1,
      prs_file = prs_file,
      ppv_missense = T2_normalize_silico_missense.out2,
      metadata = metadata_tsv
  }

  call T5_log_reg {
    input:
      input_tsv=aggregate_tsvs.out1
  }

  #call normalize_analysis5 {
  #  input: analysis5_tsv=analysis5_tsv
  #}

  #call normalize_step10 {
  #  input: step10_tsv=step10_tsv
  #}


  #call normalize_homozygosity {
  #  input: homoz_tsv=homozygosity_tsv,
  #         gene_of_interest=gene_of_interest
  #}

  # call normalize_svs {
  #   input: sv_tsv=sv_tsv
  # }

  #call aggregate {
  #  input:
  #    saige_normalized=normalize_saige.saige_normalized,
  #    analysis5_normalized=normalize_analysis5.analysis5_normalized,
  #    step10_normalized=normalize_step10.step10_normalized,
  #    prs_normalized=normalize_prs.prs_normalized,
  #    metadata_normalized=normalize_metadata.metadata_normalized,
  #    homozygosity_normalized=normalize_homozygosity.homozygosity_normalized,
  #    gene_of_interest=gene_of_interest,
  #    cohort_name=cohort_name
  #}

}

task T1_normalize_saige {
  input {
    File saige_results_tsv
    File sample_list
  }

  command <<<
  set -euxo pipefail
  python3 <<CODE
  import pandas as pd

  # --- Load data ---
  df = pd.read_csv("~{saige_results_tsv}", sep="\t", index_col=False)

  # --- Basic cleanup ---
  df = df[['original_id', 'gene']].dropna().drop_duplicates()

  # --- Pivot to wide format ---
  wide = (
      df.assign(value=1)
        .pivot_table(index='original_id', columns='gene', values='value', fill_value=0)
        .reset_index()
  )
  wide['original_id'] = wide['original_id'].astype(str)

  # --- Load sample list ---
  with open("~{sample_list}", 'r') as f:
      all_samples = [str(line.strip()) for line in f if line.strip()]

  # --- Ensure all samples are present ---
  existing_samples = set(wide['original_id'])
  missing_samples = [s for s in all_samples if s not in existing_samples]

  if missing_samples:
      # Create a DataFrame of zeros for missing samples
      zero_df = pd.DataFrame(0, index=range(len(missing_samples)), columns=wide.columns)
      zero_df['original_id'] = missing_samples
      wide = pd.concat([wide, zero_df], ignore_index=True)

  # --- Fill any remaining NaNs and ensure correct type ---
  wide = wide.fillna(0).astype({col: int for col in wide.columns if col != 'original_id'})

  # --- Save ---
  wide.to_csv("saige_normalized.tsv", sep="\t", index=False)
  CODE
  >>>

  output {
    File out1 = "saige_normalized.tsv"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}

task T2_normalize_missense_burden {
  input {
    File missense_results_tsv
    File step_10_output_file
    File cosmic_ufc
    File sample_list
    String cancer_type
  }

  command <<<
  set -euxo pipefail
  python3 <<CODE
  import pandas as pd

  # --- Load data ---
  df = pd.read_csv("~{missense_results_tsv}", sep="\t", index_col=False)
  df = df[df['cancer'] == "~{cancer_type}"]
  df = df.loc[df.groupby('cancer')['p_value'].idxmin()].reset_index(drop=True)
  criteria = df[:,0]['Pathogenic_Threshold']

  cosmic_df = pd.read_csv("~{cosmic_ufc}",sep='\t',index=False)
  cosmic_gene_list = cosmic_df[cosmic_df['cancer_types'].isin(['~{cancer_type}','all'])]['gene'].tolist()
  CODE
  >>>

  output {
    File out1 = "saige_normalized.tsv"
  }

  runtime {                                                                                                                
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}

# Structural variants normalization task - commented out per request
# task normalize_svs {
#   input {
#     File sv_tsv
#   }
#   command <<~PY
#     set -euo pipefail
#     python - <<'PYCODE'
# import pandas as pd
# df = pd.read_csv("~{sv_tsv}", sep="\t", dtype=str)
# # normalize to sample,gene,sv_type,info
# rows=[]
# for _,r in df.iterrows():
#     sample = r.get("sample","")
#     gene = r.get("gene","")
#     svtype = r.get("sv_type", r.get("type",""))
#     rows.append({"sample":sample,"gene":gene,"svtype":svtype,"info":r.to_json()})
# pd.DataFrame(rows).to_csv("svs_normalized.tsv", sep="\t", index=False)
# PYCODE
#   PY
#   output {
#     File svs_normalized = "svs_normalized.tsv"
#   }
#   runtime {
#     docker: "vanallenlab/g2c_pipeline"
#     memory: "8G"
#     cpu: 1
#   }
# }

task aggregate_tsvs {
  input {
    File saige_tsv
    File prs_file
    File ppv_missense
    File metadata
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  df1 = pd.read_csv("~{saige_tsv}",sep='\t',index_col=False) #original_id
  df2 = pd.read_csv("~{prs_file}",sep='\t',index_col=False) #sample
  df3 = pd.read_csv("~{metadata}",sep='\t',index_col=False) #original_id
  df4 = pd.read_csv("~{ppv_missense}",sep='\t',index_col=False) #original_id
  # --- Convert IDs to strings ---
  df1['original_id'] = df1['original_id'].astype(str)
  df2['sample'] = df2['sample'].astype(str)
  df3['original_id'] = df3['original_id'].astype(str)
  df4['original_id'] = df4['original_id'].astype(str)

  # --- Subset metadata to relevant columns ---
  keep_cols = ['original_id', 'PC1', 'PC2', 'PC3', 'PC4', 'inferred_sex','original_dx']
  df3 = df3[[c for c in keep_cols if c in df3.columns]]

  # --- Standardize naming to match for merge ---
  df2 = df2.rename(columns={'sample': 'original_id'})

  # --- Merge all dataframes on 'original_id' ---
  merged = (
      df1.merge(df2, on='original_id', how='outer')
         .merge(df3, on='original_id', how='outer')
         .merge(df4, on='original_id', how='outer')
  )

  # --- Save final combined table ---
  merged.to_csv("combined_risk_factors.tsv", sep='\t', index=False)

  CODE
  >>>
  output {
    File out1 = "combined_risk_factors.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}

task T5_log_reg {
  input{
    File input_tsv
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  import statsmodels.api as sm
  import numpy as np

  def nagelkerke_r2(model):
      """
      Compute Nagelkerke pseudo-R² for a fitted statsmodels Logit model.
      """
      llf = model.llf        # log-likelihood of fitted model
      llnull = model.llnull  # log-likelihood of null model
      n = model.nobs

      # McFadden/Nagelkerke style
      numerator = 1 - np.exp((2 * (llnull - llf)) / n)
      denominator = 1 - np.exp((2 * llnull) / n)
      r2_nagelkerke = numerator / denominator if denominator != 0 else np.nan

      return r2_nagelkerke

  def attributable_fraction(df, outcome_col, predictor_cols=[], covariates=[]):
      """
      Compute fraction of variance explained by predictor_col on outcome_col
      df: DataFrame with all variables
      outcome_col: binary 0/1 column
      predictor_col: the variable you want to evaluate
      covariates: list of covariate column names to include
      """
      # Prepare design matrices
      X_full = df[predictor_cols + covariates]
      X_full.to_csv("test.tsv",sep='\t',index=False)
      X_full = sm.add_constant(X_full)
      X_reduced = df[covariates] if covariates else pd.DataFrame({'const': [1]*len(df)})
      X_reduced = sm.add_constant(X_reduced, has_constant='add')
      y = df[outcome_col]

      # Fit logistic models
      full_model = sm.Logit(y, X_full).fit(disp=0)
      reduced_model = sm.Logit(y, X_reduced).fit(disp=0)
      # Compute Nagelkerke R2
      r2_full = nagelkerke_r2(full_model)
      r2_reduced = nagelkerke_r2(reduced_model)

      # Fraction of variance explained by predictor
      fraction_explained = r2_full - r2_reduced

      return {
          'R2_full': r2_full,
          'R2_reduced': r2_reduced,
          'fraction_explained': fraction_explained
      }

  df = pd.read_csv("~{input_tsv}",sep='\t',index_col=False)

  # Initialize covariates
  covariates = ['PC1', 'PC2', 'PC3', 'PC4']

  # Check if inferred_sex has more than 1 unique value
  if df['inferred_sex'].nunique() > 1:
      # Create sex_binary: 1 = male, 0 = female
      df['sex_binary'] = df['inferred_sex'].str.lower().map({'male': 1, 'female': 0})
      covariates.append('sex_binary')

  df['case'] = df['original_dx'].apply(lambda x: 0 if x.lower() == 'control' else 1)
  predictor_cols = ['PGS','PPV_Missense']
  #predictor_cols = [col for col in df.columns if col not in covariates + ['original_id', 'case','inferred_sex','original_dx']]
  result = attributable_fraction(df, outcome_col='case', predictor_cols=predictor_cols, covariates=covariates)

  # Convert result dict to DataFrame for easier saving
  results_df = pd.DataFrame([result])
  results_df['predictors'] = ','.join(predictor_cols)

  # Reorder columns
  results_df = results_df[['predictors','R2_full','R2_reduced','fraction_explained']]

  # Save to TSV
  results_df.to_csv("attributable_fraction_results.tsv", sep='\t', index=False)

  CODE
  >>>
  output {
    File out1 = "attributable_fraction_results.tsv"
  }
  runtime{
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task T2_normalize_silico_missense {
  input {
    File step_10_cpg_output
    File cosmic_ufc
    String tier # 3 or 4
    String AF #001 or 0001
    File sample_list
    String cancer_type
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  # --- Load data ---
  df = pd.read_csv("~{step_10_cpg_output}", sep="\t",index_col=False)

  # Assume first column contains the variant criteria (like "ATM_REVEL_050_001")
  criteria_col = df.columns[0]
  df = df.set_index(criteria_col)

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
  tier = "~{tier}"
  AF = "~{AF}"

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
