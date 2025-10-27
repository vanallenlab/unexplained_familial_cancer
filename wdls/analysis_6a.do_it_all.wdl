version 1.0

workflow ANALYSIS_6A_DO_IT_ALL {
  input {
    String cancer_type = "thyroid"
    String workspace_bucket = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
    File sample_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/thyroid/thyroid.list"
    # SAIGE_RESULTS
    File analysis_3b_saige_results = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_3_SAIGE_GENE/results/thyroid.patient_report"

    # PRS metrics
    #Float analysis_4_p_cutoff = 0.05
    File prs_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/ADJUSTED_PRS/thyroid.PGS000795.pgs"

    # Damaging Missense
    #Float analysis_5_p_cutoff = 0.05
    #File step_10_cpg_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_10_VISUALIZE_VEP/v2/ufc.cpg.variant_counts.tsv.gz""
    #File cosmic_tsv = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cosmic_ufc.tsv"
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

  call T1_normalize_saige {
    input: 
      saige_results_tsv=analysis_3b_saige_results,
      sample_list = sample_list
  }

  call aggregate_tsvs {
    input:
      saige_tsv=T1_normalize_saige.out1,
      prs_file = prs_file,
      metadata = metadata_tsv
  }
  #call normalize_analysis5 {
  #  input: analysis5_tsv=analysis5_tsv
  #}

  #call normalize_step10 {
  #  input: step10_tsv=step10_tsv
  #}

  #call normalize_prs {
  #  input: prs_tsv=prs_tsv
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
    File metadata
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  df1 = pd.read_csv("~{saige_tsv}",sep='\t',index_col=False) #original_id
  df2 = pd.read_csv("~{prs_file}",sep='\t',index_col=False) #sample
  df3 = pd.read_csv("~{metadata}",sep='\t',index_col=False) #original_id

  # --- Convert IDs to strings ---
  df1['original_id'] = df1['original_id'].astype(str)
  df2['sample'] = df2['sample'].astype(str)
  df3['original_id'] = df3['original_id'].astype(str)

  # --- Subset metadata to relevant columns ---
  keep_cols = ['original_id', 'PC1', 'PC2', 'PC3', 'PC4', 'inferred_sex']
  df3 = df3[[c for c in keep_cols if c in df3.columns]]

  # --- Standardize naming to match for merge ---
  df2 = df2.rename(columns={'sample': 'original_id'})

  # --- Merge all dataframes on 'original_id' ---
  merged = (
      df1.merge(df2, on='original_id', how='outer')
         .merge(df3, on='original_id', how='outer')
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
