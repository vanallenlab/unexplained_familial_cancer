# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_4F_PRS_PENETRANCE {
  input{
    String analysis_4b_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/ADJUSTED_PRS/"
    String analysis_4f_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/PRS_PENETRANCE/"
    File phenotype_data = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/dfci-ufc.aou.phenos.v2.tsv.gz" 

    File PGS_ID = "PGS000783"
    File cancer_type = "breast"
  }
  String sample_data = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/"
  File specific_cohort_sample_data = sample_data + cancer_type + "/" + cancer_type + ".metadata"
  String adjusted_prs_file = analysis_4b_output_dir + cancer_type + "." + PGS_ID + ".pgs"

  call T1_prs_grouping_for_boxplot {
    input:
      prs_file = adjusted_prs_file,
      metadata_file = specific_cohort_sample_data,
      cancer_type = cancer_type,
      phenotype_data = phenotype_data
  }
}

task T1_prs_grouping_for_boxplot {
  input {
    File prs_file      # e.g. breast.PGS000783.pgs
    File metadata_file
    File phenotype_data
    String cancer_type
  }

  command <<<
  python3 <<CODE
  import pandas as pd

  # Load input
  prs = pd.read_csv("~{prs_file}", sep="\t",index_col=False)
  meta = pd.read_csv("~{metadata_file}", sep="\t",index_col=False)
  phenotype = pd.read_csv("~{phenotype_data}",sep='\t',index_col=False)

  # Ensure consistent sample IDs
  prs['sample'] = prs['sample'].astype(str)
  meta['sample'] = meta['original_id'].astype(str)
  phenotype['sample'] = phenotype['Sample'].astype(str)
 
  merged = meta.merge(prs, on='sample')
  merged = merged.merge(phenotype, on=['sample','original_dx'])

  # Initialize group
  merged['group'] = "NA"

  # Normalize to lowercase for safer matching
  cancer = "~{cancer_type}".lower()

  # Controls
  merged.loc[merged['original_dx'].str.lower() == "control", 'group'] = "Control"

  # Sporadic: original_dx contains cancer_type, family_dx does not
  mask_sporadic = (
      merged['original_dx'].fillna("").str.lower().str.contains(cancer) &
      ~merged['family_dx'].fillna("").str.lower().str.contains(cancer)
  )
  merged.loc[mask_sporadic, 'group'] = f"Sporadic_{cancer.capitalize()}"

  # Familial: original_dx contains cancer_type, family_dx does too
  mask_familial = (
      merged['original_dx'].fillna("").str.lower().str.contains(cancer) &
      merged['family_dx'].fillna("").str.lower().str.contains(cancer)
  )
  merged.loc[mask_familial, 'group'] = f"Familial_{cancer.capitalize()}"

  # Keep only needed columns for plotting
  out_df = merged[['sample', 'PGS', 'group']]

  # Write output
  out_df.to_csv(f"prs_groups_{cancer}.tsv", sep="\t", index=False)


  CODE
  >>>

  output {
    File prs_groups = "prs_groups_~{cancer_type}.tsv"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "8G"
    preemptible: 3
  }
}

