# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_4E_PRS_HISTORY {
  input {
    Array[String] ufc_cancer_type = ["breast","breast","breast","breast","bladder","bladder","bladder","cervix","cervix","colorectal","colorectal","colorectal","colorectal","uterus","uterus","uterus","kidney","kidney","kidney","leukemia","lung","lung","lung","lung","melanoma","melanoma","melanoma","non-hodgkins","non-hodgkins","ovary","ovary","ovary","ovary","pancreas","pancreas","pancreas","prostate","prostate","prostate","prostate","thyroid","brain","esophagus"]
    Array[String] aou_cancer_type = ["breast","breast","breast","breast","bladder","bladder","bladder","cervix","cervix","colorectal","colorectal","colorectal","colorectal","uterus","uterus","uterus","kidney","kidney","kidney","blood_soft_tissue","lung","lung","lung","lung","skin","skin","skin","blood_soft_tissue","blood_soft_tissue","ovary","ovary","ovary","ovary","pancreas","pancreas","pancreas","prostate","prostate","prostate","prostate","thyroid","brain","esophagus"]
    Array[String] PGS_IDS = ["PGS000783","PGS003380","PGS004242","PGS004688","PGS004241","PGS000782","PGS004687","PGS000784","PGS003389","PGS000785","PGS003386","PGS004243","PGS004689","PGS000786","PGS003381","PGS004244","PGS000787","PGS004690","PGS004245","PGS000788","PGS000789","PGS003391","PGS004246","PGS004691","PGS000790","PGS003382","PGS004247","PGS000791","PGS004248","PGS000793","PGS003385","PGS004249","PGS004692","PGS000794","PGS004250","PGS004693","PGS000795","PGS003383","PGS004251","PGS004694","PGS000797","PGS003384","PGS003388"]
    File analysis_4_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/"
  }
  Int negative_shards = 42
  scatter(i in range(length(PGS_IDS) - negative_shards)){
    File metadata = analysis_4_dir + "UFC_REFERENCE_FILES/" + aou_cancer_type + "_family/" + aou_cancer_type + "_family.metadata" 
    File prs_file = analysis_4_dir  + PGS_IDS[i] + ".raw.pgs"

    call T1_analyze_history_and_regress {
      input:
        metadata = metadata,
        prs_file = prs_file,
        ufc_cancer_type = ufc_cancer_type[i],
        aou_cancer_type = aou_cancer_type[i],
        PGS_ID = PGS_IDS[i]
    }
    call Tasks.copy_file_to_storage {
      input:
        text_file = T1_analyze_history_and_regress.out1,
        output_dir = analysis_4_dir + "FAMILIAL_4E/"
    }
  }

}

task T1_analyze_history_and_regress {
  input {
    File prs_file
    File metadata
    String PGS_ID
    String ufc_cancer_type
    String aou_cancer_type
  }
  command <<<
  set -euxo pipefail

  python3 <<CODE
  import pandas as pd
  from scipy.stats import zscore
  import statsmodels.api as sm

  # --- Load files ---
  prs = pd.read_csv("~{prs_file}", sep="\t", index_col=False)
  metadata = pd.read_csv("~{metadata}",sep='\t', usecols = ['original_id','inferred_sex'],index_col=False)
  metadata.rename(columns={'original_id':'sample'}, inplace=True)

  # --- Ensure consistent sample IDs ---
  prs['sample'] = prs['sample'].astype(str)
  metadata['sample'] = metadata['sample'].astype(str)

  # --- Merge with family history ---
  merged = prs.merge(metadata, left_on="sample", right_on="Sample", how="inner")

  # --- Initialize group ---
  merged['group'] = "control"

  # --- Normalize cancer string ---
  ufc_cancer = "~{ufc_cancer_type}"
  aou_cancer = "~{aou_cancer_type}"

  # --- Define groups ---

  # Define excluded-only cancer types
  excluded = ["unspecific", "basal_cell", "squamous_cell"]

  mask_not_familial = (
      ~merged['original_dx'].fillna("").str.lower().str.contains('|'.join(['control', ufc_cancer])) &
      # Keep only those who have *something else* besides the excluded terms
      merged['original_dx'].fillna("").str.lower().apply(
          lambda x: any(term not in excluded and term != "" for term in x.split(';'))
      ) &
      (merged['maternal_family_dx'].fillna("").str.lower().str.contains(aou_cancer) | merged['paternal_family_dx'].fillna("").str.lower().str.contains(aou_cancer))
  )
  merged.loc[mask_not_familial, 'group'] = f"Not-Inherited {ufc_cancer.capitalize()}"

  # Familial: original_dx contains breast, family_dx contains breast
  mask_familial = (
      merged['original_dx'].fillna("").str.lower().str.contains(ufc_cancer) &
      (merged['maternal_family_dx'].fillna("").str.lower().str.contains(aou_cancer) | merged['paternal_family_dx'].fillna("").str.lower().str.contains(aou_cancer))
  )
  merged.loc[mask_familial, 'group'] = f"Familial {ufc_cancer.capitalize()}"

  merged = merged[merged['group'].isin([
      f"Familial {ufc_cancer.capitalize()}",
      f"Not-Inherited {ufc_cancer.capitalize()}"
  ])]

  merged[['PC1','PC2','PC3','PC4','age']] = merged[['PC1','PC2','PC3','PC4','age']].apply(zscore)
  merged['sex_binary'] = merged['inferred_sex'].apply(lambda x: 1 if x == 'male' else 0)

  covariates = ['PC1', 'PC2', 'PC3', 'PC4']
  if merged['sex_binary'].nunique() > 1:
    covariates.append('sex_binary')
  # Set up regression: regress PGS ~ PC1 + PC2 + PC3 + PC4 + male_sex
  X = merged[covariates]
  X = sm.add_constant(X)
  y = merged['PGS']

  model = sm.OLS(y, X).fit()
  merged['PGS'] = model.resid

  # Calculate mean and std using only controls
  control_pgs = merged.loc[merged['group'] == f"control", 'PGS']
  control_mean = control_pgs.mean()
  control_std = control_pgs.std()

  # Standardize PGS for all samples using control-derived mean and std
  merged['PGS'] = (merged['PGS'] - control_mean) / control_std

  # --- Print to stdout ---
  merged[['PGS','group']].to_csv("~{ufc_cancer_type}.~{PGS_ID}.analysis_4e.tsv",sep='\t',index=False)
  CODE
  >>>

  output {
    File out1 = "~{ufc_cancer_type}.~{PGS_ID}.analysis_4e.tsv"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "2G"
    preemptible: 3
  }
}

