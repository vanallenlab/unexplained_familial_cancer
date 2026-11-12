
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow STEP_R5_ANCESTRY_ADJUST {
  input {
    File metadata_tsv
    File prs_tsv = "gs://dataproc-staging-wb-cordial-diamond-9893/AoU_Cancer_PRS/CDRv9_prs.tsv.gz"
  }

  call T1_ResidualizePRS {
    input:
      metadata_tsv = metadata_tsv,
      prs_tsv = prs_tsv
  }
  call Tasks.copy_file_to_storage {
    input:
      text_file = T1_ResidualizePRS.out1,
      output_dir = "gs://dataproc-staging-wb-cordial-diamond-9893/AoU_Cancer_PRS/"
  }
}

task T1_ResidualizePRS {
  input {
    File metadata_tsv
    File prs_tsv
  }

  command <<<
  python3 <<CODE
  import pandas as pd
  import numpy as np
  import statsmodels.api as sm

  metadata = pd.read_csv("~{metadata_tsv}", sep="\t", dtype={"person_id": str})

  prs = pd.read_csv("~{prs_tsv}", sep="\t", dtype={"person_id": str})

  covariates = ["PC1", "PC2", "PC3", "PC4", "PC5"]

  merged_df = metadata[["person_id"] + covariates].merge(prs, on="person_id", how="inner")

  pgs_columns = [col for col in prs.columns if col.startswith("PGS")]

  adjusted_df = merged_df[["person_id"]].copy()
  summary = []

  for pgs in pgs_columns:
    valid = merged_df[covariates + [pgs]].notna().all(axis=1)

    X = sm.add_constant(merged_df.loc[valid, covariates])
    y = merged_df.loc[valid, pgs]

    model = sm.OLS(y, X).fit()

    adjusted_df[pgs] = np.nan
    adjusted_df.loc[valid, pgs] = model.resid

  adjusted_df.to_csv("CDRv9_prs.pc5_adjusted.tsv.gz", sep="\t", index=False, compression="gzip")
  CODE
  >>>

  output {
    File out1 = "CDRv9_prs.pc5_adjusted.tsv.gz"
  }

  runtime {
    docker: "vanallenlab/gwc_pipeline"
    memory: "8 GB"
  }
}
