# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_4C_PRS {
  input {
    Array[String] ufc_cancer_type = ["breast","breast","breast","breast","bladder","bladder","bladder","cervix","cervix","colorectal","colorectal","colorectal","colorectal","uterus","uterus","uterus","kidney","kidney","kidney","leukemia","lung","lung","lung","lung","melanoma","melanoma","melanoma","non-hodgkins","non-hodgkins","ovary","ovary","ovary","ovary","pancreas","pancreas","pancreas","prostate","prostate","prostate","prostate","thyroid","brain","esophagus"]
    Array[String] aou_cancer_type = ["breast","breast","breast","breast","bladder","bladder","bladder","cervix","cervix","colorectal","colorectal","colorectal","colorectal","uterus","uterus","uterus","kidney","kidney","kidney","blood_soft_tissue","lung","lung","lung","lung","skin","skin","skin","blood_soft_tissue","blood_soft_tissue","ovary","ovary","ovary","ovary","pancreas","pancreas","pancreas","prostate","prostate","prostate","prostate","thyroid","brain","esophagus"]
    Array[String] PGS_IDS = ["PGS000783","PGS003380","PGS004242","PGS004688","PGS004241","PGS000782","PGS004687","PGS000784","PGS003389","PGS000785","PGS003386","PGS004243","PGS004689","PGS000786","PGS003381","PGS004244","PGS000787","PGS004690","PGS004245","PGS000788","PGS000789","PGS003391","PGS004246","PGS004691","PGS000790","PGS003382","PGS004247","PGS000791","PGS004248","PGS000793","PGS003385","PGS004249","PGS004692","PGS000794","PGS004250","PGS004693","PGS000795","PGS003383","PGS004251","PGS004694","PGS000797","PGS003384","PGS003388"]
    File analysis_4_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/"
    File phenotype_data = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/dfci-ufc.aou.phenos.v2.tsv.gz"  
  }

  scatter(i in range(length(PGS_IDS) - 0)){
    File prs_file = analysis_4_dir + "ADJUSTED_PRS/" + ufc_cancer_type[i] + "." + PGS_IDS[i] + ".pgs"
    File metadata = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + ufc_cancer_type[i] + "/" + ufc_cancer_type[i] + ".metadata"

    call T1_analyze_inheritance {
      input:
        phenotype_data = phenotype_data,
        prs_file = prs_file,
        metadata = metadata,
        ufc_cancer_type = ufc_cancer_type[i],
        aou_cancer_type = aou_cancer_type[i],
        PGS_ID = PGS_IDS[i]
    }
    call Tasks.copy_file_to_storage {
      input:
        text_file = T1_analyze_inheritance.out1,
        output_dir = analysis_4_dir + "FAMILIAL/"
    }
  }

}

task T1_analyze_inheritance {
  input {
    File metadata
    File phenotype_data
    File prs_file
    String PGS_ID
    String ufc_cancer_type
    String aou_cancer_type
  }
  command <<<
  set -euxo pipefail

  python3 <<CODE
  import pandas as pd
  import matplotlib.pyplot as plt

  # --- Load files ---
  prs = pd.read_csv("~{prs_file}", sep="\t", index_col=False)
  meta = pd.read_csv("~{metadata}", sep="\t", index_col=False)
  #phenos = pd.read_csv("~{phenotype_data}", sep="\t", usecols=['Sample','family_dx'], index_col=False)

  # --- Ensure consistent sample IDs ---
  prs['sample'] = prs['sample'].astype(str)
  meta['original_id'] = meta['original_id'].astype(str)
  #phenos['Sample'] = phenos['Sample'].astype(str)

  # --- Merge PRS with metadata ---
  merged = prs.merge(meta, left_on="sample", right_on="original_id", how="inner")

  # --- Merge with family history ---
  #merged = merged.merge(phenos, left_on="sample", right_on="Sample", how="left")

  # --- Initialize group ---
  merged['group'] = "NA"

  # --- Normalize cancer string ---
  ufc_cancer = "~{ufc_cancer_type}"
  aou_cancer = "~{aou_cancer_type}"

  # --- Define groups ---
  # Controls
  merged.loc[merged['original_dx'].str.lower() == "control", 'group'] = "Control"

  # Sporadic: original_dx contains breast, family_dx does not
  mask_sporadic = (
      merged['original_dx'].fillna("").str.lower().str.contains(ufc_cancer) &
      ~merged['family_dx'].fillna("").str.lower().str.contains(aou_cancer)
  )
  merged.loc[mask_sporadic, 'group'] = f"Isolated {ufc_cancer.capitalize()}"

  # Familial: original_dx contains breast, family_dx contains breast
  mask_familial = (
      merged['original_dx'].fillna("").str.lower().str.contains(ufc_cancer) &
      merged['family_dx'].fillna("").str.lower().str.contains(aou_cancer)
  )
  merged.loc[mask_familial, 'group'] = f"Familial {ufc_cancer.capitalize()}"

  # --- Keep only needed columns ---
  out_df = merged[['sample', 'PGS', 'group']]

  # --- Keep only PGS and group ---
  subset = merged[['PGS', 'group']]

  # --- Print to stdout ---
  subset.to_csv("~{ufc_cancer_type}.~{PGS_ID}.analysis_4c.tsv",sep='\t',index=False)
  CODE
  >>>

  output {
    File out1 = "~{ufc_cancer_type}.~{PGS_ID}.analysis_4c.tsv"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "2G"
    preemptible: 3
  }
}

