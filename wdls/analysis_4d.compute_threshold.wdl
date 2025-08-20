# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_4D_COMPUTE_THRESHOLDS {
  input {
    String PGS_ID = "PGS000795"
    String google_bucket = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
    String analysis_4_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/"
    String cancer_type = "prostate"
    File phenos_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/dfci-ufc.aou.phenos.v2.tsv.gz"
  }
  File pgs_file = google_bucket + "/ANALYSIS_4_PRS/" + cancer_type + "." + PGS_ID + ".pgs"

  call T1_compute_thresholds {
    input:
      pgs_file = pgs_file,
      phenos_file = phenos_file,
      cancer_type = cancer_type,
      PGS_ID = PGS_ID
  }
  call Tasks.copy_file_to_storage as copy0{
   input:
      text_file = T1_compute_thresholds.or2_thresholds,
      output_dir = analysis_4_output_dir
  }
  call Tasks.copy_file_to_storage as copy1{
   input:
      text_file = T1_compute_thresholds.or5_thresholds,
      output_dir = analysis_4_output_dir
  }
  call Tasks.copy_file_to_storage as copy2{
   input:
      text_file = T1_compute_thresholds.or2_cases,
      output_dir = analysis_4_output_dir
  }
  call Tasks.copy_file_to_storage as copy3{
   input:
      text_file = T1_compute_thresholds.or2_controls,
      output_dir = analysis_4_output_dir
  }
  call Tasks.copy_file_to_storage as copy4{
   input:
      text_file = T1_compute_thresholds.or5_cases,
      output_dir = analysis_4_output_dir
  }
  call Tasks.copy_file_to_storage as copy5{
   input:
      text_file = T1_compute_thresholds.or5_controls,
      output_dir = analysis_4_output_dir
  }

}

task T1_compute_thresholds {
    input {
      File pgs_file
      File phenos_file
      String PGS_ID
      String cancer_type
    }

    command <<<
    set -euo pipefail

    python3 <<CODE
    import pandas as pd
    import numpy as np
    from scipy.stats import fisher_exact
    from scipy import stats

    # Load data
    pgs = pd.read_csv("~{pgs_file}", sep='\t')
    phenos = pd.read_csv("~{phenos_file}", sep='\t')

    # Ensure string match
    pgs['sample'] = pgs['sample'].astype(str)
    phenos['Sample'] = phenos['Sample'].astype(str)

    # Merge
    df = pd.merge(pgs, phenos, left_on='sample', right_on='Sample', how="left")

    def find_thresholds_vs_middle(pgs_col):
        # Drop missing values
        prs_values = df[pgs_col].dropna()

        # Calculate middle 20% range (40th–60th percentile)
        control_prs_values = df.loc[df['cancer'].str.lower() == "control", pgs_col].dropna()
        lower_ref = np.percentile(control_prs_values, 40)
        upper_ref = np.percentile(control_prs_values, 60)

        # Define reference group
        ref_group = df[(df[pgs_col] >= lower_ref) & (df[pgs_col] <= upper_ref)]

        thresholds = sorted(prs_values.unique(), reverse=True)

        best_or2 = None
        best_or5 = None
        thresh_or2 = None
        thresh_or5 = None

        rows = []
        for t in thresholds:
            # High group = above threshold
            high = df[df[pgs_col] > t]

            # Skip if overlap with reference group (to avoid same individuals)
            if high.empty or ref_group.empty:
                continue

            # Counts for 2x2 table
            a = sum(high['cancer'] != 'control')  # cases in high PRS
            b = sum(high['cancer'] == 'control')  # controls in high PRS
            c = sum(ref_group['cancer'] != 'control')  # cases in reference
            d = sum(ref_group['cancer'] == 'control')  # controls in reference
            if min(a, b, c, d) == 0:
                continue  # avoid divide by zero

            OR = (a / b) / (c / d)

            # Fisher's exact p-value
            table = [[a, b], [c, d]]
            _, pval = stats.fisher_exact(table)

            rows.append((t, OR, pval))

            # Track thresholds where OR is closest to 2 and 5
            if best_or2 is None or abs(OR - 2) < abs(best_or2 - 2):
                best_or2 = OR
                thresh_or2 = t

            if best_or5 is None or abs(OR - 5) < abs(best_or5 - 5):
                best_or5 = OR
                thresh_or5 = t

        # Save results
        df_out = pd.DataFrame(rows, columns=["threshold", "odds_ratio", "p_value"])
        df_out.to_csv("prs_thresholds_vs_middle.txt", sep="\t", index=False)

        return thresh_or2, thresh_or5


    thresholds_or2 = {}
    thresholds_or5 = {}
    samples_or2 = {}
    samples_or5 = {}

    pgs_cols = [col for col in df.columns if col.startswith("PGS")]

    for pgs_col in pgs_cols:
        t2, t5 = find_thresholds_vs_middle(pgs_col)
        thresholds_or2[pgs_col] = t2
        thresholds_or5[pgs_col] = t5
        samples_or2[pgs_col] = df[(df[pgs_col] > t2)]['sample'].tolist()
        samples_or5[pgs_col] = df[(df[pgs_col] > t5)]['sample'].tolist()

    # Save results
    pd.DataFrame.from_dict(thresholds_or2, orient='index', columns=['threshold']).to_csv("~{cancer_type}.~{PGS_ID}.thresholds_or2.tsv", sep='\t')
    pd.DataFrame.from_dict(thresholds_or5, orient='index', columns=['threshold']).to_csv("~{cancer_type}.~{PGS_ID}.thresholds_or5.tsv", sep='\t')

    def write_case_control_files(samples_dict, label, df, cancer_col="cancer"):
        cases_file = f"~{cancer_type}.~{PGS_ID}.{label}.cases.tsv"
        controls_file = f"~{cancer_type}.~{PGS_ID}.{label}.controls.tsv"

        with open(cases_file, 'w') as f_cases, open(controls_file, 'w') as f_controls:
            for score, ids in samples_dict.items():
                for i in ids:
                    status = df.loc[df['sample'] == i, cancer_col].values[0]
                    if status.lower() != "control":   # anything not 'control'
                        f_cases.write(f"{i}\n")
                    else:
                        f_controls.write(f"{i}\n")

    # Write separate files for OR2 and OR5
    write_case_control_files(samples_or2, "OR2", df)
    write_case_control_files(samples_or5, "OR5", df)

    CODE
    >>>

    output {
      File or2_thresholds = "~{cancer_type}.~{PGS_ID}.thresholds_or2.tsv"
      File or5_thresholds = "~{cancer_type}.~{PGS_ID}.thresholds_or5.tsv"
      File or2_cases = "~{cancer_type}.~{PGS_ID}.OR2.cases.tsv"
      File or2_controls = "~{cancer_type}.~{PGS_ID}.OR2.controls.tsv"
      File or5_cases = "~{cancer_type}.~{PGS_ID}.OR5.cases.tsv"
      File or5_controls = "~{cancer_type}.~{PGS_ID}.OR5.controls.tsv"
      File out1 = "prs_thresholds_vs_middle.txt"
    }

    runtime {
      docker: "vanallenlab/pydata_stack"
      memory: "4G"
    }
}
