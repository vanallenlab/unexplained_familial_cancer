version 1.0

workflow identify_pgs_thresholds {
  input {
    String PGS_ID = ""
    String google_bucket = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
    String cancer_type
    File phenos_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/dfci-ufc.aou.phenos.v2.tsv.gz"
  }
  File pgs_file = google_bucket + "/ANALYSIS_4_PRS/" + cancer_type + "." + PGS_ID + ".pgs"

  call T1_compute_thresholds {
    input:
      pgs_file = pgs_file,
      phenos_file = phenos_file
  }

  output {
    File or2_thresholds = compute_thresholds.or2_thresholds
    File or5_thresholds = compute_thresholds.or5_thresholds
    File or2_samples = compute_thresholds.or2_samples
    File or5_samples = compute_thresholds.or5_samples
  }
}

task T1_compute_thresholds {
    input {
      File pgs_file
      File phenos_file
    }

    command <<<
    set -euo pipefail

    python3 <<EOF
    import pandas as pd
    import numpy as np
    from scipy.stats import fisher_exact

    # Load data
    pgs = pd.read_csv("~{pgs_file}", sep='\t')
    phenos = pd.read_csv("~{phenos_file}", sep='\t')

    # Ensure string match
    pgs['Sample'] = pgs['Sample'].astype(str)
    phenos['original_id'] = phenos['original_id'].astype(str)

    # Merge
    df = pd.merge(pgs, phenos, left_on='Sample', right_on='original_id')

    # Only keep if cancer is 'case' or 'control'
    df = df[df['cancer'].isin(['case', 'control'])]

    def find_thresholds(pgs_col):
        thresholds = sorted(df[pgs_col].dropna().unique())
        best_or2 = None
        best_or5 = None
        thresh_or2 = None
        thresh_or5 = None

        for t in thresholds:
            high = df[df[pgs_col] > t]
            low = df[df[pgs_col] <= t]

            a = sum(high['cancer'] == 'case')
            b = sum(high['cancer'] == 'control')
            c = sum(low['cancer'] == 'case')
            d = sum(low['cancer'] == 'control')

            if min(a, b, c, d) == 0:
                continue  # avoid divide by zero

            OR = (a / b) / (c / d)

            if best_or2 is None or abs(OR - 2) < abs(best_or2 - 2):
                best_or2 = OR
                thresh_or2 = t

            if best_or5 is None or abs(OR - 5) < abs(best_or5 - 5):
                best_or5 = OR
                thresh_or5 = t

        return thresh_or2, thresh_or5

    thresholds_or2 = {}
    thresholds_or5 = {}
    samples_or2 = {}
    samples_or5 = {}

    pgs_cols = [col for col in df.columns if col.startswith("PGS")]

    for pgs_col in pgs_cols:
        t2, t5 = find_thresholds(pgs_col)
        thresholds_or2[pgs_col] = t2
        thresholds_or5[pgs_col] = t5
        samples_or2[pgs_col] = df[(df[pgs_col] > t2) & (df['phenotype'] != 'control')]['Sample'].tolist()
        samples_or5[pgs_col] = df[(df[pgs_col] > t5) & (df['phenotype'] != 'control')]['Sample'].tolist()

    # Save results
    pd.DataFrame.from_dict(thresholds_or2, orient='index', columns=['threshold']).to_csv("thresholds_or2.tsv", sep='\t')
    pd.DataFrame.from_dict(thresholds_or5, orient='index', columns=['threshold']).to_csv("thresholds_or5.tsv", sep='\t')

    with open("samples_or2.tsv", 'w') as f:
        for score, ids in samples_or2.items():
            for i in ids:
                f.write(f"{score}\t{i}\n")

    with open("samples_or5.tsv", 'w') as f:
        for score, ids in samples_or5.items():
            for i in ids:
                f.write(f"{score}\t{i}\n")
  EOF
  >>>

  output {
    File or2_thresholds = "thresholds_or2.tsv"
    File or5_thresholds = "thresholds_or5.tsv"
    File or2_samples = "samples_or2.tsv"
    File or5_samples = "samples_or5.tsv"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "4G"
    cpu: 1
  }
}
