# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_4B_PGS_INHERITANCE {
  input {
    File fam_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/fam_files/modified_ceph.fam"  # .fam file (6-column PLINK format)
    File analysis_4_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/"
  }
  String cohort = "ceph"

  call Tasks.list_files_from_directory {
    input:
      dir = analysis_4_dir,
      suffix = ".raw.pgs"
  }

  scatter(i in range(length(list_files_from_directory.out1) - 25)){
    call T1_analyze_inheritance {
      input:
        fam_file = fam_file,
        pgs_file = list_files_from_directory.out1[i],
        cohort = cohort
    }
  }

  call Tasks.concatenateFiles {
    input:
      files = T1_analyze_inheritance.out1,
      output_name = "CEPH_transmissions"
  }

  call Tasks.copy_file_to_storage {
    input:
      text_file = concatenateFiles.out1,
      output_dir = analysis_4_dir
  }
}

task T1_analyze_inheritance {
  input {
    File fam_file
    File pgs_file
    String cohort
  }
  String output_file = basename(pgs_file, ".raw.pgs") + ".transmission"
  command <<<
  set -euxo pipefail

  python3 <<CODE
  import pandas as pd
  import numpy as np
  from scipy.stats import ttest_1samp

  # Load .fam and PGS files
  fam = pd.read_csv("~{fam_file}", delim_whitespace=True)
  pgs = pd.read_csv("~{pgs_file}", sep='\t')

  # Merge PGS info
  fam = fam.merge(pgs, left_on="SUBJECT_ID", right_on="sample", how="left")

  # Identify complete trios
  trios = fam.dropna(subset=["FATHER", "MOTHER"])
  trios = trios[trios["FATHER"] != "0"]
  trios = trios[trios["MOTHER"] != "0"]

  # Get proband + parent PGS
  pgs_dict = dict(zip(pgs['sample'], pgs['PGS']))
  records = []

  for _, row in trios.iterrows():
      child = row['SUBJECT_ID']
      father = row['FATHER']
      mother = row['MOTHER']
    
      if child in pgs_dict and father in pgs_dict and mother in pgs_dict:
          pgs_child = pgs_dict[child]
          pgs_parents_mean = (pgs_dict[father] + pgs_dict[mother]) / 2
          delta = pgs_child - pgs_parents_mean
          records.append((child, pgs_child, pgs_dict[father], pgs_dict[mother], delta))

  # Analyze
  df = pd.DataFrame(records, columns=["proband", "pgs_proband", "pgs_father", "pgs_mother", "delta"])
  n_trios = len(df)
  mean_delta = df['delta'].mean()
  std_delta = df['delta'].std()
  t_stat, p_value = ttest_1samp(df['delta'], 0)

  # Write summary
  with open("~{output_file}", "w") as f:
      f.write(f"PGS Inheritance Analysis\\n")
      f.write(f"Number of complete trios analyzed: {n_trios}\\n")
      f.write(f"Mean delta (proband - midparent PGS): {mean_delta:.4f}\\n")
      f.write(f"Std dev of delta: {std_delta:.4f}\\n")
      f.write(f"T-test against 0: t = {t_stat:.4f}, p = {p_value:.4e}\\n")
  CODE
  >>>

  output {
    File out1 = "~{output_file}"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "2G"
    preemptible: 3
    disk: "local-disk 4 HDD"
  }
}

