# The Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0

struct CancerPair {
  String cancer_A
  String prs_id_A
  String cancer_B
  String prs_id_B
}


# TODO Remove Explained Samples
workflow ANALYSIS_4E_PRS {
  input {
    String google_workspace = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
    File metadata_tsv = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/dfci-ufc.aou.phenos.v2.tsv.gz"
    File explained_samples = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/samples_with_pvs.list"
    String metadata_sample_col = "sample_id"
    String metadata_cancer_col = "cancer"
    String prs_sample_col = "sample_id"
    String prs_score_col = "PGS"
    Float prs_threshold = 3.0
    Array[CancerPair] cancer_pairs
  }


  scatter(p in cancer_pairs){
    File prs_path_A = google_workspace + "/ANALYSIS_4_PRS/" + p.cancer_A + "." + p.prs_id_A + ".pgs"
    File prs_path_B = google_workspace + "/ANALYSIS_4_PRS/" + p.cancer_B + "." + p.prs_id_B + ".pgs"

    File prsA_threshold = google_workspace + "/ANALYSIS_4_PRS/" + p.cancer_A + "." + p.prs_id_A + ".thresholds_or5.tsv"
    File prsB_threshold = google_workspace + "/ANALYSIS_4_PRS/" + p.cancer_B + "." + p.prs_id_B + ".thresholds_or5.tsv"

    call T1_filter_metadata_for_both_cancers {
      input:
        metadata_tsv = metadata_tsv,
        explained_samples = explained_samples,
        cancer_A = p.cancer_A,
        cancer_B = p.cancer_B
    }

    call extract_prs_scores as extract_prs_scores_A {
      input:
        prs_tsv = prs_path_A,
        sample_list = T1_filter_metadata_for_both_cancers.out1,
        prs_sample_col = prs_sample_col,
        prs_score_col = prs_score_col
    }

    call extract_prs_scores as extract_prs_scores_B {
      input:
        prs_tsv = prs_path_B,
        sample_list = T1_filter_metadata_for_both_cancers.out1,
        prs_sample_col = prs_sample_col,
        prs_score_col = prs_score_col
      }
      call compute_fisher {
        input:
          prsA_scores_tsv = extract_prs_scores_A.out1,
          prsB_scores_tsv = extract_prs_scores_B.out1,
          prsA_threshold = prsA_threshold,
          prsB_threshold = prsB_threshold
      }
  }





  output {
    #File matched_samples = filter_metadata_for_both_cancers.matched_samples_tsv
    #File prsA_extracted = extract_prs_scores_A.scores_tsv
    #File prsB_extracted = extract_prs_scores_B.scores_tsv
    #File contingency_table = compute_fisher.contingency_tsv
    #Float fisher_pvalue = compute_fisher.fisher_pvalue
    #Float fisher_odds_ratio = compute_fisher.fisher_oddsratio
    #File fisher_report = compute_fisher.report_txt
  }
}


task T1_filter_metadata_for_both_cancers {
  input {
    File metadata_tsv
    File explained_samples
    String cancer_A
    String cancer_B
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  df = pd.read_csv("~{metadata_tsv}", sep="\t", index_col=False)

  # Load explained samples (assume one sample per line)
  explained = pd.read_csv("~{explained_samples}", header=None)[0].astype(str).tolist()

  # Select rows that contain both cancers in 'original_dx'
  mask = df['original_dx'].str.contains("~{cancer_A}", case=False, regex=False, na=False) & \
         df['original_dx'].str.contains("~{cancer_B}", case=False, regex=False, na=False)

  samples = df.loc[mask, 'Sample'].astype(str)

  # Remove explained samples
  filtered = samples[~samples.isin(explained)]

  # Save as a one-column list file
  filtered.to_csv("matched_samples.list", index=False, header=False)
  CODE
  >>>
  runtime { 
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
  output {
    File out1 = "matched_samples.list"
  }
}

task extract_prs_scores {
  input {
    File prs_tsv
    File sample_list
    String prs_sample_col
    String prs_score_col
  }
  command <<<
  python3 <<CODE
  import pandas as pd

  # Inputs
  prs = pd.read_csv("~{prs_tsv}", sep="\t")
  samples = pd.read_csv("~{sample_list}", header=None, names=["sample"], sep="\t")

  # Subset to the relevant samples
  subset = prs[prs["sample"].isin(samples["sample"])][["sample", "~{prs_score_col}"]]

  # Save
  subset.to_csv("extracted_scores.tsv", sep="\t", index=False)
  print(f"Wrote extracted_scores.tsv with {len(subset)} rows")

  CODE
  >>>
  runtime { 
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
  output {
    File out1 = "extracted_scores.tsv"
  }
}

task compute_fisher {
  input {
    File prsA_scores_tsv
    File prsB_scores_tsv
    File prsA_threshold
    File prsB_threshold
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  from scipy.stats import fisher_exact

  thresholdA_df = pd.read_csv("~{prsA_threshold}", sep="\t", index_col=0)
  thresholdA = df.loc["PGS", "threshold"]
  thresholdB_df = pd.read_csv("~{prsB_threshold}", sep="\t", index_col=0)
  thresholdB = df.loc["PGS", "threshold"]

  A = pd.read_csv("~{prsA_scores_tsv}", sep="\t", index_col=False)
  B = pd.read_csv("~{prsB_scores_tsv}", sep="\t", index_col=False)
  # Assume both files have sample_id and score columns (first two columns)
  A.columns = [c for c in A.columns]
  B.columns = [c for c in B.columns]
  # normalize columns: first is sample id, second is score
  A = A.iloc[:,0:2]
  B = B.iloc[:,0:2]
  A.columns = ["sample","scoreA"]
  B.columns = ["sample","scoreB"]
  # merge
  merged = pd.merge(A, B, on="sample", how="inner")
  merged["highA"] = merged["scoreA"] > ~{prsA_threshold}
  merged["highB"] = merged["scoreB"] > ~{prsB_threshold}
  # build contingency:
  #            highB   not_highB
  # highA       a         b
  # not_highA   c         d
  a = int(((merged["highA"]) & (merged["highB"])).sum())
  b = int(((merged["highA"]) & (~merged["highB"])).sum())
  c = int(((~merged["highA"]) & (merged["highB"])).sum())
  d = int(((~merged["highA"]) & (~merged["highB"])).sum())
  cont = pd.DataFrame([["highA_highB", a],
                       ["highA_notHighB", b],
                       ["notHighA_highB", c],
                       ["notHighA_notHighB", d]],
                      columns=["cell", "count"])
  cont.to_csv("contingency_long.tsv", sep="\\t", index=False)
  # fisher test
  table = [[a, b],[c, d]]
  oddsratio, pvalue = fisher_exact(table, alternative="two-sided")
  with open("fisher_report.txt","w") as fh:
      fh.write("n_samples_merged\\t{0}\\n".format(len(merged)))
      fh.write("contingency_table\\n")
      fh.write("        highB not_highB\\n")
      fh.write(f"highA   {a} {b}\\n")
      fh.write(f"notHighA {c} {d}\\n")
      fh.write(f"odds_ratio\\t{oddsratio}\\n")
      fh.write(f"pvalue\\t{pvalue}\\n")
  # also write wide contingency as TSV
  import csv
  with open("contingency_wide.tsv","w") as f:
      w = csv.writer(f, delimiter="\\t")
      w.writerow(["", "highB", "not_highB"])
      w.writerow(["highA", a, b])
      w.writerow(["notHighA", c, d])
  print("WROTE contingency and fisher_report")
  # expose values via stdout for WDL to capture
  print(pvalue)
  print(oddsratio)
  CODE
  >>>
  runtime { 
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
  output {
    File contingency_tsv = "contingency_wide.tsv"
    File report_txt = "fisher_report.txt"
    #Float fisher_pvalue = read_float(stdout(), 1)   # first printed float
    #Float fisher_oddsratio = read_float(stdout(), 2) # second printed float
  }
}
