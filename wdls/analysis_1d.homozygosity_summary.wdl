# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0
version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_1D_Homozygosity_Summary {
  input {
    String cancer_type
    String haplotype
    String extra_haplotype
    String chr
    Array[File] vcfs
    Array[File] vcf_idxs
    File analysis_1b_output = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/analysis_1b_output.tsv.gz"
  }
  File metadata = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + '.metadata'
  File roh_tsv = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr" + chr + ".txt"

  call Tasks.ConcatVcfs {
    input:
      vcfs = vcfs,
      vcf_idxs = vcf_idxs,
      callset_name = "out"
  }
  call T1_Prepare_VCF {
    input:
      vcf = ConcatVcfs.merged_vcf,
      vcf_idx = ConcatVcfs.merged_vcf_idx,
      haplotype = extra_haplotype
  }
  call T2_Prepare_Tsv {
    input:
      transposed_matrix = T1_Prepare_VCF.out1,
      metadata = metadata,
      output_name = cancer_type + "." + haplotype + ".snp_homozygosity.1d.tsv"
  }
  call T3_Prepare_AUC {
    input:
      analysis_1b_output = analysis_1b_output,
      metadata = metadata,
      haplotype = haplotype
  }
  call T4_locus_of_ROH {
    input:
      roh_tsv = roh_tsv,
      extra_haplotype = extra_haplotype,
      metadata = metadata,
      output_name = cancer_type + "." + haplotype + ".roh_homozygosity.1d.tsv"
  }
}

task T3_Prepare_AUC {
  input {
    File analysis_1b_output
    File metadata
    String haplotype
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  from sklearn.metrics import roc_auc_score
  metadata = pd.read_csv("~{metadata}",sep='\t',index_col=False)
  hap_matrix = pd.read_csv("~{analysis_1b_output}",sep='\t',index_col=False)
  haplotype = "~{haplotype}".replace(':','_').replace('-','_')  #chr8:105100001-105600000
  print(haplotype)
  hap_matrix = hap_matrix[hap_matrix["GeneName"] == haplotype]

  # Drop the first two columns: Genomic_Region and GeneName
  hap_matrix = hap_matrix.iloc[:, 2:]
  hap_matrix = hap_matrix.T
  print(hap_matrix)

  hap_matrix = (
      hap_matrix
      .reset_index()   # index → column
      .rename(columns={
          "index": "original_id",
          hap_matrix.columns[0]: "percent_homozygosed"
      })
  )

  # Ensure original_id is string
  hap_matrix["original_id"] = hap_matrix["original_id"].astype(str)

  # Ensure metadata original_id is also string
  metadata["original_id"] = metadata["original_id"].astype(str)
  merged = metadata.merge(
    hap_matrix,
    on="original_id",
    how="left"
  )

  merged["cancer_status"] = (merged["original_dx"] != "control").astype(int)
  roc_df = merged[["original_id", "cancer_status", "percent_homozygosed"]].dropna()
  
  auroc = roc_auc_score(
    roc_df["cancer_status"],
    roc_df["percent_homozygosed"]
  )
  print(f"AUROC (percent_homozygosed vs cancer_status): {auroc:.4f}")
  roc_df.to_csv(
      "~{haplotype}.auroc.tsv",
      sep="\t",
      index=False
  )
  CODE   
  >>>
  output {
    #File out1 = "~{haplotype}.auroc.tsv"
  }
  runtime {
    memory: "16GB"
    disks: "local-disk 10 HDD"
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task T1_Prepare_VCF {
  input {
    File vcf
    File vcf_idx
    String haplotype
  }

  command <<<
  set -euxo pipefail
  bcftools view -r ~{haplotype} ~{vcf} -Oz -o tmp1.vcf.gz
  bcftools view -i 'AF>0.01 & AF<0.99' tmp1.vcf.gz -Oz -o tmp2.vcf.gz

  # 1. extract sample names
  samples=$(bcftools query -l tmp2.vcf.gz | tr '\n' '\t')

  # 2. write header
  echo -e "CHROM\tPOS\tREF\tALT\t${samples}" > genotype_matrix.tsv

  # 3. append variants with converted GT values
  bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
    tmp2.vcf.gz |
  awk 'BEGIN{FS=OFS="\t"} {
    for (i=5; i<=NF; i++) {
      if ($i ~ /^[0][\/|][0]$/) $i=0;
      else if ($i ~ /^[0][\/|][1]$/ || $i ~ /^[1][\/|][0]$/) $i=1;
      else if ($i ~ /^[1][\/|][1]$/) $i=2;
      else $i=".";
    }
    print
  }' >> genotype_matrix.tsv

  python3 <<CODE
  import pandas as pd

  # Load genotype matrix
  df = pd.read_csv("genotype_matrix.tsv", sep="\t")

  # Create merged variant ID: chr_pos_ref_alt
  df["variant_id"] = df["CHROM"].astype(str) + "_" + \
                    df["POS"].astype(str) + "_" + \
                    df["REF"] + "_" + df["ALT"]

  # Drop the original columns
  df = df.drop(columns=["CHROM", "POS", "REF", "ALT"])

  # Move variant_id to index
  df = df.set_index("variant_id")

  # Transpose → samples become rows
  df_t = df.T

  # Save
  df_t.to_csv("genotype_matrix_transposed.tsv", sep="\t")
  CODE
  >>>

  output {
    File out1 = "genotype_matrix_transposed.tsv"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
    memory: "4G"
  }
}

task T2_Prepare_Tsv {
  input {
    File transposed_matrix
    File metadata
    String output_name
  }
  command <<<
  sed '1s/^/original_id/' ~{transposed_matrix} > tmp.tsv
  python3 <<CODE
  import pandas as pd
  import numpy as np

  # ------------------------------------------------------
  # Load data
  # ------------------------------------------------------
  geno = pd.read_csv("tmp.tsv", sep="\t", dtype=str)
  meta = pd.read_csv("~{metadata}", sep="\t", dtype=str)

  # Make IDs strings
  geno["original_id"] = geno["original_id"].astype(str)
  meta["original_id"] = meta["original_id"].astype(str)

  # Merge metadata + genotypes
  df = meta.merge(geno, on="original_id", how="left")

  # ------------------------------------------------------
  # Identify case / control samples
  # ------------------------------------------------------
  df["is_control"] = df["original_dx"].str.lower().eq("control")
  df["is_case"] = ~df["is_control"]

  # ------------------------------------------------------
  # Identify SNP columns (all genotype columns except metadata)
  # ------------------------------------------------------
  snp_cols = [c for c in df.columns if c.startswith("chr") and not c.endswith("_ploidy")]
  meta_cols = [c for c in df.columns if c not in snp_cols]

  # ------------------------------------------------------
  # Function: percent homozygous at a SNP
  # ------------------------------------------------------
  def percent_homozygous(series):
      """Return % homozygous (0 or 2) ignoring missing values."""
      vals = series.replace(".", np.nan).astype(float)
      vals = vals.dropna()
      if len(vals) == 0:
          return np.nan
      hom = ((vals == 0) | (vals == 2)).sum()
      return hom / len(vals)

  # ------------------------------------------------------
  # Compute % homozygous for every SNP
  # ------------------------------------------------------
  results = []

  for snp in snp_cols:
      case_vals = df.loc[df["is_case"], snp]
      ctrl_vals = df.loc[df["is_control"], snp]

      case_pct = percent_homozygous(case_vals)
      ctrl_pct = percent_homozygous(ctrl_vals)

      results.append([snp, case_pct, ctrl_pct])

  # ------------------------------------------------------
  # Create results dataframe
  # ------------------------------------------------------
  result_df = pd.DataFrame(
      results,
      columns=["snp", "case_percent_homozygosed", "control_percent_homozygosed"]
  )

  # Save output
  result_df.to_csv("~{output_name}", sep="\t", index=False)
  gzip ~{output_name}
  CODE
  >>>
  output {
    File out1 = "~{output_name}.gz"
  }
  runtime {
    docker:"vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task T4_locus_of_ROH {
  input {
    File roh_tsv
    String extra_haplotype
    File metadata
    String output_name
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  import numpy as np
  # ------------------------------------------------------------
  # 1. Load ROH file
  #    Each row = one ROH segment for one sample
  # ------------------------------------------------------------
  roh_df = pd.read_csv(
      "~{roh_tsv}",
      sep="\t",
      index_col=False,
      header=None,
      comment="#",
      names=["RG", "original_id", "chr", "start", "end", "length", "num_markers", "Quality"]
  )

  # Keep only relevant columns
  roh_df = roh_df[["original_id", "chr", "start", "end", "length"]]

  # ------------------------------------------------------------
  # 2. Load metadata
  # ------------------------------------------------------------
  metadata = pd.read_csv("~{metadata}", sep="\t", index_col=False)

  # Restrict ROHs to samples present in metadata
  sample_list = metadata["original_id"].astype(str).tolist()
  roh_df["original_id"] = roh_df["original_id"].astype(str)
  roh_df = roh_df[roh_df["original_id"].isin(sample_list)]

  # Define cases and controls
  cases_list = metadata[metadata["original_dx"] != "control"]["original_id"].astype(str).tolist()
  controls_list = metadata[metadata["original_dx"] == "control"]["original_id"].astype(str).tolist()

  num_cases = len(cases_list)
  num_controls = len(controls_list)

  # ------------------------------------------------------------
  # 3. Parse haplotype coordinates
  #    Format assumed: chrX.start-end
  # ------------------------------------------------------------
  haplotype_start = int("~{extra_haplotype}".split(":")[1].split("-")[0])
  haplotype_end   = int("~{extra_haplotype}".split(":")[1].split("-")[1])

  # Restrict ROHs to those overlapping the haplotype
  roh_df = roh_df[
      (roh_df["start"] <= haplotype_end) &
      (roh_df["end"] >= haplotype_start)
  ].copy()

  # ------------------------------------------------------------
  # 4. Initialize per-base-pair signal arrays
  #    One value per base in the haplotype window
  # ------------------------------------------------------------
  hap_length = haplotype_end - haplotype_start + 1

  case_signal = [0.0] * hap_length
  control_signal = [0.0] * hap_length

  # Normalized increment per ROH
  case_increment = 1.0 / num_cases if num_cases > 0 else 0.0
  control_increment = 1.0 / num_controls if num_controls > 0 else 0.0

  # ------------------------------------------------------------
  # 5. Add ROH signal per base pair
  #    Overlapping ROHs stack additively
  # ------------------------------------------------------------
  for _, row in roh_df.iterrows():

      # Clip ROH to haplotype boundaries
      start = max(row["start"], haplotype_start)
      end   = min(row["end"], haplotype_end)

      if start > end:
          continue

      # Convert genomic positions to array indices
      start_idx = start - haplotype_start
      end_idx   = end - haplotype_start

      # Add normalized signal
      if row["original_id"] in cases_list:
          for i in range(start_idx, end_idx + 1):
              case_signal[i] += case_increment

      elif row["original_id"] in controls_list:
          for i in range(start_idx, end_idx + 1):
              control_signal[i] += control_increment

  # ------------------------------------------------------------
  # 6. Build output dataframe
  # ------------------------------------------------------------
  positions = list(range(haplotype_start, haplotype_end + 1))

  signal_df = pd.DataFrame({
      "chr": roh_df.iloc[0]["chr"] if len(roh_df) > 0 else None,
      "pos": positions,
      "case_roh_fraction": np.round(case_signal,2),
      "control_roh_fraction": np.round(control_signal,2)
  })

  # ------------------------------------------------------------
  # 7. Write output file (ready for plotting)
  # ------------------------------------------------------------
  signal_df.to_csv(
      "~{output_name}",
      sep="\t",
      index=False
  )

  CODE
  gzip ~{output_name}
  >>>
  output {
    File out1 = "~{output_name}.gz"
  }
  runtime {
    preemptible:3
    docker:"vanallenlab/g2c_pipeline"
  }
}

