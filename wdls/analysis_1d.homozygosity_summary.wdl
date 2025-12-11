# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0
version 1.0

workflow ANALYSIS_1D_Homozygosity_Summary {
  input {
    String cancer_type
    String haplotype
    File vcf
  }
  File metadata = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + '.metadata'
  File vcf_idx = vcf + ".tbi"

  call T1_Prepare_VCF {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      haplotype = haplotype
  }
  call T2_Prepare_Tsv {
    input:
      transposed_matrix = T1_Prepare_VCF.out1,
      metadata = metadata,
      output_name = cancer_type + "." + haplotype + "snp_homozygosity.tsv"
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
  CODE
  >>>
  output {
    File out1 = "~{output_name}"
  }
  runtime {
    docker:"vanallenlab/pydata_stack"
    preemptible:3
  }
}

