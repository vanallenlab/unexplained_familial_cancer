# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g
# Adapted from Ryan Collins code:
# https://github.com/vanallenlab/ped_germline_SV/blob/main/gatksv_scripts/get_kinship_and_pcs.py`

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow EXTRA_1_MAYBE_PGV {
  input {
    File sv_tsv = "gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/aacr_sv_files/ufc.prelim_plp_rare_svs.autosomal.remapped_ids.tsv.gz"
    Array[String] cancer_types = ["basal_cell","breast","colorectal","thyroid","bladder","uterus","ovary","hematologic","non-hodgkin","melanoma","prostate","squamous_cell","kidney","lung","sarcoma","neuroendocrine","brain"] 
  }
  scatter (cancer_type in cancer_types){
    call T0_find_samples {
      input:
        step_10_output = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/FINAL_FILES/step_10_output.dec10_2025.tsv.gz", 
        cosmic_tsv = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cosmic_ufc.v3.tsv",
        cancer_type = cancer_type,
        sv_tsv = sv_tsv
    }
    call Tasks.copy_file_to_storage {
      input:
        text_file = T0_find_samples.out1,
        output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/maybe_explained_samples/"
    }
  }
}

task T0_find_samples {
  input {
    File step_10_output
    File cosmic_tsv
    File sv_tsv
    File already_explained = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/samples_with_pvs.nov7.list"
    String cancer_type
  }
  command <<<
  python3 <<CODE
  import pandas as pd

  # ------------------------------------------------------------
  # 1. Load main variant matrix
  # ------------------------------------------------------------
  df = pd.read_csv("~{step_10_output}", sep="\t", index_col=False)
  sv_df = pd.read_csv("~{sv_tsv}",sep='\t',index_col=False)
  # ------------------------------------------------------------
  # 2. Keep only desired Tier impacts
  # ------------------------------------------------------------
  valid_tiers = ("Tier3_0001")
  df = df[
      df["gene_impact"]
      .astype(str)
      .str.endswith(valid_tiers)
  ].copy()

  # ------------------------------------------------------------
  # 3. Load COSMIC gene list and filter by cancer type
  # ------------------------------------------------------------
  cosmic_df = pd.read_csv(
      "~{cosmic_tsv}",
      sep="\t",
      header=None,
      names=["gene", "cancers"]
  )
  cosmic_df = cosmic_df[
      cosmic_df["cancers"]
      .astype(str)
      .str.contains("~{cancer_type}", case=False, na=False)
  ]

  # Extract unique COSMIC genes
  gene_set = set(
      cosmic_df["gene"]
      .dropna()
      .astype(str)
      .unique()
  )

  # ------------------------------------------------------------
  # 4. Extract gene name from gene_impact column
  #    (gene is the first token before "_")
  # ------------------------------------------------------------
  df["gene"] = (
      df["gene_impact"]
      .astype(str)
      .str.split("_")
      .str[0]
  )

  # ------------------------------------------------------------
  # 5. Filter variant matrix to COSMIC genes only
  # ------------------------------------------------------------
  df_filtered = df[df["gene"].isin(gene_set)].copy()
  print(df_filtered['gene_impact'].tolist())
  # ------------------------------------------------------------
  # 6. Identify columns with ≥1 signal
  #     (numeric columns only)
  # ------------------------------------------------------------
  numeric_df = df_filtered.select_dtypes(include="number")

  cols_with_signal = numeric_df.columns[
      (numeric_df >= 1).any()
  ].tolist()
  numeric_df.to_csv("halfway.tsv",sep='\t',index=False)
  # ------------------------------------------------------------
  # 8. Incorporate structural variants (SVs)
  #    If any gene in sv_df["genes"] overlaps gene_set,
  #    add the corresponding #sample to cols_with_signal
  # ------------------------------------------------------------

  # Ensure cols_with_signal is a set to avoid duplicates
  cols_with_signal = set(cols_with_signal)

  for _, row in sv_df.iterrows():
      if pd.isna(row["genes"]):
          continue

      # Split comma-delimited gene list and strip whitespace
      sv_genes = {
          g.strip()
          for g in str(row["genes"]).split(",")
          if g.strip()
      }

      # Check for overlap with COSMIC gene set
      if sv_genes & gene_set:
          cols_with_signal.add(str(row["#sample"]))

  # Convert back to sorted list if needed
  cols_with_signal = sorted(cols_with_signal)

  # ------------------------------------------------------------
  # 7. Write column names with signal to file
  # ------------------------------------------------------------
  with open("tmp.list", "w") as out:
      for col in cols_with_signal:
          out.write(f"{col}\n")
  CODE
  cat ~{already_explained} tmp.list | sort -u > ~{cancer_type}.sensitivity_ppv.list
  >>>
  output {
    File out1 = "~{cancer_type}.sensitivity_ppv.list"
    File out2 = "halfway.tsv"
  }
  runtime {
    docker:"vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task T1_Convert_To_TSV {
  input {
    File vcf
    File cpg_list
    File samples_of_interest
    File cpg_bed_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/riaz_genes.coordinates.bed"
    File rare_variants_001 = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/rare_001.tsv"
  }
  String output_file = basename(vcf, ".vcf.bgz") + ".tsv"
  command <<<
  set -euxo pipefail
  
  # Filter to just Samples of interest
  bcftools view -S ~{samples_of_interest} -Oz -o tmp1.vcf.gz ~{vcf}
  bcftools index -t tmp1.vcf.gz
  bcftools view -R ~{cpg_bed_file} tmp1.vcf.gz -Oz -o tmp2.vcf.gz

  bcftools view --include ID==@~{rare_variants_001} tmp2.vcf.gz -O z -o tmp3.vcf.gz
  rm ~{rare_variants_001}
 
  # Get the header started
  echo -e "CHROM\tPOS\tID\tREF\tALT\tIMPACT\tSYMBOL\tclinvar_clnsig\tBiotype\tConsequence\tSAMPLES" > ~{output_file}


  bcftools +split-vep tmp3.vcf.gz -i 'GT="alt"' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%IMPACT\t%SYMBOL\t%clinvar_CLNSIG\t%BIOTYPE\t%Consequence\t[%SAMPLE,]\n' -d > tmp.txt
  rm tmp1.vcf.gz

  # Filter the File to only include CPGS
  awk 'NR==FNR {cpg[$1]; next} $7 in cpg' <(cut -f1 ~{cpg_list}) tmp.txt > tmp3.txt || touch tmp3.txt

  sort -u < tmp3.txt >> ~{output_file}

  python3 <<CODE
  import pandas as pd
  df = pd.read_csv("~{output_file}",sep='\t',index_col=False)
  df = df[
      (
          (df["clinvar_clnsig"].isin(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"]))
      )
      |
      (
          (df["IMPACT"] == "HIGH") &
          (df["clinvar_clnsig"].isin(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic", "."]))
      )
  ]

  #df = df[(df['IMPACT'] == "MODERATE") | (df['IMPACT'] == "HIGH")]
  #df = df[(df['clinvar_clnsig'] == "Pathogenic") | (df['clinvar_clnsig'] == "Likely_pathogenic") | (df['clinvar_clnsig'] == "Pathogenic/Likely_pathogenic") | (df['clinvar_clnsig'] == ".")]
  df.to_csv("~{output_file}",sep='\t',index=False)
  CODE
  >>>
  output {
    File out1 = "~{output_file}"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    disks: "local-disk 20 HDD"
    memory: "8GB"
    preemptible: 1
  }
}

task concatenateFiles {
  input {
    Array[File] files
    String callset_name
  }

  command <<<
  # Put Header Down
  head -n -1 ~{files[0]} > ~{callset_name}.tsv

  # Add files
  for f in ~{sep=" " files}; do
    tail -n +2 "$f"
  done > ~{callset_name}.tsv
  >>>

  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }

  output {
    File out = "~{callset_name}.tsv"
  }
}

