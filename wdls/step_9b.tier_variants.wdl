# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# # Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow STEP_9B_TIER_VARIANTS {
  input {
    String step_9_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/sharded_vcfs"
    File subjects_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cohorts/ufc_subjects.aou.list"
    String step_9b_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/"
  }
  
  Int positive_shards = 2960
  call Tasks.gather_positive_vcfs {
    input:
      dir = step_9_output_dir,
      positive_shards = positive_shards
  }

  call Tasks.sort_vcf_list {
    input:
      unsorted_vcf_list = gather_positive_vcfs.vcf_list
  }

  Int negative_shards = 0
  scatter (i in range(length(sort_vcf_list.vcf_arr)-negative_shards)){
    call T0_Initial_Filter{
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        subjects_list = subjects_list
    }

    #call T0_Filter_Tier0 as T0_Filter_Tier0_001{
    #  input:
    #    rare_variants = T0_Initial_Filter.out1,
    #    vcf = sort_vcf_list.vcf_arr[i]
    #}
    #call T0_Filter_Tier0 as T0_Filter_Tier0_0001{
    #  input:
    #    rare_variants = T0_Initial_Filter.out2,
    #    vcf = sort_vcf_list.vcf_arr[i]
    #}
    #call T1_Filter_Tier1 as T1_Filter_Tier1_001{
    #  input:
    #    rare_variants = T0_Initial_Filter.out1,
    #    vcf = sort_vcf_list.vcf_arr[i]
    #}
    #call T1_Filter_Tier1 as T1_Filter_Tier1_0001{
    #  input:
    #    rare_variants = T0_Initial_Filter.out2,
    #    vcf = sort_vcf_list.vcf_arr[i]
    #}
    #call T2_Filter_Tier2 as T2_Filter_Tier2_001{
    #  input:
    #    rare_variants = T0_Initial_Filter.splice_out1,
    #    vcf = sort_vcf_list.vcf_arr[i]
    #}
    #call T2_Filter_Tier2 as T2_Filter_Tier2_0001{
    #  input:
    #    rare_variants = T0_Initial_Filter.splice_out2,
    #    vcf = sort_vcf_list.vcf_arr[i]
    #}
    call T34_Filter_Tier34 as T34_Filter_Tier34_001{
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        rare_variants = T0_Initial_Filter.out1,
    }
    call T34_Filter_Tier34 as T34_Filter_Tier34_0001{
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        rare_variants = T0_Initial_Filter.out2,
    }
    #call T6_Filter_Tier6 as T6_Filter_Tier6_001{
    #  input:
    #    vcf = sort_vcf_list.vcf_arr[i],
    #    rare_variants = T0_Initial_Filter.out1,
    #}
    #call T6_Filter_Tier6 as T6_Filter_Tier6_0001{
    #  input:
    #    vcf = sort_vcf_list.vcf_arr[i],
    #    rare_variants = T0_Initial_Filter.out2,
    #}

  }

  call Tasks.concatenateFiles as Concat_Rare_001 {
    input:
      files = T0_Initial_Filter.out1,
      output_name = "rare_001"
  }
  call Tasks.concatenateFiles as Concat_Rare_0001 {
    input:
      files = T0_Initial_Filter.out2,
      output_name = "rare_0001"
  }

  call Tasks.copy_file_to_storage as copy_rare_001{
    input:
      text_file = Concat_Rare_001.out2,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy_rare_0001{
    input:
      text_file = Concat_Rare_0001.out2,
      output_dir = step_9b_output_dir
  }

  #call Tasks.concatenateFiles as Tier0_Concat_001 {
  #  input:
  #    files = T0_Filter_Tier0_001.out1,
  #    output_name = "tier0_001"
  #}
  #call Tasks.concatenateFiles as Tier0_Concat_0001 {
  #  input:
  #    files = T0_Filter_Tier0_0001.out1,
  #    output_name = "tier0_0001"
  #}
  #call Tasks.copy_file_to_storage as copy0_001{
  #  input:
  #    text_file = Tier0_Concat_001.out2,
  #    output_dir = step_9b_output_dir
  #}
  #call Tasks.copy_file_to_storage as copy0_0001{
  #  input:
  #    text_file = Tier0_Concat_0001.out2,
  #    output_dir = step_9b_output_dir
  #}
  #call Tasks.concatenateFiles as Tier1_Concat_001 {
  #  input:
  #    files = T1_Filter_Tier1_001.out1,
  #    output_name = "tier1_001"
  #}
  #call Tasks.concatenateFiles as Tier1_Concat_0001 {
  #  input:
  #    files = T1_Filter_Tier1_0001.out1,
  #    output_name = "tier1_0001"
  #}
  #call Tasks.copy_file_to_storage as copy1_001{
  #  input:
  #    text_file = Tier1_Concat_001.out2,
  #    output_dir = step_9b_output_dir
  #}
  #call Tasks.copy_file_to_storage as copy1_0001{
  #  input:
  #    text_file = Tier1_Concat_0001.out2,
  #    output_dir = step_9b_output_dir
  #}

  #call Tasks.concatenateFiles as Tier2_Concat_001 {
  #  input:
  #    files = T2_Filter_Tier2_001.out1,
  #    output_name = "tier2_001"
  #}
  #call Tasks.concatenateFiles as Tier2_Concat_0001 {
  #  input:
  #    files = T2_Filter_Tier2_0001.out1,
  #    output_name = "tier2_0001"
  #}
  #call Tasks.copy_file_to_storage as copy2_001{
  #  input:
  #    text_file = Tier2_Concat_001.out2,
  #    output_dir = step_9b_output_dir
  #}
  #call Tasks.copy_file_to_storage as copy2_0001{
  #  input:
  #    text_file = Tier2_Concat_0001.out2,
  #    output_dir = step_9b_output_dir
  #}

  call Tasks.concatenateFiles as Tier3_Concat_001 {
    input:
      files = T34_Filter_Tier34_001.out1,
      output_name = "tier3_001"
  }

  call Tasks.concatenateFiles as Tier3_Concat_0001 {
    input:
      files = T34_Filter_Tier34_0001.out1,
      output_name = "tier3_0001"
  }

  call Tasks.copy_file_to_storage as copy3_001{
    input:
      text_file = Tier3_Concat_001.out2,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy3_0001{
    input:
      text_file = Tier3_Concat_0001.out2,
      output_dir = step_9b_output_dir
  }
  call Tasks.concatenateFiles as Tier4_Concat_001 {
    input:
      files = T34_Filter_Tier34_001.out2,
      output_name = "tier4_001"
  }
  call Tasks.concatenateFiles as Tier4_Concat_0001 {
    input:
      files = T34_Filter_Tier34_0001.out2,
      output_name = "tier4_0001"
  }
  call Tasks.copy_file_to_storage as copy4_001{
    input:
      text_file = Tier4_Concat_001.out2,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy4_0001{
    input:
      text_file = Tier4_Concat_0001.out2,
      output_dir = step_9b_output_dir
  }
  #call Tasks.concatenateFiles as Tier5_Concat_001 {
  #  input:
  #    files = T34_Filter_Tier34_001.out3,
  #    output_name = "tier5_001"
  #}
  #call Tasks.concatenateFiles as Tier5_Concat_0001 {
  #  input:
  #    files = T34_Filter_Tier34_0001.out3,
  #    output_name = "tier5_0001"
  #}
  #call Tasks.copy_file_to_storage as copy5_001{
  #  input:
  #    text_file = Tier5_Concat_001.out2,
  #    output_dir = step_9b_output_dir
  #}
  #call Tasks.copy_file_to_storage as copy5_0001{
  #  input:
  #    text_file = Tier5_Concat_0001.out2,
  #    output_dir = step_9b_output_dir
  #}

  #call Tasks.concatenateFiles as Tier6_Concat_001 {
  #  input:
  #    files = T6_Filter_Tier6_001.out1,
  #    output_name = "tier6_001"
  #}
  #call Tasks.concatenateFiles as Tier6_Concat_0001 {
  #  input:
  #    files = T6_Filter_Tier6_0001.out1,
  #    output_name = "tier6_0001"
  #}

  #call Tasks.copy_file_to_storage as copy6_001{
  #  input:
  #    text_file = Tier6_Concat_001.out2,
  #    output_dir = step_9b_output_dir
  #}
  #call Tasks.copy_file_to_storage as copy6_0001{
  #  input:
  #    text_file = Tier6_Concat_0001.out2,
  #    output_dir = step_9b_output_dir
  #}


  #call Tasks.copy_file_to_storage as copy7_001{
  #  input:
  #    text_file = combine_missense_001.out1,
  #    output_dir = step_9b_output_dir
  #}
  #call Tasks.copy_file_to_storage as copy7_0001{
  #  input:
  #    text_file = combine_missense_0001.out1,
  #    output_dir = step_9b_output_dir
  #}

  #call make_group_file as make_group_file_001 {
  #  input:
  #    tier1 = Tier1_Concat_001.out2,
  #    tier2 = Tier2_Concat_001.out2,
  #    tier3 = Tier3_Concat_001.out2,
  #    tier4 = Tier4_Concat_001.out2,
  #    tier5 = Tier5_Concat_001.out2,
  #    tier6 = Tier6_Concat_001.out2,
  #    af = "001"
  #}
  #call make_group_file as make_group_file_0001 {
  #  input:
  #    tier1 = Tier1_Concat_0001.out2,
  #    tier2 = Tier2_Concat_0001.out2,
  #    tier3 = Tier3_Concat_0001.out2,
  #    tier4 = Tier4_Concat_0001.out2,
  #    tier5 = Tier5_Concat_0001.out2,
  #    tier6 = Tier6_Concat_0001.out2,
  #    af = "0001"
  #}
}

task make_group_file {
  input {
    File tier1
    File tier2
    File tier3
    File tier4
    File tier5
    File tier6
    String af
  }
  command <<<
  set -euxo pipefail
  python3 <<CODE
  import pandas as pd
  from collections import defaultdict

  # Ordered from highest to lowest
  tiers = ["~{tier1}", "~{tier2}", "~{tier3}", "~{tier4}", "~{tier5}", "~{tier6}"]
  tier_labels = {f"tier{i}": f"T{i}" for i in range(1, 7)}

  # Dictionary: gene -> {variant: tier_label}
  gene_to_variants = defaultdict(dict)
  print("chek1")

  # Process each tier file (highest first)
  for rank, tier_file in enumerate(tiers, start=1):
      tier_name = f"tier{rank}"
      tier_label = tier_labels[tier_name]
      print(tier_file, tier_name)
      df = pd.read_csv(tier_file, sep="\t", names=["GENE", "VARIANT"], comment="#", dtype=str)
      for _, row in df.iterrows():
          gene = row["GENE"].strip()
          var = row["VARIANT"].strip().replace('_',':')

          # Only assign if not already labeled by a higher-priority tier
          if var not in gene_to_variants[gene]:
              gene_to_variants[gene][var] = tier_label
  print("check2")

  # --- Write output ---
  with open("ufc_~{af}.groupfile", "w") as out:
      for gene, variants in sorted(gene_to_variants.items()):
          vars_sorted = sorted(variants.keys())
          annos = [variants[v] for v in vars_sorted]
          out.write(f"{gene}\tvar\t{' '.join(vars_sorted)}\n")
          out.write(f"{gene}\tanno\t{' '.join(annos)}\n")
  CODE
  >>>
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible:3
  }
  output {
    File out1 = "ufc_~{af}.groupfile"
  }
}
task combine_missense{
  input{
    File tier3
    File tier4
  }
  command <<<
  set -euxo pipefail
  cat ~{tier3} ~{tier4} | sort -u > tier3and4.tsv
  >>>
  output {
    File out1 = "tier3and4.tsv"
  }
  runtime {
    docker: "ubuntu:latest"
    preemptible:3
  }
}
task T0_Initial_Filter{
  input{
    File vcf
    File subjects_list
  }
  command <<<
  set -euxo pipefail

  basename=$(basename "~{vcf}")
  chr=$(echo "$basename" | sed -E 's/.*(chr[0-9XYM]+).*/\1/')
  shard=$(echo "$basename" | sed -E 's/.*finalrun\.([0-9]+)\..*/\1/')

  bcftools view -S ~{subjects_list} ~{vcf} -Oz -o tmp.vcf.gz
  bcftools index -t tmp.vcf.gz

  # Extract relevant fields from VEP-annotated VCF
  echo -e 'ID\tAFR_AF\tEAS_AF\tAMR_AF\tFIN_AF\tNFE_AF\tSAS_AF\tMID_AF\tAMI_AF\tASJ_AF\tOTH_AF\tAF\tBIOTYPE\n' > extracted_variants.txt
  bcftools +split-vep tmp.vcf.gz \
  -f '%ID\t%gnomAD_AF_non_cancer_afr\t%gnomAD_AF_non_cancer_eas\t%gnomAD_AF_non_cancer_amr\t%gnomAD_AF_non_cancer_fin\t%gnomAD_AF_non_cancer_nfe\t%gnomAD_AF_non_cancer_sas\t%gnomAD_AF_non_cancer_mid\t%gnomAD_AF_non_cancer_ami\t%gnomAD_AF_non_cancer_asj\t%gnomAD_AF_non_cancer_oth\t%AF\t%BIOTYPE\n' \
  -d >> extracted_variants.txt || true

  python3 <<CODE
  import pandas as pd
  import numpy as np
  import sys

  # --- Load file ---
  #df = pd.read_csv("extracted_variants.txt", sep="\t")

  try:
      df = pd.read_csv("extracted_variants.txt", sep="\t")
  except pd.errors.EmptyDataError:
      print("No data found in extracted_variants.txt — skipping.")
      df = pd.DataFrame(columns=["ID","AFR_AF","EAS_AF","AMR_AF","FIN_AF","NFE_AF","SAS_AF",
                                 "MID_AF","AMI_AF","ASJ_AF","OTH_AF","AF","BIOTYPE"])

  if df.empty:
      print("No variants to process — exiting cleanly.")
      for fn in [
          "filtered_variants.protein_coding.001.txt",
          "filtered_variants.protein_coding.0001.txt",
          "filtered_variants.all_biotypes.001.txt",
          "filtered_variants.all_biotypes.0001.txt",
      ]:
          open(fn, "w").close()
      sys.exit(0)


  # --- Define population groupings ---
  major_cols = ["AFR_AF", "EAS_AF", "AMR_AF", "FIN_AF", "NFE_AF", "SAS_AF"]
  minor_cols = ["MID_AF", "AMI_AF", "ASJ_AF", "OTH_AF"]
  cohort_col = "AF"
  biotype_col = "BIOTYPE"

  # --- Convert missing values "." → 0 ---
  freq_cols = major_cols + minor_cols + [cohort_col]
  df[freq_cols] = df[freq_cols].replace(".", 0).apply(pd.to_numeric, errors="coerce").fillna(0)

  # --- Tier 1 thresholds ---
  mask_major_t1 = (df[major_cols] < 0.01).all(axis=1)
  mask_minor_t1 = (df[minor_cols] < 0.10).all(axis=1)
  mask_cohort_t1 = df[cohort_col] < 0.01
  mask_keep_t1 = mask_major_t1 & mask_minor_t1 & mask_cohort_t1

  # --- Tier 2 (stricter) thresholds ---
  mask_major_t2 = (df[major_cols] < 0.001).all(axis=1)
  mask_minor_t2 = (df[minor_cols] < 0.01).all(axis=1)
  mask_cohort_t2 = df[cohort_col] < 0.001
  mask_keep_t2 = mask_major_t2 & mask_minor_t2 & mask_cohort_t2

  # --- Define biotype mask ---
  mask_protein = df[biotype_col] == "protein_coding"

  # --- Tier 1 outputs ---
  # With protein_coding restriction
  df.loc[mask_keep_t1 & mask_protein, "ID"].drop_duplicates().to_csv(
      "filtered_variants.protein_coding.001.txt", sep="\t", index=False, header=False
  )

  # With protein_coding restriction
  df.loc[mask_keep_t2 & mask_protein, "ID"].drop_duplicates().to_csv(
      "filtered_variants.protein_coding.0001.txt", sep="\t", index=False, header=False
  )

  # Without protein_coding restriction
  df.loc[mask_keep_t1, "ID"].drop_duplicates().to_csv(
      "filtered_variants.all_biotypes.001.txt", sep="\t", index=False, header=False
  )

  # Without protein_coding restriction
  df.loc[mask_keep_t2, "ID"].drop_duplicates().to_csv(
      "filtered_variants.all_biotypes.0001.txt", sep="\t", index=False, header=False
  )

  CODE
  >>>
  output {
    File out1 = "filtered_variants.protein_coding.001.txt"
    File out2 = "filtered_variants.protein_coding.0001.txt"
    File splice_out1 = "filtered_variants.all_biotypes.001.txt"
    File splice_out2 = "filtered_variants.all_biotypes.0001.txt"
  }
  runtime {
    docker: 'vanallenlab/g2c_pipeline'
    #memory:"8GB"
    #disks: "local-disk 24 HDD"
    preemptible: 3
  }
}

task T0_Filter_Tier0 {
  input {
    File vcf
    File rare_variants
  }

  command <<<
  set -euxo pipefail
  # Get necessary information
  bcftools view --include ID==@~{rare_variants} ~{vcf} -G -O z -o tmp.vcf.gz
  echo -e 'GENE\tID\tCLINVAR\n' > tmp0.txt
  bcftools +split-vep tmp.vcf.gz -f '%SYMBOL\t%ID\t%clinvar_CLNSIG\n' -d >> tmp0.txt

  python3 <<CODE
  import pandas as pd

  df = pd.read_csv("tmp0.txt",sep='\t',index_col=False)
  df = df[df['CLINVAR'] == 'Uncertain_significance']
  df[["GENE",'ID']].drop_duplicates().to_csv("filtered_variants.tier0.txt",header=False,index=False,sep='\t')
  CODE
  >>>
  output {
    File out1 = "filtered_variants.tier0.txt"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task T1_Filter_Tier1 {
  input {
    File vcf
    File rare_variants
  }

  command <<<
  set -euxo pipefail
  # Get necessary information
  bcftools view --include ID==@~{rare_variants} ~{vcf} -G -O z -o tmp.vcf.gz
  echo -e 'GENE\tID\tIMPACT\tCLINVAR\tCONSEQUENCE\tBIOTYPE\tLOFTEE\n' > tmp0.txt
  bcftools +split-vep tmp.vcf.gz -f '%SYMBOL\t%ID\t%IMPACT\t%clinvar_CLNSIG\t%Consequence\t%BIOTYPE\t%LoF\n' -d >> tmp0.txt
  
  python3 <<CODE
  import pandas as pd
  
  df = pd.read_csv("tmp0.txt",sep='\t',index_col=False)
  df = df[df['BIOTYPE'] == "protein_coding"]
  df = df[~df['CLINVAR'].str.contains('benign',case=False,na=False)]
  df = df[((df['IMPACT'] == 'HIGH') & (df['LOFTEE'] == "HC"))]
  #df = df[((df['IMPACT'] == 'HIGH') & (df['LOFTEE'] == "HC")) | (df['CLINVAR'] == "Pathogenic") | (df['CLINVAR'] == "Likely_pathogenic") | (df['CLINVAR'] == "Pathogenic/Likely_pathogenic")]
  df[['GENE','ID']].drop_duplicates().to_csv("filtered_variants.txt",header=False,index=False,sep='\t')
  CODE
  >>>
  output {
    File out1 = "filtered_variants.txt"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task T2_Filter_Tier2 {
  input {
    File vcf
    File rare_variants
  }

  command <<<
  set -euxo pipefail

  # Get necessary information
  # Extract necessary VEP + SpliceAI info
  bcftools view --include ID==@~{rare_variants} ~{vcf} -G -O z -o tmp.vcf.gz
  echo -e 'GENE\tID\tIMPACT\tCLINVAR\tCONSEQUENCE\tBIOTYPE\tSPLICE_AG\tSPLICE_AL\tSPLICE_DG\tSPLICE_DL\tSPLICE_GENE\tLOFTEE\n' > tmp0.txt
  bcftools +split-vep tmp.vcf.gz -f '%SYMBOL\t%ID\t%SpliceAI_pred_DS_AG\t%SpliceAI_pred_DS_AL\t%SpliceAI_pred_DS_DG\t%SpliceAI_pred_DS_DL\t%SpliceAI_pred_SYMBOL\t%LoF\n' -d >> tmp0.txt

  python3 <<CODE
  import pandas as pd
  df = pd.read_csv("tmp0.txt",sep='\t',index_col=False)

  # Convert SpliceAI score columns to numeric (coerce invalids to NaN)
  score_cols = ["SPLICE_AG", "SPLICE_AL", "SPLICE_DG", "SPLICE_DL"]
  df[score_cols] = df[score_cols].apply(pd.to_numeric, errors="coerce")

  # Filter: max score > 0.5 and matching gene symbol
  df_filtered = df[(df[score_cols].max(axis=1) >= 0.5) | (df['LOFTEE'] == "LC")]

  # Write filtered results back
  df_filtered[['SPLICE_GENE','ID']].drop_duplicates().to_csv("splice_variants.txt", sep='\t', index=False,header=False)
  CODE

  >>>
  output {
    File out1 = "splice_variants.txt"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}
task T34_Filter_Tier34 {
  input {
    File vcf
    File rare_variants
  }

  command <<<
  set -euxo pipefail

  # Get necessary information  
  bcftools view --include ID==@~{rare_variants} ~{vcf} -G -O z -o tmp.vcf.gz
  echo -e 'GENE\tID\tIMPACT\tCLINVAR\tCONSEQUENCE\tBIOTYPE\tAM_CLASS\tPRIMATE_AI_PRED\tREVEL_SCORE\n' > tmp0.txt
  bcftools +split-vep tmp.vcf.gz \
    -f '%SYMBOL\t%ID\t%IMPACT\t%clinvar_CLNSIG\t%Consequence\t%BIOTYPE\t%am_class\t%PrimateAI_pred\t%REVEL_score\n' -d >> tmp0.txt
  # CODE to parse VEP REVEL Scores
  python3 <<CODE
  import pandas as pd
  import numpy as np

  def process_revel_score(score_field: str) -> float:
      parts = score_field.split("&")
      numeric_values = [float(p) for p in parts if p != "."]
      return max(numeric_values) if numeric_values else None

  df = pd.read_csv("tmp0.txt", sep="\t")
  df["max_revel_score"] = df["REVEL_SCORE"].apply(process_revel_score)
  df = df[~df['CLINVAR'].str.contains("benign",case=False,na=False)]
  df = df[df['BIOTYPE'] == "protein_coding"]
  df = df[df['IMPACT'] == "MODERATE"]
  df = df[df['CONSEQUENCE'].str.contains("missense_variant")]

  # --- Define criteria ---
  criteria = [
      df["max_revel_score"] >= 0.75,
      df["PRIMATE_AI_PRED"].eq("D"),
      df["AM_CLASS"].eq("likely_pathogenic")
  ]
  clinvar_criteria = df["CLINVAR"].isin(["Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic"])

  # Count how many of the three are true
  df["criteria_met"] = sum(criteria)
  # --- Filter: keep variants meeting >=2 of the 3 criteria ---
  filtered = df[(df["criteria_met"] >= 2) | clinvar_criteria].copy()
  filtered[['GENE','ID']].drop_duplicates().to_csv("filtered_variants.tier3.txt",index=False,header=False,sep='\t')

  # --- Define criteria 4 ---
  criteria = [
      df["max_revel_score"] >= 0.5,
      df["PRIMATE_AI_PRED"].eq("D"),
      df["AM_CLASS"].eq("likely_pathogenic")
  ]

  # Count how many of the three are true
  df["criteria_met"] = sum(criteria)
  # --- Filter: keep variants meeting >=2 of the 3 criteria ---
  filtered = df[(df["criteria_met"] >= 1) | clinvar_criteria].copy()
  filtered[['GENE','ID']].drop_duplicates().to_csv("filtered_variants.tier4.txt",index=False,header=False,sep='\t')

  # --- Define criteria 5 ---
  criteria = [
      df["max_revel_score"] >= 0.5,
      df["PRIMATE_AI_PRED"].eq("D"),
      df["AM_CLASS"].eq("likely_pathogenic")
  ]

  # Count how many of the three are true
  df["criteria_met"] = sum(criteria)
  # --- Filter: keep variants meeting >=2 of the 3 criteria ---
  filtered = df[(df["criteria_met"] == 0) & ~clinvar_criteria].copy()
  filtered[['GENE','ID']].drop_duplicates().to_csv("filtered_variants.tier5.txt",index=False,header=False,sep='\t')
  CODE
  >>>
  output {
    File out1 = "filtered_variants.tier3.txt"
    File out2 = "filtered_variants.tier4.txt"
    File out3 = "filtered_variants.tier5.txt"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task T6_Filter_Tier6 {
  input {
    File vcf
    File rare_variants
  }

  command <<<
  set -euxo pipefail

  # Get necessary information
  bcftools view --include ID==@~{rare_variants} ~{vcf} -G -O z -o tmp.vcf.gz
  echo -e 'GENE\tID\tIMPACT\tCLINVAR\tCONSEQUENCE\tBIOTYPE\n' > tmp0.txt
  bcftools +split-vep tmp.vcf.gz \
    -f '%SYMBOL\t%ID\t%IMPACT\t%clinvar_CLNSIG\t%Consequence\t%BIOTYPE\n' -d >> tmp0.txt
  # CODE to parse VEP REVEL Scores
  python3 <<CODE
  import pandas as pd
  import numpy as np


  df = pd.read_csv("tmp0.txt", sep="\t")
  #df = df[~df['CLNVAR'].str.contains("benign",case=False,na=False)]
  df = df[df['BIOTYPE'] == "protein_coding"]
  df = df[df['IMPACT'] == "LOW"]
  df = df[df['CONSEQUENCE'].str.contains("synonymous_variant")]
  df[['GENE','ID']].drop_duplicates().to_csv("filtered_variants.tier6.txt",index=False,header=False,sep='\t')
  CODE
  >>>
  output {
    File out1 = "filtered_variants.tier6.txt"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}
