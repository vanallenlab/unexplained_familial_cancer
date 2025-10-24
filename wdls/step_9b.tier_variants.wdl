# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# # Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow STEP_9_TIER_VARIANTS {
  input {
    String step_9_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/sharded_vcfs"
    File subjects_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cohorts/ufc_subjects.aou.list"
    String step_9b_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/"
  }

  call Tasks.gather_vcfs {
    input:
      dir = step_9_output_dir
  }

  call Tasks.sort_vcf_list {
    input:
      unsorted_vcf_list = gather_vcfs.vcf_list
  }

  Int negative_shards = 2958
  scatter (i in range(length(sort_vcf_list.vcf_arr)-negative_shards)){
    call T0_Initial_Filter{
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        subjects_list = subjects_list
    }
    call T0_Filter_Tier0{
      input:
        rare_variants = T0_Initial_Filter.out1,
        vcf = sort_vcf_list.vcf_arr[i]
    }
    call T1_Filter_Tier1{
      input:
        rare_variants = T0_Initial_Filter.out1,
        vcf = sort_vcf_list.vcf_arr[i]
    }
    call T34_Filter_Tier34 {
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        rare_variants = T0_Initial_Filter.out1,
    }
    call T6_Filter_Tier6 {
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        rare_variants = T0_Initial_Filter.out1,
    }
  }

  call Tasks.concatenateFiles as Tier0_Concat {
    input:
      files = T1_Filter_Tier1.out1,
      output_name = "tier1"
  }

  call Tasks.copy_file_to_storage as copy0{
    input:
      text_file = Tier0_Concat.out2,
      output_dir = step_9b_output_dir
  }

  call Tasks.concatenateFiles as Tier1_Concat {
    input:
      files = T1_Filter_Tier1.out1,
      output_name = "tier1"
  }

  call Tasks.copy_file_to_storage as copy1{
    input:
      text_file = Tier1_Concat.out2,
      output_dir = step_9b_output_dir
  }

  call Tasks.concatenateFiles as Tier3_Concat {
    input:
      files = T34_Filter_Tier34.out1,
      output_name = "tier3"
  }

  call Tasks.copy_file_to_storage as copy3{
    input:
      text_file = Tier3_Concat.out2,
      output_dir = step_9b_output_dir
  }

  call Tasks.concatenateFiles as Tier4_Concat {
    input:
      files = T34_Filter_Tier34.out2,
      output_name = "tier4"
  }

  call Tasks.copy_file_to_storage as copy4{
    input:
      text_file = Tier4_Concat.out2,
      output_dir = step_9b_output_dir
  }
  call Tasks.concatenateFiles as Tier5_Concat {
    input:
      files = T34_Filter_Tier34.out3,
      output_name = "tier5"
  }

  call Tasks.copy_file_to_storage as copy5{
    input:
      text_file = Tier5_Concat.out2,
      output_dir = step_9b_output_dir
  }

  call Tasks.concatenateFiles as Tier6_Concat {
    input:
      files = T6_Filter_Tier6.out1,
      output_name = "tier6"
  }

  call Tasks.copy_file_to_storage as copy6{
    input:
      text_file = Tier6_Concat.out2,
      output_dir = step_9b_output_dir
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
  bcftools +split-vep tmp.vcf.gz \
  -f '%ID\t%gnomAD_AF_non_cancer_afr\t%gnomAD_AF_non_cancer_eas\t%gnomAD_AF_non_cancer_amr\t%gnomAD_AF_non_cancer_fin\t%gnomAD_AF_non_cancer_nfe\t%gnomAD_AF_non_cancer_sas\t%gnomAD_AF_non_cancer_mid\t%gnomAD_AF_non_cancer_ami\t%gnomAD_AF_non_cancer_asj\t%gnomAD_AF_non_cancer_oth\t%AF\t%BIOTYPE\n' \
  -d 2>/dev/null | sort -u > extracted_variants.txt

  # Define thresholds
  awk -F'\t' '
  BEGIN {
      OFS="\t";
      # major populations (<1%)
      major_idx[1]=2; major_idx[2]=3; major_idx[3]=4; major_idx[4]=5; major_idx[5]=6; major_idx[6]=7;
      # minor populations (<10%)
      minor_idx[1]=8; minor_idx[2]=9; minor_idx[3]=10; minor_idx[4]=11;
      cohort_idx = 12; # cohort AF (<1%)
      biotype_idx = 13;
  }
  {
      # convert missing AFs (".") to 0
      for (i=2; i<=12; i++) {
          if ($i==".") $i=0;
          $i += 0;
      }

      # only consider protein-coding
      if ($biotype_idx != "protein_coding") next;

      # check major populations
      keep=1
      for (i in major_idx) {
          if ($major_idx[i] >= 0.01) { keep=0; break }
      }
      # check minor populations
      if (keep) {
          for (i in minor_idx) {
              if ($minor_idx[i] >= 0.10) { keep=0; break }
          }
      }
      # check cohort
      if (keep && $cohort_idx >= 0.01) keep=0;

      if (keep) print $1
  }' extracted_variants.txt > filtered_variants.txt

  # Filter to relevant variants
  #if [ -s filtered_variants.txt ]; then
  #  bcftools view --include ID==@filtered_variants.txt ~{vcf} -O z -o $basename
  #else
  #  bcftools view -h ~{vcf} -Oz -o $basename
  #fi
  >>>
  output {
    File out1 = "filtered_variants.txt"
  }
  runtime {
    docker: 'vanallenlab/bcftools'
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
  df[['ID']].drop_duplicates().to_csv("filtered_variants.tier0.txt",header=False,index=False)
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
  echo -e 'GENE\tID\tIMPACT\tCLINVAR\tCONSEQUENCE\tBIOTYPE\n' > tmp0.txt
  bcftools +split-vep tmp.vcf.gz -f '%SYMBOL\t%ID\t%IMPACT\t%clinvar_CLNSIG\t%Consequence\t%BIOTYPE\n' -d >> tmp0.txt
  
  python3 <<CODE
  import pandas as pd
  
  df = pd.read_csv("tmp0.txt",sep='\t',index_col=False)
  df = df[df['BIOTYPE'] == "protein_coding"]
  df = df[~df['CLINVAR'].str.contains('benign',case=False,na=False)]
  df = df[(df['IMPACT'] == 'HIGH') | (df['CLINVAR'].str.contains("Pathogenic",na=False)) | (df['CLINVAR'].str.contains("Likely_pathogenic",na=False))]
  df[['ID']].drop_duplicates().to_csv("filtered_variants.txt",header=False,index=False)
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
  echo -e 'GENE\tID\tIMPACT\tCLINVAR\tCONSEQUENCE\tBIOTYPE\tSPLICE_AG\tSPLICE_DS\tSPLICE_DL\nSPLICE_GENE\n'
  bcftools +split-vep tmp.vcf.gz -f '%SYMBOL\t%ID\t%IMPACT\t%clinvar_CLNSIG\t%Consequence\t%BIOTYPE\t%SpliceAI_pred_DS_AG\t%SpliceAI_pred_DS_AL\t%SpliceAI_pred_DS_DG\t%SpliceAI_pred_DS_DL\t%SpliceAI_pred_SYMBOL\n' -d > tmp0.txt

  # Filter variants and compute max SpliceAI DS
  awk -F '\t' 'BEGIN { OFS="\t" }
  {
      symbol = $1
      impact = $3
      clin = $4
      consequence = $5
      biotype = $6
      splice_ai_gene = $11

      # SpliceAI DS fields
      ds[1] = ($7=="" ? 0 : $7)
      ds[2] = ($8=="" ? 0 : $8)
      ds[3] = ($9=="" ? 0 : $9)
      ds[4] = ($10=="" ? 0 : $10)

      # max DS
      max_ds = ds[1]
      for(i=2;i<=4;i++) if(ds[i] > max_ds) max_ds = ds[i]

      # HIGH-impact variant condition
      high_impact = (biotype=="protein_coding" && impact=="HIGH" && consequence !~ /NMD_transcript_variant/i && clin !~ /benign/i)

      # SpliceAI condition
      splice_ai = (max_ds > 0.5 && symbol==splice_ai_gene && clin !~ /benign/i)

      # Pathogenic variant condition
      pathogenic = (clin ~ /pathogenic/i && consequence !~ /LOW|MODIFIER/)
    
      if (high_impact || splice_ai || pathogenic) {
        print $2
      }
  }' tmp0.txt | sort -u > filtered_variants.txt


  # Filter to relevant variants
  #if [ -s filtered_variants.txt ]; then
  #  bcftools view --include ID==@filtered_variants.txt ~{vcf} -O z -o $output_vcf
  #  bcftools index -t "$output_vcf"
  #else
  #  bcftools view -h ~{vcf} -Oz -o $output_vcf
  #  bcftools index -t "$output_vcf"
  #fi
  >>>
  output {
    File out1 = "filtered_variants.txt"
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
  # Count how many of the three are true
  df["criteria_met"] = sum(criteria)
  # --- Filter: keep variants meeting >=2 of the 3 criteria ---
  filtered = df[df["criteria_met"] >= 2].copy()
  filtered[['ID']].drop_duplicates().to_csv("filtered_variants.tier3.txt",index=False,header=False)

  # --- Define criteria 4 ---
  criteria = [
      df["max_revel_score"] >= 0.5,
      df["PRIMATE_AI_PRED"].eq("D"),
      df["AM_CLASS"].eq("likely_pathogenic")
  ]

  # Count how many of the three are true
  df["criteria_met"] = sum(criteria)
  # --- Filter: keep variants meeting >=2 of the 3 criteria ---
  filtered = df[df["criteria_met"] >= 1].copy()
  filtered[['ID']].drop_duplicates().to_csv("filtered_variants.tier4.txt",index=False,header=False)

  # --- Define criteria 5 ---
  criteria = [
      df["max_revel_score"] >= 0.5,
      df["PRIMATE_AI_PRED"].eq("D"),
      df["AM_CLASS"].eq("likely_pathogenic")
  ]

  # Count how many of the three are true
  df["criteria_met"] = sum(criteria)
  # --- Filter: keep variants meeting >=2 of the 3 criteria ---
  filtered = df[df["criteria_met"] == 0].copy()
  filtered[['ID']].drop_duplicates().to_csv("filtered_variants.tier5.txt",index=False,header=False)
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
  df[['ID']].drop_duplicates().to_csv("filtered_variants.tier6.txt",index=False,header=False)
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
