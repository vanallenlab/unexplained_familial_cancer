# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# # Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow STEP_R1_TIER_VARIANTS {
  input {
    String CDRv9_cpg_vcfs = "gs://dataproc-staging-wb-cordial-diamond-9893/vep_annotated_cpgs_sept2026/CDRv9"
    String CDRv9_cpg_genotyped_vcfs = "gs://dataproc-staging-wb-cordial-diamond-9893/unannotated_genotyped_cpgs_sept2026/CDRv9"
    String step_r1_output_dir = "gs://dataproc-staging-wb-cordial-diamond-9893/step_r1_output/"
    Array[String] genes = ["BRCA1","BRCA2","CHEK2","PALB2","ATM","CDH1","STK11","TP53","PTEN","BARD1","RAD51C","RAD51D","NF1","VHL","TSC1","TSC2","FH","FLCN","BAP1","SDHA","SDHB","SDHC","SDHD"]
  }
  
  scatter (gene in genes){

    # File to look up
    File genotyped_vcf = CDRv9_cpg_genotyped_vcfs + "/" + gene + ".gt.vcf.bgz"

    call gather_vcf {
      input:
        dir = CDRv9_cpg_vcfs,
        gene = gene
    }

    call T0_Initial_Filter{
      input:
        vcf = gather_vcf.out1
    }

    call T1_Filter_Tier1 as T1_Filter_Tier1_0001{
      input:
        rare_variants = T0_Initial_Filter.out1,
        vcf = gather_vcf.out1,
        gene = gene
    }

    call TX_Lookup as TX_Lookup_T1 {
      input:
        variants = T1_Filter_Tier1_0001.out1,
        vcf = genotyped_vcf,
        gene = gene,
        tier = "1"
    }

    call Tasks.copy_file_to_storage as copy1{
      input:
        text_file = TX_Lookup_T1.out1,
        output_dir = step_r1_output_dir
    }

    call T2_Filter_Tier2 as T2_Filter_Tier2_0001{
      input:
        rare_variants = T0_Initial_Filter.out1,
        vcf = gather_vcf.out1,
        gene = gene
    }

    call TX_Lookup as TX_Lookup_T2 {
      input:
        variants = T2_Filter_Tier2_0001.out1,
        vcf = genotyped_vcf,
        gene = gene,
        tier = "2"
    }

    call Tasks.copy_file_to_storage as copy2{
      input:
        text_file = TX_Lookup_T2.out1,
        output_dir = step_r1_output_dir
    }

    call T34_Filter_Tier34 as T34_Filter_Tier34_0001{
      input:
        vcf = gather_vcf.out1,
        rare_variants = T0_Initial_Filter.out1,
        gene = gene
    }

    call TX_Lookup as TX_Lookup_T3 {
      input:
        variants = T34_Filter_Tier34_0001.out1,
        vcf = genotyped_vcf,
        gene = gene,
        tier = "3"
    }

    call Tasks.copy_file_to_storage as copy3{
      input:
        text_file = TX_Lookup_T3.out1,
        output_dir = step_r1_output_dir
    }

    call TX_Lookup as TX_Lookup_T4 {
      input:
        variants = T34_Filter_Tier34_0001.out2,
        vcf = genotyped_vcf,
        gene = gene,
        tier = "4"
    }

    call Tasks.copy_file_to_storage as copy4{
      input:
        text_file = TX_Lookup_T4.out1,
        output_dir = step_r1_output_dir
    }

    call TX_Lookup as TX_Lookup_T5 {
      input:
        variants = T34_Filter_Tier34_0001.out3,
        vcf = genotyped_vcf,
        gene = gene,
        tier = "5"
    }

    call Tasks.copy_file_to_storage as copy5{
      input:
        text_file = TX_Lookup_T5.out1,
        output_dir = step_r1_output_dir
    }

    call T6_Filter_Tier6 as T6_Filter_Tier6_0001{
      input:
        vcf = gather_vcf.out1,
        rare_variants = T0_Initial_Filter.out1,
        gene = gene
    }

    call TX_Lookup as TX_Lookup_T6 {
      input:
        variants = T6_Filter_Tier6_0001.out1,
        vcf = genotyped_vcf,
        gene = gene,
        tier = "6"
    }

    call Tasks.copy_file_to_storage as copy6{
      input:
        text_file = TX_Lookup_T6.out1,
        output_dir = step_r1_output_dir
    }
  }

  #call Tasks.copy_file_to_storage as copy1{
  #  input:
  #    text_file = T1_Filter_Tier1_0001.out1,
  #    output_dir = step_r1_output_dir
  #}
  #call Tasks.copy_file_to_storage as copy3{
  #  input:
  #    text_file = T34_Filter_Tier34_0001.out1,
  #    output_dir = step_r1_output_dir
  #}
  #call Tasks.copy_file_to_storage as copy4{
  #  input:
  #    text_file = T34_Filter_Tier34_0001.out1,
  #    output_dir = step_r1_output_dir
  #}
  #call Tasks.copy_file_to_storage as copy5{
  #  input:
  #    text_file = T34_Filter_Tier34_0001.out3,
  #    output_dir = step_r1_output_dir
  #}
  #call Tasks.copy_file_to_storage as copy6{
  #  input:
  #    text_file = T6_Filter_Tier6_0001.out1,
  #    output_dir = step_9b_output_dir
  #}
}

task gather_vcf {
  input {
    String dir
    String gene
  }
  command <<<
  gsutil ls ~{dir}/*.vcf.bgz | grep -E "/~{gene}([_.]|$)" > vcf.list
  >>>
  output {
    File out1 = read_string("vcf.list")
  }
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:588.0.0"
    preemptible: 3
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
  }
  command <<<
  set -euxo pipefail

  basename=$(basename "~{vcf}")
  chr=$(echo "$basename" | sed -E 's/.*(chr[0-9XYM]+).*/\1/')
  shard=$(echo "$basename" | sed -E 's/.*finalrun\.([0-9]+)\..*/\1/')

  bcftools view ~{vcf} -Oz -o tmp.vcf.gz
  bcftools index -t tmp.vcf.gz

  # Extract relevant fields from VEP-annotated VCF
  echo -e 'ID\tCHROM\tPOS\tREF\tALT\tAFR_AF\tEAS_AF\tAMR_AF\tFIN_AF\tNFE_AF\tSAS_AF\tMID_AF\tAMI_AF\tASJ_AF\tOTH_AF\tAF\tBIOTYPE\tclinvar_clnsig\n' > extracted_variants.txt
  bcftools +split-vep tmp.vcf.gz \
  -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%gnomAD_AF_non_cancer_afr\t%gnomAD_AF_non_cancer_eas\t%gnomAD_AF_non_cancer_amr\t%gnomAD_AF_non_cancer_fin\t%gnomAD_AF_non_cancer_nfe\t%gnomAD_AF_non_cancer_sas\t%gnomAD_AF_non_cancer_mid\t%gnomAD_AF_non_cancer_ami\t%gnomAD_AF_non_cancer_asj\t%gnomAD_AF_non_cancer_oth\t%AF\t%BIOTYPE\t%clinvar_CLNSIG\n' \
  -d >> extracted_variants.txt || true

  python3 <<CODE
  import pandas as pd
  import numpy as np
  import sys

  # --- Load file ---
  df = pd.read_csv("extracted_variants.txt", sep="\t")
  df = df[df["BIOTYPE"] == "protein_coding"]

  # --- Define population groupings ---
  major_cols = ["AFR_AF", "EAS_AF", "AMR_AF", "FIN_AF", "NFE_AF", "SAS_AF"]
  minor_cols = ["MID_AF", "AMI_AF", "ASJ_AF", "OTH_AF"]
  cohort_col = "AF"
  gnomad_cols = major_cols + minor_cols

  # Identify SNVs vs indels
  is_snv = (df["REF"].str.len() == 1) & (df["ALT"].str.len() == 1)
  is_indel = ~is_snv

  # --- Require gnomAD frequencies for indels before converting missing values ---
  indel_has_gnomad = is_indel & (df[gnomad_cols] != ".").all(axis=1)
  indel_has_clinvar = is_indel & (df["clinvar_clnsig"] != ".")
  mask_ok = is_snv | indel_has_gnomad | indel_has_clinvar

  # --- Convert missing values "." → 0 ---
  freq_cols = major_cols + minor_cols + [cohort_col]
  df[freq_cols] = df[freq_cols].replace(".", 0).apply(pd.to_numeric, errors="coerce").fillna(0)

  # --- Ultra Rare (AF < 0.1%) Variants ---
  mask_major_t1 = (df[major_cols] < 0.001).all(axis=1)
  mask_minor_t1 = (df[minor_cols] < 0.01).all(axis=1)
  mask_cohort_t1 = df[cohort_col] < 0.001
  mask_keep_t1 = mask_major_t1 & mask_minor_t1 & mask_cohort_t1

  # With protein_coding restriction
  df.loc[mask_keep_t1, "ID"].drop_duplicates().to_csv(
      "filtered_variants.protein_coding.0001.txt", sep="\t", index=False, header=False
  )

  CODE
  >>>
  output {
    File out1 = "filtered_variants.protein_coding.0001.txt"
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
  df = df[(df['CLINVAR'] == 'Uncertain_significance')]
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
    String gene
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
  df[['ID']].drop_duplicates().to_csv("~{gene}.tier1_0001.txt",header=False,index=False,sep='\t')
  CODE
  >>>
  output {
    File out1 = "~{gene}.tier1_0001.txt"
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
    String gene
  }

  command <<<
  set -euxo pipefail

  # Get necessary information
  # Extract necessary VEP + SpliceAI info
  bcftools view --include ID==@~{rare_variants} ~{vcf} -G -O z -o tmp.vcf.gz
  echo -e 'GENE\tID\tSPLICE_AG\tSPLICE_AL\tSPLICE_DG\tSPLICE_DL\tSPLICE_GENE\tLOFTEE\tCLINVAR\tCONSEQUENCE\tBIOTYPE\n' > tmp0.txt
  bcftools +split-vep tmp.vcf.gz -f '%SYMBOL\t%ID\t%SpliceAI_pred_DS_AG\t%SpliceAI_pred_DS_AL\t%SpliceAI_pred_DS_DG\t%SpliceAI_pred_DS_DL\t%SpliceAI_pred_SYMBOL\t%LoF\t%clinvar_CLNSIG\t%Consequence\t%BIOTYPE\n' -d >> tmp0.txt

  python3 <<CODE
  import pandas as pd
  df = pd.read_csv("tmp0.txt",sep='\t',index_col=False)
  df = df[df['SPLICE_GENE'] == "~{gene}"]

  # Identify IDs to remove:
  # those that have ANY row where (missense_variant AND protein_coding)
  ids_to_remove = (
      df.groupby("ID")
        .filter(lambda g: ((g["CONSEQUENCE"].str.contains("missense_variant")) &
                           (g["BIOTYPE"] == "protein_coding")).any())
        ["ID"].unique()
  )

  # Drop those IDs from df
  df = df[~df["ID"].isin(ids_to_remove)]

  # Convert SpliceAI score columns to numeric (coerce invalids to NaN)
  score_cols = ["SPLICE_AG", "SPLICE_AL", "SPLICE_DG", "SPLICE_DL"]
  df[score_cols] = df[score_cols].apply(pd.to_numeric, errors="coerce")

  # Filter: max score > 0.5 and matching gene symbol
  df = df[~df['CLINVAR'].str.contains('benign',case=False,na=False)]
  df_filtered = df[~(df['LOFTEE'] == "HC") & ((df[score_cols].max(axis=1) >= 0.5) | (df['LOFTEE'] == "LC") | ((df[score_cols].max(axis=1) >= 0.2) & (df['CLINVAR'].isin(['Pathogenic','Likely_pathogenic','Pathogenic/Likely_pathogenic']))))]
  df_filtered = df_filtered[df_filtered['SPLICE_GENE'] != "."]

  # Write filtered results back
  df_filtered[['SPLICE_GENE','ID']].drop_duplicates().to_csv("~{gene}.tier2_0001.txt", sep='\t', index=False,header=False)

  CODE

  >>>
  output {
    File out1 = "~{gene}.tier2_0001.txt"
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
    String gene
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
  filtered[['ID']].drop_duplicates().to_csv("~{gene}.tier3_0001.txt",index=False,header=False,sep='\t')

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
  filtered[['ID']].drop_duplicates().to_csv("~{gene}.tier4_0001.txt",index=False,header=False,sep='\t')

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
  filtered[['ID']].drop_duplicates().to_csv("~{gene}.tier5_0001.txt",index=False,header=False,sep='\t')
  CODE
  >>>
  output {
    File out1 = "~{gene}.tier3_0001.txt"
    File out2 = "~{gene}.tier4_0001.txt"
    File out3 = "~{gene}.tier5_0001.txt"
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
    String gene
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
  df[['GENE','ID']].drop_duplicates().to_csv("~{gene}.tier6_0001.txt",index=False,header=False,sep='\t')
  CODE
  >>>
  output {
    File out1 = "~{gene}.tier6_0001.txt"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task TX_Lookup {
  input {
    File vcf
    File variants
    String gene
    String tier
  } 
  String vcf_basename=basename("~{vcf}")
  command <<<
  set -euxo pipefail
  if [ $(wc -l < ~{variants}) -eq 0 ]; then
    echo "No Variants To Keep"
    touch ~{gene}.tier~{tier}_0001.list
    exit 0
  fi
  bcftools view --include ID==@~{variants} ~{vcf} -O z -o tmp.vcf.gz
  bcftools query -f '[%SAMPLE\t%GT\n]' tmp.vcf.gz | awk '$2 !~ /^(\.\/\.|\.\|\.)$/ && $2 !~ /^0[\/|]0$/ {print $1}' | sort -u > ~{gene}.tier~{tier}_0001.list
  >>>
  output {
    File out1 = "~{gene}.tier~{tier}_0001.list"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
  }
} 
