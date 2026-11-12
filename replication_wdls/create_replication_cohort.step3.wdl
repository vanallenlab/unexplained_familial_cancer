# The Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow CREATE_REPLICATION_COHORT_STEP_3{
    input {
        String AoU_version = "CDRv9"
        String Google_Project = "wb-cordial-diamond-9893"
        Int positive_shards = 154
    }

    String output_name = "oct11_2026_pgv_carriers.~{AoU_version}.list"
    String vepped_cpgs_repository = "gs://dataproc-staging-wb-cordial-diamond-9893/vep_annotated_cpgs_sept2026/" + AoU_version
    String genotyped_cpgs_repository = "gs://dataproc-staging-wb-cordial-diamond-9893/unannotated_genotyped_cpgs_sept2026/" + AoU_version

    call gather_vcfs as T1_gather_annotated_vcfs{
        input:
          dir = vepped_cpgs_repository,
          positive_shards = positive_shards
    }

    call Tasks.sort_vcf_list as T2_sort_vcf_list{
        input:
            unsorted_vcf_list = T1_gather_annotated_vcfs.vcf_list
    }

    scatter(i in range(length(T2_sort_vcf_list.vcf_arr))){
        call T3_find_rare_pathogenic_variants {
            input:
                vcf = T2_sort_vcf_list.vcf_arr[i]
        }
    }	
    
    call Tasks.concatenateFiles as concatenateFiles1{
      input:
        files = T3_find_rare_pathogenic_variants.out1,
        output_name = "PGVs_in_CPGs." + AoU_version + ".list"
    }

    call gather_vcfs as T4_gather_genotyped_vcfs {
       input:
         dir = genotyped_cpgs_repository,
         positive_shards = positive_shards
    }

    call Tasks.sort_vcf_list as T5_sort_vcf_list{
        input:
            unsorted_vcf_list = T4_gather_genotyped_vcfs.vcf_list
    }

    scatter(i in range(length(T5_sort_vcf_list.vcf_arr))){
        call Tasks.Filter_VCF_USING_ID {
            input:
                vcf = T5_sort_vcf_list.vcf_arr[i],
                variant_list = concatenateFiles1.out2
        }
        call T7_get_pgv_carriers {
            input:
                vcf = Filter_VCF_USING_ID.out1
        }
    }

    call Tasks.concatenateFiles as concatenateFiles2{
      input:
        files = T7_get_pgv_carriers.out1,
        output_name = output_name
    }

    call Tasks.copy_file_to_storage as copy1{
      input:
        text_file = concatenateFiles1.out2,
        output_dir = "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/"
    }

    call Tasks.copy_file_to_storage as copy2{
      input:
        text_file = concatenateFiles2.out2,
        output_dir = "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/" 
    }   
}

task gather_vcfs {
  input {
    String dir
    Int positive_shards
  }
  command <<<
  gsutil ls ~{dir}/*.vcf.bgz | grep -Ev 'TSTD2|MYO1A|CCDC27|OLFML2B' | head -n ~{positive_shards} > vcf.list
  >>>
  output {
    File vcf_list = "vcf.list"
  }
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task T3_find_rare_pathogenic_variants{
  input{
    File vcf
  }

  String filename = basename(vcf)
  String gene = sub(filename, "\\.vep\\.vcf\\.bgz$", "")

  command <<<
  set -euxo pipefail

  # Annotate VCF with ID
  bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Oz -o annotated.vcf.bgz ~{vcf}
  tabix -p vcf annotated.vcf.bgz

  # Extract relevant fields from VEP-annotated VCF
  echo -e 'ID\tCHROM\tPOS\tREF\tALT\tAC\tAFR_AF\tEAS_AF\tAMR_AF\tFIN_AF\tNFE_AF\tSAS_AF\tMID_AF\tAMI_AF\tASJ_AF\tOTH_AF\tAF\tclinvar_clnsig\tIMPACT\n' > extracted_variants.txt
  bcftools +split-vep annotated.vcf.bgz \
  -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%AC\t%gnomAD_AF_non_cancer_afr\t%gnomAD_AF_non_cancer_eas\t%gnomAD_AF_non_cancer_amr\t%gnomAD_AF_non_cancer_fin\t%gnomAD_AF_non_cancer_nfe\t%gnomAD_AF_non_cancer_sas\t%gnomAD_AF_non_cancer_mid\t%gnomAD_AF_non_cancer_ami\t%gnomAD_AF_non_cancer_asj\t%gnomAD_AF_non_cancer_oth\t%AF\t%clinvar_CLNSIG\t%IMPACT\n' \
  -d >> extracted_variants.txt || true

  python3 <<CODE
  import pandas as pd
  import numpy as np
  import sys

  # --- Load file ---
  df = pd.read_csv("extracted_variants.txt", sep="\t")
  df = df[(df["clinvar_clnsig"].isin(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"])) | ((df["IMPACT"] == "HIGH") & (~df["clinvar_clnsig"].isin(["Benign", "Likely_benign", "Benign/Likely_benign"])))]

  # --- Define population groupings ---
  major_cols = ["AFR_AF", "EAS_AF", "AMR_AF", "FIN_AF", "NFE_AF", "SAS_AF"]
  minor_cols = ["MID_AF", "AMI_AF", "ASJ_AF", "OTH_AF"]
  cohort_col = "AF"
  freq_cols = major_cols + minor_cols + [cohort_col]
  gnomad_cols = major_cols + minor_cols

  # Identify SNVs vs indels
  is_snv = (df["REF"].str.len() == 1) & (df["ALT"].str.len() == 1)
  is_indel = ~is_snv

  # Identify positions containing multiple indel rows
  indel_counts = df.loc[is_indel].groupby(["CHROM", "POS"]).size()
  multiallelic_indel_positions = indel_counts[indel_counts > 1].index

  is_multiallelic_indel = is_indel & df.set_index(["CHROM", "POS"]).index.isin(multiallelic_indel_positions)

  # --- Require gnomAD frequencies for indels before converting missing values ---
  indel_has_gnomad = is_indel & (df[gnomad_cols] != ".").all(axis=1)

  # For SNVs, missing AF is treated as 0
  df[freq_cols] = df[freq_cols].replace(".", 0).apply(pd.to_numeric, errors="coerce")

  # --- Tier 1 thresholds ---
  mask_major_t1 = (df[major_cols] < 0.01).all(axis=1)
  mask_minor_t1 = (df[minor_cols] < 0.10).all(axis=1)
  mask_cohort_t1 = df[cohort_col] < 0.01
  mask_keep_t1 = mask_major_t1 & mask_minor_t1 & mask_cohort_t1

  # Exclude only multiallelic indels that are absent from gnomAD
  #exclude_indel = is_multiallelic_indel & ~indel_has_gnomad
  #mask_keep_t1 &= ~exclude_indel

  # Without protein_coding restriction
  df.loc[mask_keep_t1, "ID"].drop_duplicates().to_csv("~{gene}.rare_variants.txt", sep="\t", index=False, header=False)
  df.loc[mask_keep_t1, ["ID", "AC"]].to_csv("~{gene}.rare_variants_with_AC.txt", sep="\t", index=False)

  CODE
  >>>
  output {
    File out1 = "~{gene}.rare_variants.txt"
    File out2 = "~{gene}.rare_variants_with_AC.txt"
  }
  runtime {
    docker: 'vanallenlab/g2c_pipeline'
    preemptible: 3
  }
}

task T7_get_pgv_carriers {
    input {
        File vcf
    }

    command <<<
    set -euxo pipefail

    bcftools query -l ~{vcf} > all_samples.txt

    bcftools query -f '[%SAMPLE\t%GT\n]' -i 'GT="alt"' ~{vcf} \
        | cut -f1 \
        | sort -u > sample_ids.txt
    >>>

    output {
        File out1 = "sample_ids.txt"
    }

    runtime {
        docker: "vanallenlab/g2c_pipeline"
        preemptible: 3
    }
}
