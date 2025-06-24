# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g
# Adapted from Ryan Collins code:
# https://github.com/vanallenlab/ped_germline_SV/blob/main/gatksv_scripts/get_kinship_and_pcs.py`

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow STEP_13_PATHOGENIC_SAMPLES {
  input {
    String step_9_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/sharded_vcfs" # Directoryto STEP_9 Output VCFs
    File cpg_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/germline_cpg_lists/gene_lists/cosmic_germline_genes.list"
    File samples_of_interest
  }

  Int negative_shards = 0
  call Tasks.gather_vcfs {
    input:
      dir = step_9_output_dir  
  }

  call Tasks.sort_vcf_list {
    input:
      unsorted_vcf_list =  gather_vcfs.vcf_list
  }
  
  scatter (i in range(length(sort_vcf_list.vcf_arr)-negative_shards)){
    call T1_Convert_To_TSV {
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        cpg_list = cpg_list,
        samples_of_interest = samples_of_interest
    }
  }

  call Tasks.concatenateFiles_noheader as concatenate_DOMINANTCPGs{
    input:
      callset_name = "dominant_samples",
      files = T1_Convert_To_TSV.out1
  }
  #call Tasks.concatenateFiles_noheader as concatenate_RECESSIVECPGs{
  #  input:
  #    callset_name = "recessive_samples",
  #    files = T1_Convert_To_TSV.out2
  #}
  #call Tasks.concatenateFiles_noheader as concatenate_OTHERCPGs{
  #  input:
  #    callset_name = "other_samples",
  #    files = T1_Convert_To_TSV.out3
  #}
}

task T1_Convert_To_TSV {
  input {
    File vcf
    File cpg_list
    File samples_of_interest
  }
  String output_file = basename(vcf, ".vcf.bgz") + ".tsv"
  command <<<
  set -euxo pipefail
  
  # Filter to just Samples of interest
  bcftools view -S ~{samples_of_interest} -Oz -o tmp1.vcf.gz ~{vcf} 
  rm ~{vcf}

  # Get the header started
  echo -e "CHROM\tPOS\tID\tREF\tALT\tIMPACT\tSYMBOL\tclinvar_clnsig\tgnomAD_AF_non_cancer_afr\tgnomAD_AF_non_cancer_amr\tgnomAD_AF_non_cancer_eas\tgnomAD_AF_non_cancer_fin\tgnomAD_AF_non_cancer_nfe\tgnomAD_AF_non_cancer_sas\tgnomAD_AF_non_cancer_ami\tgnomAD_AF_non_cancer_asj\tgnomAD_AF_non_cancer_mid\tgnomAD_AF_non_cancer_oth\tSAMPLES" > ~{output_file}

  # If gnomad annotations are not found, this is not a regions enriched w/ genes. So we are not concerned about it here.
  if ! grep -q '^##INFO=<ID=gnomAD_AF_non_cancer_afr' <(bcftools view -h tmp1.vcf.gz); then
      echo "Error: gnomAD annotations not found in VCF header." >&2
      touch dominant.tsv
      exit 0
  fi

  bcftools +split-vep tmp1.vcf.gz -i 'GT="alt" && AF <= 0.01' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%IMPACT\t%SYMBOL\t%clinvar_clnsig\t%gnomAD_AF_non_cancer_afr\t%gnomAD_AF_non_cancer_amr\t%gnomAD_AF_non_cancer_eas\t%gnomAD_AF_non_cancer_fin\t%gnomAD_AF_non_cancer_nfe\t%gnomAD_AF_non_cancer_sas\t%gnomAD_AF_non_cancer_ami\t%gnomAD_AF_non_cancer_asj\t%gnomAD_AF_non_cancer_mid\t%gnomAD_AF_non_cancer_oth\t[%SAMPLE,]\n' -d > tmp.txt
  rm tmp1.vcf.gz
  grep -E 'HIGH|Pathogenic|Likely_pathogenic' tmp.txt > tmp2.txt || touch tmp2.txt
  rm tmp.txt

  # Filter the File to only include CPGS
  awk 'NR==FNR {cpg[$1]; next} $7 in cpg' <(grep AD ~{cpg_list} | grep -v RAD | cut -f1 ) tmp2.txt > tmp3.txt

  rm tmp2.txt
  sort -u < tmp3.txt >> ~{output_file}


  python3 <<CODE
  import pandas as pd
  df = pd.read_csv("~{output_file}",sep='\t',index_col=False)
  df = df[df['clinvar_clnsig'].isin(["Pathogenic","Pathogenic/Likely_pathogenic", "Likely_pathogenic","."])]

  # Identify relevant columns
  gnomad_major_cols = ['gnomAD_AF_non_cancer_afr','gnomAD_AF_non_cancer_amr','gnomAD_AF_non_cancer_eas','gnomAD_AF_non_cancer_fin','gnomAD_AF_non_cancer_nfe','gnomAD_AF_non_cancer_sas']
  gnomad_minor_cols = ['gnomAD_AF_non_cancer_mid','gnomAD_AF_non_cancer_ami','gnomAD_AF_non_cancer_asj','gnomAD_AF_non_cancer_oth']
  gnomad_cols = gnomad_major_cols + gnomad_minor_cols

  # Replace '.' with 0 and convert to float
  df[gnomad_cols] = df[gnomad_cols].replace('.', 0).astype(float)

  # Build mask for gnomAD filter
  gnomad_ok = (df[gnomad_major_cols] <= 0.01).all(axis=1) & (df[gnomad_minor_cols] <= 0.1).all(axis=1)

  # Build mask for allowed ClinVar significance
  clinvar_ok = df["clinvar_clnsig"].isin(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"])

  # Combine logic: keep row if either condition is satisfied
  filtered_df = df[gnomad_ok | clinvar_ok]
  print("Check X")
  cpg_df = pd.read_csv("~{cpg_list}",sep='\t',index_col=False,header=None, names=['GENE','PHENOTYPE_RELATIONSHIP','Role_in_Cancer'])

  # Merge on gene symbol
  filtered_df = filtered_df.merge(cpg_df, left_on='SYMBOL', right_on='GENE', how='left')

  # Final TSG/ClinVar logic
  # Keep if: Role_in_Cancer includes 'tsg' OR clinvar_clnsig is not "."
  is_tsg = filtered_df['Role_in_Cancer'].str.contains('TSG')
  has_valid_clinvar = ~filtered_df['clinvar_clnsig'].isin(['.'])

  filtered_df = filtered_df[is_tsg | has_valid_clinvar]

  filtered_df.to_csv("dominant.tsv",sep='\t',index=False)
  CODE
  >>>
  output {
    File out1 = "dominant.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    disks: "local-disk 20 HDD"
    preemptible: 3
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

