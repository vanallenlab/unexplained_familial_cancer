# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g
# Adapted from Ryan Collins code:
# https://github.com/vanallenlab/ped_germline_SV/blob/main/gatksv_scripts/get_kinship_and_pcs.py`

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_0_PATHOGENIC_SAMPLES {
  input {
    String step_9_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/sharded_vcfs" # Directoryto STEP_9 Output VCFs
    File cpg_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/riaz_genes.list"
    File samples_of_interest = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cohorts/ufc_subjects.aou.list"
  }

  Int negative_shards = 0
  Int positive_shards = 260
  call Tasks.gather_positive_vcfs {
    input:
      dir = step_9_output_dir,
      positive_shards = positive_shards  
  }

  call Tasks.sort_vcf_list {
    input:
      unsorted_vcf_list =  gather_positive_vcfs.vcf_list
  }
  
  scatter (i in range(length(sort_vcf_list.vcf_arr)-negative_shards)){
    call T1_Convert_To_TSV {
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        cpg_list = cpg_list,
        samples_of_interest = samples_of_interest
    }
  }

  call Tasks.concatenateFiles_noheader as concatenate_CPGs{
    input:
      callset_name = "missed_samples",
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
    File cpg_bed_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/riaz_genes.coordinates.bed"
    File tier1_variants = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier1.tsv"
  }
  String output_file = basename(vcf, ".vcf.bgz") + ".tsv"
  command <<<
  set -euxo pipefail
  
  # Filter to just Samples of interest
  bcftools view -S ~{samples_of_interest} -Oz -o tmp1.vcf.gz ~{vcf}
  bcftools index -t tmp1.vcf.gz
  bcftools view -R ~{cpg_bed_file} tmp1.vcf.gz -Oz -o tmp2.vcf.gz 
  bcftools view --include ID==@~{tier1_variants} tmp2.vcf.gz -O z -o tmp3.vcf.gz
 
  # Get the header started
  echo -e "CHROM\tPOS\tID\tREF\tALT\tIMPACT\tSYMBOL\tclinvar_clnsig\tBiotype\tConsequence\tSAMPLES" > ~{output_file}


  bcftools +split-vep tmp3.vcf.gz -i 'GT="alt"' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%IMPACT\t%SYMBOL\t%clinvar_CLNSIG\t%BIOTYPE\t%Consequence\t[%SAMPLE,]\n' -d > tmp.txt
  rm tmp1.vcf.gz

  # Filter the File to only include Autosomal Dominant CPGS
  awk 'NR==FNR {cpg[$1]; next} $7 in cpg' <(cut -f1 ~{cpg_list}) tmp.txt > tmp3.txt || touch tmp3.txt

  sort -u < tmp3.txt >> ~{output_file}

  python3 <<CODE
  import pandas as pd
  df = pd.read_csv("~{output_file}",sep='\t',index=False)
  df = df[(df['IMPACT'] == "MODERATE") | (df['IMPACT'] == "HIGH")]
  >>>
  output {
    File out1 = "~{output_file}"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    disks: "local-disk 20 HDD"
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

