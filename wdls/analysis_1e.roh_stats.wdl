# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_1C_LOGISTIC_REGRESSION {
  input { 
    String output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/stats/"
    File phenotype_data = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/dfci-ufc.aou.phenos.v2.tsv.gz"
    File sample_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cohorts/"
  }
  Int negative_shards = 0

  call T1_Concat_Files {
  }
  call T2_Get_Stats {
    input:
      genome_roh_file = T1_Concat_Files.out1,
      phenotype_data = phenotype_data
  }
}

task T1_Concat_Files {
  input {
  }

  command <<<
  set -euxo pipefail

  gsutil ls gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr{1..22}.txt \
  | gsutil cat - \
  | grep -v '^#' \
  > genome_roh.tsv

  >>>

  output {
    Array[File] out1 = "genome_roh.tsv.gz"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}


task T2_Get_Stats {
  input {
    File genome_roh_file
    File phenotype_data
    File sample_list
  }

  command <<<
  set -euxo pipefail
  /opt/conda/bin/python3 /opt/roh_tools.py \
    --phenotype-data ~{phenotype_basename_unzipped} \
    --samples ~{sample_list} \
    --roh-file ~{genome_roh_file}

  gzip roh_stats.zip roh_stats/
  >>>
  output {
    File out1 = "roh_stats.zip"
  }
  runtime {
    docker: "vanallenlab/g2c_ufc"
    preemptible: 3
  }
}