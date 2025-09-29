# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0

workflow STEP_12_FILTER_COHORT {
  input {
    String complex_logic
    String lower_case_name = "hboc"
    Int max_age = 200
    Int min_age = 0
    String cohorts = "aou"
    File explained_samples = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/samples_with_pvs.list"
    File pca = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_11_GENETIC_RELATEDNESS/ufc.eigenvec"
    File kinship = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_11_GENETIC_RELATEDNESS/ufc.kin0"
    File patient_info = "gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz"
    File phenotype_data = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/dfci-ufc.aou.phenos.v2.tsv.gz"
    File subjects_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_4_RANDOM_FOREST/inputs/ufc_subjects.list"
    String use_original_dx = "True"
  }
  call T1_Filter {
    input:
      metadata = patient_info,
      phenotype_data = phenotype_data,
      sample_list = subjects_list,
      pca = pca,
      kinship = kinship,
      min_age = min_age,
      max_age = max_age,
      cohorts = cohorts,
      exclude_samples = explained_samples,
      outfile = lower_case_name + ".cohort.metadata",
      log_file = lower_case_name + ".log",
      use_original_dx = use_original_dx,
      complex_logic = complex_logic
  }
  call copy_to_storage {
    input:
      metadata = T1_Filter.out1,
      log_file = T1_Filter.out2
  }
}
task T1_Filter {
  input {
    File metadata
    File phenotype_data
    File sample_list
    File pca
    File kinship
    File exclude_samples
    Int min_age
    Int max_age
    String cohorts
    String use_original_dx
    String outfile
    String log_file
    String complex_logic
  }
  String metadata_basename = basename(metadata)
  String metadata_basename_unzipped = basename(metadata,".gz")
  String phenotype_basename = basename(phenotype_data)
  String phenotype_basename_unzipped = basename(phenotype_data, ".gz") 
  command <<<
  set -euxo pipefail
  mv ~{metadata} .
  gunzip ~{metadata_basename}
  mv ~{phenotype_data} .
  gunzip ~{phenotype_basename}

  /opt/conda/bin/python3 /opt/step_12.downsample_controls_per_disease.ufc.py \
    --metadata ~{metadata_basename_unzipped} \
    --phenotype-data ~{phenotype_basename_unzipped} \
    --sample-list ~{sample_list} \
    --pca ~{pca} \
    --kinship ~{kinship} \
    --min-age ~{min_age} \
    --max-age ~{max_age} \
    --exclude-samples ~{exclude_samples} \
    --cohorts ~{cohorts} \
    --outfile ~{outfile} \
    --log-file ~{log_file} \
    --use-original-dx ~{use_original_dx} \
    --complex-logic "~{complex_logic}"
  >>>
  runtime {
    docker: "vanallenlab/g2c_ufc"
    preemptible: 3
    disks: "local-disk 2 HDD"
  }
  output {
    File out1 = "~{outfile}"
    File out2 = "~{log_file}"
  }
}

task copy_to_storage {
  input {
    File metadata
    File log_file
  }
  String safe_cancer_type = basename(log_file, ".log")
  command <<<

  tail -n +2 ~{metadata} | cut -f2 | sort -u > "~{safe_cancer_type}.list"
  gsutil -m cp "~{safe_cancer_type}.list" \
          "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/~{safe_cancer_type}/"
  gsutil -m cp ~{metadata} \
          "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/~{safe_cancer_type}/~{safe_cancer_type}.metadata"
  gsutil -m cp ~{log_file} \
          "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/~{safe_cancer_type}/~{safe_cancer_type}.cohort.log"
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preepmtible: 3
    disks: "local-disk 2 HDD"
  }
}
