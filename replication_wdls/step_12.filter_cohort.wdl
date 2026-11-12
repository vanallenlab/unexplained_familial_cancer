# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0

workflow STEP_12_FILTER_COHORT {
  input {
    #Array[String] cancer_types = ["Breast"] #["Bladder","Breast","Colorectal","Leukemia","Kidney","Lung","Melanoma","Non-Hodgkins","Ovary","Prostate","Thyroid","endometrial"]
    Array[String] cancer_types = ["breast"] #["bladder","breast","colorectal","leukemia","kidney","lung","melanoma","non-hodgkins","ovary","prostate","thyroid","endometrial"]
    Array[String] sex_arr = ["XX"] #["XX,XY","XX","XX,XY","XX,XY","XX,XY","XX,XY","XX,XY","XX,XY","XX","XY","XX,XY","XX"]
    Int familial = 0
    File explained_samples = "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/ufc_subjects.aou.pgv_carriers_in_154_cpgs.list"
    File pca = "gs://rw-migration-aou-rw-cf3e7774/STEP_11_GENETIC_RELATEDNESS/ufc.eigenvec"
    File kinship = "gs://rw-migration-aou-rw-cf3e7774/STEP_11_GENETIC_RELATEDNESS/ufc.kin0"
    File phenotype_data = "gs://rw-migration-aou-rw-84a0039b/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz"
    File metadata = "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/CDRv7.ufc_original.metadata"
    File subjects_list = "gs://rw-migration-aou-rw-cf3e7774/UFC_REFERENCE_FILES/cohorts/ufc_subjects.aou.list"
  }
  scatter (i in range(length(cancer_types))){
    String log_file = cancer_types[i] + ".log"
    String output_file = cancer_types[i] + ".cohort.metadata"  
    call T1_Filter {
      input:
        phenotype_data = phenotype_data,
        metadata = metadata,
        sample_list = subjects_list,
        cancer_type = cancer_types[i],
        pca = pca,
        sex_karyotypes = sex_arr[i],
        kinship = kinship,
        exclude_samples = explained_samples,
        familial = familial,
        outfile = output_file,
        log_file = log_file,
    }
    # call copy_to_storage {
    #   input:
    #     metadata = T1_Filter.out1,
    #     log_file = T1_Filter.out2,
    #     cancer_type = cancer_types[i]
    # }
  }
}
task T1_Filter {
  input {
    File metadata
    File phenotype_data
    File sample_list
    String cancer_type
    File pca
    File kinship
    File exclude_samples
    String? sex_karyotypes
    Int familial = 0
    String outfile
    String log_file
  }
  String metadata_basename = basename(metadata)
  String metadata_basename_unzipped = basename(metadata,".gz")
  String phenotype_basename = basename(phenotype_data)
  String phenotype_basename_unzipped = basename(phenotype_data, ".gz") 
  command <<<

  /opt/conda/bin/python3 /opt/step_12.downsample_controls_per_disease.revised_ufc.py \
    --phenotype-data ~{phenotype_data} \
    --metadata ~{metadata_basename_unzipped} \
    --sample-list ~{sample_list} \
    --cancer-type '~{cancer_type}' \
    --pca ~{pca} \
    --kinship ~{kinship} \
    --sex_karyotypes ~{sex_karyotypes} \
    --exclude-samples ~{exclude_samples} \
    --familial ~{familial} \
    --outfile ~{outfile} \
    --log-file ~{log_file}
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
    String cancer_type
  }
  String safe_cancer_type = basename(log_file, ".log")
  command <<<

  tail -n +2 ~{metadata} | cut -f2 | sort -u > "~{safe_cancer_type}.list"
  gsutil -m cp "~{cancer_type}.list" \
          "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/analysis/~{cancer_type}/~{cancer_type}.list"
  gsutil -m cp ~{metadata} \
          "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/analysis/~{cancer_type}/~{cancer_type}.metadata"
  gsutil -m cp ~{log_file} \
          "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/analysis/~{cancer_type}/~{cancer_type}.cohort.log"
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
    disks: "local-disk 2 HDD"
  }
}
