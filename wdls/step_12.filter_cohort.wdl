### TEST

version 1.0

workflow STEP_12_FILTER_COHORT {
  input {
    String cancer_type
    Int max_age = 200
    Int min_age = 0
    String cohorts = "aou"
    String sex = "XX,XY"
    File explained_samples
    File pca
    File kinship
    File output_file = "output.tsv"
    File patient_info
    File phenotype_data
    File subjects_list
    String log_file
    String use_original_dx
  }
  call T1_Filter {
    input:
      metadata = patient_info,
      phenotype_data = phenotype_data,
      sample_list = subjects_list,
      cancer_type = cancer_type,
      pca = pca,
      sex_karyotypes = sex,
      kinship = kinship,
      min_age = min_age,
      max_age = max_age,
      cohorts = cohorts,
      exclude_samples = explained_samples,
      outfile = output_file,
      log_file = log_file,
      use_original_dx = use_original_dx
  }
  #call copy_to_storage {
  #  input:
  #    cancer_type = cancer_type,
  #    metadata = T1_Filter.out1,
  #    log_file = T1_Filter.out2
  #}
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
    Int min_age
    Int max_age
    String? apparent_aneuploidies
    String? sex_karyotypes
    String cohorts
    String use_original_dx
    String outfile
    String log_file
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
    --cancer-subtype ~{cancer_type} \
    --pca ~{pca} \
    --kinship ~{kinship} \
    --min-age ~{min_age} \
    --max-age ~{max_age} \
    --sex ~{sex_karyotypes} \
    --exclude-samples ~{exclude_samples} \
    --cohorts ~{cohorts} \
    --outfile ~{outfile} \
    --log-file ~{log_file} \
    --use-original-dx ~{use_original_dx}
  >>>
  runtime {
    docker: "vanallenlab/g2c_ufc"
    preemptible: 3
  }
  output {
    File out1 = "~{outfile}"
    File out2 = "~{log_file}"
  }
}

task copy_to_storage {
  input {
    String cancer_type
    File metadata
    File log_file
  }
  command <<<
  tail -n +2 ~{metadata} | cut -f2 | sort -u > ~{cancer_type}.list
  gsutil -m cp ~{cancer_type}.list \
          gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/~{cancer_type}/
  gsutil -m cp ~{metadata} \
          gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/~{cancer_type}/~{cancer_type}.metadata
  gsutil -m cp ~{log_file} \
          gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/~{cancer_type}/~{cancer_type}.cohort.log
  >>>
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preepmtible: 3
  }
}
task T1_Initial_Filter {
  input {
    String? cancer_type
    File patient_info
    File subjects_list
    File hard_remove_subjects_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_0_IMPORTANT_DATA/ufc_samples_with_plp_cpg.list"
    String? sex
    String? cohorts
    Int? min_age
    Int? max_age
  }
  command <<<
  gunzip -c ~{patient_info} > patient_data.tsv
  python3 <<CODE
  import pandas as pd
  df = pd.read_csv("patient_data.tsv",sep='\t',index_col=False)

  if "~{cohorts}".strip():
      df = df[df['cohort'].str.contains("~{cohorts}",case=False,na=False)]

  with open("~{subjects_list}", "r") as file:
      subjects_set = {line.strip() for line in file}  # Read each line, strip whitespace, and store in a set

  with open("~{hard_remove_subjects_list}", "r" as file:
      subjects_to_remove_set = {line.strip() for line in file}

  removed_subjects = subjects_set & subjects_to_remove_set
  percentage_removed = (len(removed_subjects) / len(subjects_set)) * 100

  subjects_set = subjects_set - subjects_to_remove_set

  df = df[df['original_id'].isin(subjects_set)]
  
  if "~{cancer_type}" != "pancancer":
      df = df[df['cancer'].str.contains("~{cancer_type}|control",case=False,na=False)]

  if "~{sex}".strip():
      df = df[df['inferred_sex'] == "~{sex}"]
  
  if "~{min_age}".strip():
      df = df[df['age'] >= int("~{min_age}")] 

  if "~{max_age}".strip():
      df = df[df['age'] <= int("~{max_age}")]

  df.to_csv("ufc.t1_initial_filter.tsv",sep='\t',index=False) 
  CODE
  >>>
  output {
    File out1 = "ufc.t1_initial_filter.tsv"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}
