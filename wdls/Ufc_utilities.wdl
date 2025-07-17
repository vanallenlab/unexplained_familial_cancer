# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Common WDL tasks shared across workflows


version 1.0

task list_files_from_directory {
  input {
    String dir
    String suffix = "*"
    String prefix = "*"
  }
  command <<<
  set -euxo pipefail
  gsutil ls ~{dir} | (
    if [[ ~{prefix} != "*" && ~{suffix} != "*" ]]; then
      grep -E "^.*/~{prefix}.*~{suffix}$"
    elif [[ ~{prefix} != "*" ]]; then
      grep -E "^.*/~{prefix}.*$"
    elif [[ ~{suffix} != "*" ]]; then
      grep -E "^.*/.*~{suffix}$"
    else
      cat
    fi
  ) | sort -u > sorted_files.list

  >>>
  output {
    Array[String] out1 = read_lines("sorted_files.list")
  }
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task gather_chromosome_level_vcfs {
  input {
    String dir
    String sex_chromosomes = "TRUE"
  }
  command <<<
  if [ "~{sex_chromosomes}" == "TRUE" ]; then
    for chr in {1..22} X Y; do
      gsutil ls "~{dir}/*.vcf.bgz" | grep -E "chr${chr}-" > "chr${chr}.list"
    done
  else
    for chr in {1..22}; do
      gsutil ls "~{dir}/*.vcf.bgz" | grep -E "chr${chr}-" > "chr${chr}.list"
    done
  fi
 
  >>>
  output {
    Array[File] out1 = glob("*.list")
  }
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task gather_a_chromosome_level_vcfs {
  input {
    String dir
    String chr_num
  }
  command <<<
  gsutil ls "~{dir}/*.vcf.bgz" | grep "chr~{chr_num}" > "chr~{chr_num}.list"

  >>>
  output {
    Array[File] out1 = glob("*.list")
  }
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task gather_chromosome22_level_vcfs {
  input {
    String dir
  }
  command <<<
  gsutil ls "~{dir}/*.vcf.bgz" | grep "chr22" > "chr${chr}.list"

  >>>
  output {
    Array[File] out1 = glob("*.list")
  }
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task gather_vcfs {
  input {
    String dir
  }
  command <<<
  set -eu -o pipefail
  gsutil ls ~{dir}/*.vcf.bgz > vcf.list
  >>>
  output {
    File vcf_list = "vcf.list"
  }
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}


task sort_vcf_list {
  input {
    File unsorted_vcf_list
  }
  command <<<
  set -eu -o pipefail
  python3 <<CODE
  import re

  # Input file containing the list of file paths
  file_path = '~{unsorted_vcf_list}'

  # Read and parse file paths
  with open(file_path, 'r') as f:
      paths = f.readlines()

  # Function to extract chromosome and shard number for sorting
  def extract_key(path):
      # Extract chromosome and shard numbers using regex
      match = re.search(r'chr([0-9XY]+)-finalrun\.([0-9]+)', path)
      if match:
          chrom_str, shard_str = match.groups()

          # Convert chromosome to integer, treating 'X' as 23 and 'Y' as 24
          chrom_num = 23 if chrom_str == 'X' else 24 if chrom_str == 'Y' else int(chrom_str)
          shard_num = int(shard_str)

          # Return a tuple with chromosome and shard for sorting
          return (chrom_num, shard_num)
      else:
          return (float('inf'), float('inf'))  # Unmatched lines go to the end
  # Sort paths using the extracted keys
  sorted_paths = sorted(paths, key=extract_key)


  # Write to vcf.sorted.list and vcf_idx.sorted.list
  with open('vcf.sorted.list', 'w') as vcf_file, open('vcf_idx.sorted.list', 'w') as vcf_idx_file:
      for path in sorted_paths:
          clean_path = path.strip()
          vcf_file.write(clean_path + '\n')
          vcf_idx_file.write(clean_path + '.tbi\n')
  CODE
  >>>
  output {
    Array[String] vcf_arr = read_lines("vcf.sorted.list")
    Array[String] vcf_idx_arr = read_lines("vcf_idx.sorted.list")
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task copy_file_to_storage {
  input {
    File text_file
    String output_dir
  }
  command <<<
  set -eu -o pipefail
  gsutil -m cp ~{text_file} ~{output_dir}
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task copy_vcfs_to_storage {
  input {
    File vcf
    File? vcf_idx
    String storage_directory
    String storage_subdirectory
  }
  command <<<
  set -eu -o pipefail
  gsutil -m cp ~{vcf} gs://~{storage_directory}/~{storage_subdirectory}/sharded_vcfs/
  if [ "~{vcf_idx}" != "" ]; then
    gsutil -m cp ~{vcf_idx} gs://~{storage_directory}/~{storage_subdirectory}/sharded_vcfs/
  fi
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}


task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String callset_name

    String bcftools_concat_options = ""

    Float mem_gb = 3.5
    Int cpu_cores = 2
    Int? disk_gb

    String bcftools_docker = "vanallenlab/bcftools"
  }

  String out_filename = callset_name + ".vcf.gz"

  Int default_disk_gb = ceil(2.5 * size(vcfs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools concat \
      ~{bcftools_concat_options} \
      --file-list ~{write_lines(vcfs)} \
      -O z \
      -o ~{out_filename} \
      --threads ~{cpu_cores}

    bcftools index -t ~{out_filename}
  >>>

  output {
    File merged_vcf = "~{out_filename}"
    File merged_vcf_idx = "~{out_filename}.tbi"
  }

  runtime {
    docker: bcftools_docker
    memory: mem_gb + " GB"
    cpu: cpu_cores
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
}

# We filter variants by ID using output from create_tp_bedFile task
task Filter_VCF_USING_ID {
  input {
    File vcf
    File variant_list  #This is to filter by ID: CHROM_POS_REF_ALT
  }
  String vcf_basename=basename("~{vcf}")
  command <<<
  set -euxo pipefail
  if [ $(wc -l < ~{variant_list}) -eq 0 ]; then
    echo "No Variants To Keep"
    exit 0
  fi
  bcftools view --include ID==@~{variant_list} ~{vcf} -O z -o ~{vcf_basename}
  bcftools index -t ~{vcf_basename}
  >>>
  output {
    File out1 = "~{vcf_basename}"
    File out2 = "~{vcf_basename}.tbi"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
  }
}

task Filter_VCF_USING_SAMPLES {
  input {
    File vcf
    File samples_list  #This is to filter by SAMPLES
  }
  String vcf_basename=basename(vcf)
  command <<<
  set -euxo pipefail
  if [ $(wc -l < ~{samples_list}) -eq 0 ]; then
    echo "No Subjects"
    exit 1
  fi
  bcftools view -S ~{samples_list} ~{vcf} -O z -o tmp1.vcf.gz
  bcftools view -i 'AC > 0' tmp1.vcf.gz -O z -o ~{vcf_basename}
  bcftools index -t ~{vcf_basename}
  >>>
  output {
    File out1 = "~{vcf_basename}"
    File out2 = "~{vcf_basename}.tbi"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
  }
}

task sum_tables_by_sample {
  input {
    Array[File] input_tables
    String output_filename = "summary_statistics.tsv" # default name
  }

  Int default_mem_gb = ceil(2.5 * size(input_tables, "GB")) + 4
  Int default_disk_gb = ceil(4 * size(input_tables, "GB")) + 4

  command <<<
  set -x pipefail

  python3 <<CODE
  import pandas as pd

  # Initialize empty master_df
  master_df = None

  # Loop through input files
  for file_path in "~{sep=' ' input_tables}".split():
      try:
        df = pd.read_csv(file_path, sep='\t',index_col=False)
      except Exception as e:
        print(f"{file_path} is empty: {e}")
        continue

      if master_df is None:
          master_df = df.copy()
      else:
          # Add numeric columns, keep 'Sample' untouched
          merge_cols = [col for col in df.columns if col != 'Sample']
          master_df[merge_cols] += df[merge_cols]

  # Save the final dataframe
  master_df.to_csv("~{output_filename}", sep='\t', index=False)
  CODE
  >>>

  output {
    File out1 = output_filename
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    memory: "~{default_mem_gb}GB"
    disks: "local-disk ~{default_disk_gb} HDD"
    preemptible: 3
  }
}


task append_covariates {
  input {
    File data
    File sample_data = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/dfci-ufc.sample_meta.gatkhc_posthoc_outliers.tsv"
    File continental_pcs = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_11_GENETIC_RELATEDNESS/ufc.eigenvec"
    String output_name
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  data = pd.read_csv("~{data}",index_col=False,sep='\t')
  pcs = pd.read_csv("~{continental_pcs}",sep='\t',index_col=False)
  sample_data = pd.read_csv("~{sample_data}",index_col=False,sep='\t')

  # Strip whitespace from key columns used in merging
  data["Sample"] = data["Sample"].astype(str).str.strip()
  sample_data["original_id"] = sample_data["original_id"].astype(str).str.strip()
  pcs["#IID"] = pcs["#IID"].astype(str).str.strip()

  merged_df = data.merge(sample_data,left_on="Sample", right_on="original_id")
  merged_df = merged_df.merge(pcs, left_on="Sample", right_on="#IID")
  merged_df.to_csv("~{output_name}",sep='\t',index=False)
  CODE
  >>>
  output {
    File out1 = "~{output_name}"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task concatenateFiles_noheader {
  input {
    Array[File] files
    String callset_name = "out"
  }

  command <<<
  # Put Header Down
  head -n -1 ~{files[0]} > ~{callset_name}.tsv

  # Add files
  for f in ~{sep=" " files}; do
    tail -n +2 "$f"
  done >> ~{callset_name}.tsv
  gzip ~{callset_name}.tsv
  >>>

  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }

  output {
    File out1 = "~{callset_name}.tsv.gz"
  }
}

task concatenateGzippedFiles_noheader {
  input {
    Array[File] files
    String callset_name = "out"
  }

  command <<<
  # Extract header from the first gzipped file
  zcat ~{files[0]} | head -n 1 > ~{callset_name}.tsv

  # Append data from all files, skipping the header line
  for f in ~{sep=" " files}; do
    zcat "$f" | tail -n +2
  done >> ~{callset_name}.tsv

  # Compress final file
  gzip ~{callset_name}.tsv
  >>>

  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }

  output {
    File out1 = "~{callset_name}.tsv.gz"
  }
}

task concatenateFiles {
  input {
    Array[File] files
    String output_name = "out"
  }

  command <<<
    cat ~{sep=' ' files} > ~{output_name}.tsv
    gzip ~{output_name}.tsv
  >>>

  output {
    File out1 = "~{output_name}.tsv.gz"
  }

  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }
} 
