# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0

workflow STEP_8_FILTER_VARIANTS {
  input {
    String step_4_output_dir       # Directory to STEP_4 Output VCFs

    #File rf_predictions # File output from hail RF in step 4
    Float tp_cutoff_snp	# Specify a tp_cutoff_snp based on results from step 5 through step 7
    Float tp_cutoff_indel # Specify a tp_cutoff_indel based on results from step 5 through step 7

    String storage_directory = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228" 
  }

  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call gather_vcfs {
    input:
      dir = step_4_output_dir
  }

  call sort_vcf_list {
    input:
      unsorted_vcf_list = gather_vcfs.vcf_list
  }

  call gather_tsvs{}
  Int negative_shards = 0
  scatter (i in range(length(gather_tsvs.out)-negative_shards)) { 
    # Create a list of variants (CHR_POS_ALT) that we will keep based off of rf_predictions
    call create_tp_bedFile {
      input:
        rf_predictions = gather_tsvs.out[i],
        tp_cutoff_snp = tp_cutoff_snp,
        tp_cutoff_indel = tp_cutoff_indel
    }
  }

  Array[File] shards = sort_vcf_list.vcf_arr

  # Loop through VCF shards to annotate them in parallel
  scatter (i in range(length(shards)-negative_shards)){
    # We annotate variant ID w/ CHR_POS_ALT
    call annotate_vcf {
      input:
        vcf = shards[i]
    }
    # We filter variants by ID using output from create_tp_bedFile task
    call filter_vcf {
      input:
        vcf = annotate_vcf.annotated_vcf,
        variant_list = create_tp_bedFile.bed_file[i]
    }

    if (defined(filter_vcf.filtered_vcf)){
      call copy_vcfs_to_storage {
        input:
          vcf = filter_vcf.filtered_vcf,
          vcf_idx = filter_vcf.filtered_vcf_idx, 
          storage_directory = storage_directory
      }
    }
  }

}


task copy_vcfs_to_storage {
  input {
    File? vcf
    File? vcf_idx
    String storage_directory
  }
  command <<<
  set -eu -o pipefail
  gsutil -m cp ~{vcf} gs://~{storage_directory}/STEP_8_FILTER_TO_TP_VARIANTS/sharded_vcfs/
  gsutil -m cp ~{vcf_idx} gs://~{storage_directory}/STEP_8_FILTER_TO_TP_VARIANTS/sharded_vcfs/
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task gather_tsvs {
  input {
  }
  command <<<
  python3 <<CODE
  prefix = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_4_RANDOM_FOREST/sharded_tsvs/rf_predictions."
  suffix = ".tsv.bgz"

  with open("rf_predictions.sharded.list", "w") as f:
    for i in range(3102):  # 0 to 3102 inclusive
      f.write(f"{prefix}{i}{suffix}\n")
  CODE
  >>>
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
  output {
    Array[String] out = read_lines("rf_predictions.sharded.list")
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


# We annotate variant ID w/ CHR_POS_ALT
task annotate_vcf {
  input {
    File vcf
  }
  String vcf_basename=basename("~{vcf}")
  command <<<
  bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' ~{vcf} -Oz -o ~{vcf_basename}
  >>>
  output {
    File annotated_vcf = "~{vcf_basename}"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
  }
}

# We filter variants by ID using output from create_tp_bedFile task
task filter_vcf {
  input {
    File vcf
    File variant_list
  }
  String vcf_basename=basename("~{vcf}")
  command <<<
  set -euxo pipefail
  if [ $(wc -l < ~{variant_list}) -eq 0 ]; then
    exit 0
  fi
  bcftools view --include ID==@~{variant_list} ~{vcf} -O z -o ~{vcf_basename}
  bcftools index -t ~{vcf_basename}
  >>>
  output {
    File? filtered_vcf = "~{vcf_basename}"
    File? filtered_vcf_idx = "~{vcf_basename}.tbi"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
  }
}


task SliceRemoteFiles {
  input {
    String vcf
    File vcf_idx
    String callset_name

    Int disk_gb = 50
    Float mem_gb = 12
    Int n_cpu = 4

    String chromosome = "chr22"

    String docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base:2023-07-28-v0.28.1-beta-e70dfbd7"
  }

  command <<<
    set -eu -o pipefail

    echo -e "\nSLICING QUERY REGIONS FROM VCF, REMOTELY\n"
    mv ~{vcf_idx} ./
    # Generate a BED file with intervals from 1 to 300,000,000 with steps of 1,000,000
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 300000000; i += 1000000) { print "~{chromosome}", i, (i + 1000000 - 1) } }' \
    | bgzip -c \
    > query.bed.gz

    echo -e "\nSLICING REMOTE ANNOTATION FILES:\n"
    mkdir remote_slices


    local_name=$( basename ~{vcf} )
    echo -e "$local_name"

    # Iterate over each interval and create a separate output file
    zcat query.bed.gz | while read -r chr start end; do
      interval="${chr}:${start}-${end}"
      output_file="remote_slices/${chr}_${start}_${end}.vcf.gz"
      echo -e "$output_file"
      export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
      tabix -h -R <(echo -e "${chr}\t${start}\t${end}") ~{vcf} | bgzip -c > "${output_file}"

      if [ $(bcftools view -H "${output_file}"| wc -l)  -eq 0 ]; then
        rm "${output_file}"
      else
        tabix -s 1 -b 2 -e 2 -f  "${output_file}"
      fi

    done

  >>>

  output {
    Array[File] remote_slices = glob("remote_slices/*gz")
    Array[File] remote_slice_idxs = glob("remote_slices/*gz.tbi")
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}

# Take the Random Forest Table and Output the Location and Alleles of the True Positives into a bed file
task create_tp_bedFile {
  input {
    File rf_predictions
    Float tp_cutoff_snp
    Float tp_cutoff_indel
  }
  command <<<
  set -euxo pipefail
  zcat ~{rf_predictions} > rf_predictions.tsv

  if [ $(wc -l < rf_predictions.tsv) -eq 0 ]; then
    touch tp_variants.bed
    exit 0
  fi

  python3 <<CODE
  import pandas as pd
  import gc
  import re
  
  # Function to extract fp_prob and tp_prob from the string
  def extract_probs(prob_str):
      # Use regular expressions to extract the numeric values for 'False' and 'True'
      fp_match = re.search(r'"False","value":([0-9.]+)', prob_str)
      tp_match = re.search(r'"True","value":([0-9.]+)', prob_str)
    
      # Extract the values, default to None if not found
      fp_prob = float(fp_match.group(1)) if fp_match else None
      tp_prob = float(tp_match.group(1)) if tp_match else None
    
      return fp_prob, tp_prob

  input_file = 'rf_predictions.tsv'
  output_file = 'tp_variants.bed'

  # Read the file into a pandas DataFrame
  data = pd.read_csv(input_file, delimiter='\t')

  # Apply the function to each row in the 'rf_probability' column
  data['fp_prob'], data['tp_prob'] = zip(*data['rf_probability'].apply(extract_probs))

  # Filter rows where the last column "rf_prediction" is "TP"
  filtered_data = data[
          ((data['tp_prob'] >= ~{tp_cutoff_snp}) & (data['allele_type'] == "snp")) | 
          ((data['tp_prob'] >= ~{tp_cutoff_indel}) & (data['allele_type'].str.contains("indel")))
  ]

  # We no longer need the original data DataFrame, so we delete it
  del data
  gc.collect()

  # Extract locus and alternate allele
  # Assuming 'locus' is a column in the data and 'alleles' is a string representation of a list
  output_data = filtered_data[['locus', 'alleles']].copy()
  output_data['alt_allele'] = output_data['alleles'].apply(lambda x: x.strip('[]').split(',')[1].strip().strip('"'))
  output_data['ref_allele'] = output_data['alleles'].apply(lambda x: x.strip('[]').split(',')[0].strip().strip('"'))

  # We no longer need the filtered_data DataFrame, so we delete it
  del filtered_data
  gc.collect()

  # Assuming your DataFrame is loaded as df
  # Define a function to generate the desired format
  def format_locus_allele(row):
      chrom, pos = row['locus'].split(':')  # Split the locus column into chromosome and position
      return f"{chrom}_{pos}_{row['ref_allele']}_{row['alt_allele']}"  # Format the output as chr_pos_alt

  # Apply the function to each row and save it to a new column
  output_data['formatted'] = output_data.apply(format_locus_allele, axis=1)

  # Save to a file or print the formatted output
  with open(output_file, "w") as f:
      for formatted_value in output_data['formatted']:
          f.write(f"{formatted_value}\n")

  CODE
  >>>
  output {
    File bed_file = "tp_variants.bed"
  }
  runtime{
    docker:"vanallenlab/g2c_pipeline"
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
