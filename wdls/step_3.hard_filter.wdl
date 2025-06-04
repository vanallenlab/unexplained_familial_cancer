# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0

workflow STEP_3_HARD_FILTER {
  input {
    String joint_genotyping_output_dir
    String storage_directory = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228" 
  }
  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call gather_vcfs {
    input:
      dir = joint_genotyping_output_dir
  }
  call sort_vcf_list {
    input:
      unsorted_vcf_list = gather_vcfs.vcf_list
  }
  scatter (vcf_string in sort_vcf_list.vcf_arr){
    call normalize {
      input:
        vcf_pathway = vcf_string
    }
    call hard_filter {
      input:
        vcf = normalize.output_vcf
    }
    call copy_vcfs_to_storage {
      input: 
        vcf = hard_filter.output_vcf,
        vcf_index = hard_filter.output_vcf_idx,
        storage_directory = storage_directory
    }
  }

  call write_report {
    input:
      variants_raw = hard_filter.num_variants_prefilter,
      snps_raw = hard_filter.num_snps_prefilter,
      indels_raw = hard_filter.num_indels_prefilter,
      rare_variants_raw = hard_filter.num_rare_variants_prefilter,
      rare_snps_raw = hard_filter.num_rare_snps_prefilter,
      rare_indels_raw = hard_filter.num_rare_indels_prefilter,
      variants_filtered = hard_filter.num_variants_postfilter,
      snps_filtered = hard_filter.num_snps_postfilter,
      indels_filtered = hard_filter.num_indels_postfilter,
      rare_variants_filtered = hard_filter.num_rare_variants_postfilter,
      rare_snps_filtered = hard_filter.num_rare_snps_postfilter,
      rare_indels_filtered = hard_filter.num_rare_indels_postfilter
  }

  call copy_to_storage {
    input:
      summary_stats = write_report.summary_stats,
      storage_directory = storage_directory
  }
}

task copy_to_storage {
  input {
    File summary_stats

    String storage_directory
  }
  command <<<
  set -eu -o pipefail
  gsutil -m cp ~{summary_stats} gs://~{storage_directory}/STEP_3_HARD_FILTER/
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task copy_vcfs_to_storage {
  input {
    File vcf
    File vcf_index

    String storage_directory
  }
  command <<<
  set -eu -o pipefail
  gsutil -m cp ~{vcf_index} gs://~{storage_directory}/STEP_3_HARD_FILTER/sharded_vcfs/
  gsutil -m cp ~{vcf} gs://~{storage_directory}/STEP_3_HARD_FILTER/sharded_vcfs/
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
  }
}

task write_report {
  input {
    Array[Int] variants_raw
    Array[Int] snps_raw
    Array[Int] indels_raw
    Array[Int] rare_variants_raw
    Array[Int] rare_snps_raw
    Array[Int] rare_indels_raw
    Array[Int] variants_filtered
    Array[Int] snps_filtered
    Array[Int] indels_filtered
    Array[Int] rare_variants_filtered
    Array[Int] rare_snps_filtered
    Array[Int] rare_indels_filtered
  }
  command <<<
  sum_variants_raw=$((~{sep='+' variants_raw}))
  sum_snps_raw=$((~{sep='+' snps_raw}))
  sum_indels_raw=$((~{sep='+' indels_raw}))
  sum_rare_variants_raw=$((~{sep='+' rare_variants_raw}))
  sum_rare_snps_raw=$((~{sep='+' rare_snps_raw}))
  sum_rare_indels_raw=$((~{sep='+' rare_indels_raw}))
  sum_variants_filtered=$((~{sep='+' variants_filtered}))
  sum_snps_filtered=$((~{sep='+' snps_filtered}))
  sum_indels_filtered=$((~{sep='+' indels_filtered}))
  sum_rare_variants_filtered=$((~{sep='+' rare_variants_filtered}))
  sum_rare_snps_filtered=$((~{sep='+' rare_snps_filtered}))
  sum_rare_indels_filtered=$((~{sep='+' rare_indels_filtered}))

  echo -e "Filter_Status\tVariants\tRare_Variants\tSNPs\tRare_SNPs\tIndels\tRare_Indels" > step_3_stats.tsv
  echo -e "Pre-hard_filter\t${sum_variants_raw}\t${sum_rare_variants_raw}\t${sum_snps_raw}\t${sum_rare_snps_raw}\t${sum_indels_raw}\t${sum_rare_indels_raw}" >> step_3_stats.tsv
  echo -e "Post-hard_filter\t${sum_variants_filtered}\t${sum_rare_variants_filtered}\t${sum_snps_filtered}\t${sum_rare_snps_filtered}\t${sum_indels_filtered}\t${sum_rare_indels_filtered}" >> step_3_stats.tsv

  # Calculate the percentages for each category
  percent_variants_filtered=$((100 - (100 * sum_variants_filtered / sum_variants_raw)))
  percent_rare_variants_filtered=$((100 - (100 * sum_rare_variants_filtered / sum_rare_variants_raw)))
  percent_snps_filtered=$((100 - (100 * sum_snps_filtered / sum_snps_raw)))
  percent_rare_snps_filtered=$((100 - (100 * sum_rare_snps_filtered / sum_rare_snps_raw)))
  percent_indels_filtered=$((100 - (100 * sum_indels_filtered / sum_indels_raw)))
  percent_rare_indels_filtered=$((100 - (100 * sum_rare_indels_filtered / sum_rare_indels_raw)))

  # Append the "Percent_filtered" row to the file
  echo -e "Percent_filtered\t${percent_variants_filtered}%\t${percent_rare_variants_filtered}%\t${percent_snps_filtered}%\t${percent_rare_snps_filtered}%\t${percent_indels_filtered}%\t${percent_rare_indels_filtered}%" >> step_3_stats.tsv
  >>>
  output {
    File summary_stats = "step_3_stats.tsv"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task gather_vcfs {
  input {
    String dir
  }
  command <<<
  set -eu -o pipefail
  gsutil ls ~{dir}/*/variant_filtered_vcf/*.vcf.gz > vcf.list
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

task normalize {
  input {
    File vcf_pathway
  }

  String vcf_basename=basename("~{vcf_pathway}")
  Int default_mem_gb = ceil(2 * size(vcf_pathway,"GB")) + 5

  command <<<
  set -eu -o pipefail
  # Normalize VCF
  bcftools norm -m -any ~{vcf_pathway} -Oz -o shard.norm.vcf.gz
  rm ~{vcf_pathway}
  
  # Annotate w/ VAF and VAF1
  bcftools +fill-tags shard.norm.vcf.gz -Oz -o ~{vcf_basename} -- -t VAF,VAF1
  >>>
  output {
    File output_vcf = "~{vcf_basename}"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    memory: "~{default_mem_gb} GB"
    preemptible: 3
  }
}

task hard_filter {
  input {
    File vcf
  }

  String vcf_basename=basename("~{vcf}")
  Int default_mem_gb = ceil(2 * size(vcf,"GB")) + 5

  command <<<
  set -eu -o pipefail

  # Collect How many variants there are for 1) the whole vcf, 2) rare variants, 3) snps, 4) rare snps, 5) indels, 6) rare_indels
  bcftools view -H ~{vcf} | wc -l > prefilter_variants.txt
  bcftools view -H -i 'AF<0.01' ~{vcf} | wc -l > prefilter_rare_variants.txt
  bcftools view -H -v snps ~{vcf} | wc -l > prefilter_snps.txt
  bcftools view -H -v indels ~{vcf} | wc -l > prefilter_indels.txt
  bcftools view -H -v snps -i 'AF<0.01' ~{vcf} | wc -l > prefilter_rare_snps.txt
  bcftools view -H -v indels -i 'AF<0.01' ~{vcf} | wc -l > prefilter_rare_indels.txt

  echo "Num of variants before filters"
  cat prefilter_variants.txt

  bcftools view -i 'FORMAT/DP>=10 & FORMAT/GQ >=20 & VAF >= 0.2' ~{vcf} -Oz -o ~{vcf_basename}

  # Collect How many variants there are for 1) the whole vcf, 2) rare variants, 3) snps, 4) rare snps, 5) indels, 6) rare_indels
  echo "Num of variants after filters"
  bcftools view -H ~{vcf_basename} | wc -l > postfilter_variants.txt
  bcftools view -H -i 'AF<0.01' ~{vcf_basename} | wc -l > postfilter_rare_variants.txt
  bcftools view -H -v snps ~{vcf_basename} | wc -l > postfilter_snps.txt
  bcftools view -H -v indels ~{vcf_basename} | wc -l > postfilter_indels.txt
  bcftools view -H -v snps -i 'AF<0.01' ~{vcf_basename} | wc -l > postfilter_rare_snps.txt
  bcftools view -H -v indels -i 'AF<0.01' ~{vcf_basename} | wc -l > postfilter_rare_indels.txt

  echo "Num of variants after filters"
  cat postfilter_variants.txt

  bcftools index -t ~{vcf_basename}
  >>>
  output {
    Int num_variants_prefilter = read_int("prefilter_variants.txt")
    Int num_variants_postfilter = read_int("postfilter_variants.txt")

    Int num_rare_variants_prefilter = read_int("prefilter_rare_variants.txt")
    Int num_rare_variants_postfilter = read_int("postfilter_rare_variants.txt")

    Int num_snps_prefilter = read_int("prefilter_snps.txt")
    Int num_snps_postfilter = read_int("postfilter_snps.txt")

    Int num_rare_snps_prefilter = read_int("prefilter_rare_snps.txt")
    Int num_rare_snps_postfilter = read_int("postfilter_rare_snps.txt")

    Int num_indels_prefilter = read_int("prefilter_indels.txt")
    Int num_indels_postfilter = read_int("postfilter_indels.txt")

    Int num_rare_indels_prefilter = read_int("prefilter_rare_indels.txt")
    Int num_rare_indels_postfilter = read_int("postfilter_rare_indels.txt")

    File output_vcf = "~{vcf_basename}"
    File output_vcf_idx = "~{vcf_basename}.tbi"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
    memory: "~{default_mem_gb} GB"
  }
}

task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String out_prefix

    String bcftools_concat_options = ""

    Float? input_mem_gb
    Int cpu_cores = 2
    Int? disk_gb

  }

  String tmp_filename = out_prefix + ".vcf.gz"
  String out_filename = out_prefix + ".norm.vcf.gz"

  Int default_disk_gb = ceil(2.5 * size(vcfs, "GB")) + 10
  Int default_mem_gb = ceil(1.5 * size(vcfs,"GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools concat \
      ~{bcftools_concat_options} \
      ~{sep=' ' vcfs} \
      -O z \
      -o ~{tmp_filename} \
      --threads ~{cpu_cores}

    bcftools index -t ~{tmp_filename}

    bcftools norm -m -any ~{tmp_filename} -Oz -o ~{out_filename}
    bcftools index -t ~{out_filename}
  >>>

  output {
    File merged_vcf = "~{out_filename}"
    File merged_vcf_idx = "~{out_filename}.tbi"
  }

  runtime {
    docker: "vanallenlab/bcftools"
    memory: select_first([input_mem_gb,default_mem_gb]) + " GB"
    cpu: cpu_cores
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
}
