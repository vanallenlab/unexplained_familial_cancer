# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# # Distributed under the terms of the GNU GPL v2.0
  
version 1.0 
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
  
workflow STEP_R1_TIER_VARIANTS {
  input {
    String CDRv9_cpg_vcfs = "gs://dataproc-staging-wb-cordial-diamond-9893/vep_annotated_cpgs_sept2026/CDRv9"
  } 
  

  Int positive_shards = 160

  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call gather_vcfs {
    input:
      dir = CDRv9_cpg_vcfs,
      positive_shards = positive_shards
  } 
    
  call sort_vcf_list {
    input:
      unsorted_vcf_list = gather_vcfs.vcf_list
  }

  Array[File] vcfs = sort_vcf_list.vcf_arr
  Array[File] vcf_idxs = sort_vcf_list.vcf_idx_arr

  scatter ( vcf_info in zip(vcfs, vcf_idxs) ) {
    File vcf = vcf_info.left
    File vcf_idx = vcf_info.right

    call annotate_vcf {
      input:
        vcf = vcf
    }

    call Tasks.copy_file_to_storage as copy1{
      input:
        text_file = annotate_vcf.out1,
        output_dir = "gs://dataproc-staging-wb-cordial-diamond-9893/vep_annotated_cpgs_sept2026/CDRv9/"
    }

    call Tasks.copy_file_to_storage as copy2{
      input:
        text_file = annotate_vcf.out2,
        output_dir = "gs://dataproc-staging-wb-cordial-diamond-9893/vep_annotated_cpgs_sept2026/CDRv9/"
    }
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
    File out1 = "vcf.sorted.list"
    File out2 = "vcf_idx.sorted.list"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task gather_vcfs {
  input {
    String dir
    Int positive_shards
  }
  command <<<
  gsutil ls ~{dir}/*.vcf.bgz | head -n ~{positive_shards} > vcf.list
  >>>
  output {
    File vcf_list = "vcf.list"
  }
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task annotate_vcf {
  input {
    File vcf
  }

  command <<<
    set -euxo pipefail

    bcftools annotate \
      --set-id '%CHROM\_%POS\_%REF\_%ALT' \
      -Oz \
      -o "~{basename(vcf, ".vep.vcf.bgz")}.vep.cdrv9.vcf.bgz" \
      "~{vcf}"

    bcftools index -t "~{basename(vcf, ".vep.vcf.bgz")}.vep.cdrv9.vcf.bgz"
  >>>

  output {
    File out1 = "~{basename(vcf, ".vep.vcf.bgz")}.vep.cdrv9.vcf.bgz"
    File out2 = "~{basename(vcf, ".vep.vcf.bgz")}.vep.cdrv9.vcf.bgz.tbi"
  }

  runtime {
    docker: "vanallenlab/bcftools"
  }
}
