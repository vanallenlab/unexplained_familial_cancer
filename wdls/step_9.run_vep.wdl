# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# WDL to run VEP on one or more input VCFs


version 1.0


#import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/Utilities.wdl" as Tasks
import "Vep_utilities/Utilities.wdl" as Tasks

workflow STEP_9_RUN_VEP {
  input {
    #Array[File] vcfs
    #Array[File] vcf_idxs

    File reference_fasta
    File vep_cache_tarball # VEP cache tarball downloaded from Ensembl
    Array[String?] gnomad_vcf_uris # URIs for gnomAD sites VCFs stored in Google buckets that can be remotely sliced using tabix for each shard
    Array[File?] gnomad_vcf_indexes # Indexes corresponding to gnomad_vcf_uris
    Array[String?] gnomad_infos # INFO keys to annotate from gnomad VCFs
    Array[String?] vep_remote_files # URIs for files stored in Google buckets that can be remotely sliced using tabix for each shard
    Array[File?] vep_remote_file_indexes # Indexes corresponding to vep_remote_files
    Array[File?] other_vep_files # All other files needed for VEP. These will be localized in full to each VM and moved to execution directory.

    Array[String] vep_options = [""]
    String vep_assembly = "GRCh38"
    Int vep_version = 110

    Boolean shard_vcfs = true
    Int records_per_shard = 50000
    Int remote_query_buffer = 2
    Boolean combine_output_vcfs = false
    String? cohort_prefix

    String bcftools_docker
    String vep_docker = "vanallenlab/g2c-vep:loftee_etc"

    String step_8_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_8_FILTER_TO_TP_VARIANTS/sharded_vcfs"
    String storage_directory = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
  }

  Int n_remote_files = length(select_all(vep_remote_files))
  Int n_gnomad_files = length(select_all(gnomad_vcf_uris))
  Boolean any_remote = n_remote_files + n_gnomad_files > 0 

  Int positive_shards = 1

  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call gather_vcfs {
    input:
      dir = step_8_output_dir,
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

    if (shard_vcfs) {
      call Tasks.ShardVcf {
        input:
          vcf = vcf,
          vcf_idx = vcf_idx,
          records_per_shard = records_per_shard,
          bcftools_docker = bcftools_docker
      }
    }
    Array[File] vcf_shards = select_first([ShardVcf.vcf_shards, [vcf]])
    Array[File] vcf_shard_idxs = select_first([ShardVcf.vcf_shard_idxs, [vcf_idx]])
    

    scatter ( shard_info in zip(vcf_shards, vcf_shard_idxs) ) {

      if (any_remote) {
        call SliceRemoteFiles {
          input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            gnomad_vcf_uris = gnomad_vcf_uris,
            gnomad_vcf_indexes = gnomad_vcf_indexes,
            gnomad_infos = gnomad_infos,
            vep_remote_files = vep_remote_files,
            vep_remote_file_indexes = vep_remote_file_indexes,
            query_buffer = remote_query_buffer
        }
      }

      Array[File?] all_other_vep_files = flatten(select_all([SliceRemoteFiles.remote_slices, 
                                                             SliceRemoteFiles.remote_slice_idxs, 
                                                             [SliceRemoteFiles.gnomad_vcf],
                                                             [SliceRemoteFiles.gnomad_vcf_idx],
                                                             other_vep_files]))
      if (defined(SliceRemoteFiles.gnomad_vcf)){
        call GrabLinesInGnomadShard {
          input:
            gnomad_vcf = SliceRemoteFiles.gnomad_vcf
        }
      }
      Int gnomad_variants = select_first([GrabLinesInGnomadShard.numlines,0])

      call RunVep {
        input:
          vcf = shard_info.left,
          vcf_idx = shard_info.right,
          reference_fasta = reference_fasta,
          vep_cache_tarball = vep_cache_tarball,
          other_vep_files = all_other_vep_files,
          gnomad_infos = select_all(gnomad_infos),
          vep_options = vep_options,
          vep_assembly = vep_assembly,
          vep_version = vep_version,
          docker = vep_docker,
          num_lines = gnomad_variants
      }
    }

    if (shard_vcfs) {
      call Tasks.ConcatVcfs as ConcatInnerShards {
        input:
          vcfs = RunVep.annotated_vcf,
          vcf_idxs = RunVep.annotated_vcf_idx,
          out_prefix = basename(vcf, ".vcf.bgz") + ".vep",
          bcftools_docker = bcftools_docker
      }
    }
    File vepped_vcf = select_first([ConcatInnerShards.merged_vcf, RunVep.annotated_vcf[0]])
    File vepped_vcf_idx = select_first([ConcatInnerShards.merged_vcf_idx, RunVep.annotated_vcf_idx[0]])

    call copy_vcfs_to_storage {
      input:
        vcf = vepped_vcf,
        vcf_idx = vepped_vcf_idx,
        storage_directory = storage_directory
    }
  }

  #if ( combine_output_vcfs ) {
  #  call Tasks.ConcatVcfs as ConcatOuterShards {
  #    input:
  #      vcfs = vepped_vcf,
  #      vcf_idxs = vepped_vcf_idx,
  #      out_prefix = select_first([cohort_prefix, "all_input_vcfs"]) + ".vep.merged",
  #      bcftools_docker = bcftools_docker
  #  }
  #}

  #output {
  #  Array[File] annotated_vcfs = vepped_vcf
  #  Array[File] annotated_vcf_idxs = vepped_vcf_idx
  #  File? combined_annotated_vcf = ConcatOuterShards.merged_vcf
  #  File? combined_annotated_vcf_idx = ConcatOuterShards.merged_vcf_idx
  #}
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

task copy_vcfs_to_storage {
  input {
    File vcf
    File vcf_idx
    String storage_directory
  }
  command <<<
  set -eu -o pipefail
  gsutil -m cp ~{vcf} gs://~{storage_directory}/STEP_9_RUN_VEP/sharded_vcfs/
  gsutil -m cp ~{vcf_idx} gs://~{storage_directory}/STEP_9_RUN_VEP/sharded_vcfs/
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task SliceRemoteFiles {
  input {
    File vcf
    File vcf_idx
    Array[String?] gnomad_vcf_uris
    Array[File?] gnomad_vcf_indexes
    Array[String?] vep_remote_files
    Array[File?] vep_remote_file_indexes
    Array[String?] gnomad_infos

    Int query_buffer = 1000
    String gnomad_outfile_prefix = "gnomad"

    Int disk_gb = 550
    Float mem_gb = 15.5
    Int n_cpu = 4

    String docker = "vanallenlab/g2c_pipeline:05aa88e"
    #String docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base:2023-07-28-v0.28.1-beta-e70dfbd7"
  }

  Int n_remote_files = length(select_all(vep_remote_files))
  Int n_gnomad_files = length(select_all(gnomad_vcf_uris))

  command <<<
    #set -euxo pipefail
    #export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
    #gcloud config get-value account > pet_account.txt

    set -euxo pipefail
    export CLOUDSDK_PYTHON=/opt/conda/envs/gatk-sv/bin/python3
    # Debugging helpers
    which python3
    which python
    ls -ltrha
    gcloud --version

    # Write access token to file to pick up failures
    gcloud auth application-default print-access-token > gcp_token.txt
    cat gcp_token.txt
    export GCS_OAUTH_TOKEN=$( cat gcp_token.txt)

    echo -e "\nSLICING QUERY REGIONS FROM VCF, REMOTELY\n"
    mv ~{vcf_idx} ./
    bcftools query --format '%CHROM\t%POS\n' ~{vcf} \
    | awk -v OFS="\t" -v buf=~{query_buffer} \
      '{ print $1, $2-buf, $2+buf }' \
    | awk -v OFS="\t" '{ if ($2<1) $2=1; print $1, $2, $3 }' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    | bgzip -c \
    > query.bed.gz
    n_queries=$( zcat query.bed.gz | wc -l )

    zcat query.bed.gz \
      | sed -E 's/^chr//' \
      | bgzip -c > splice_ai.query.bed.gz

    #if [ ~{n_remote_files} -gt 0 ]; then
    #  echo -e "\nSLICING REMOTE ANNOTATION FILES:\n"
    #  mkdir remote_slices
    #  mv ~{sep=" " select_all(vep_remote_file_indexes)} ./

    #  while read uri; do
    #    local_name=$( basename $uri )
    #    echo -e "$local_name"
    #    if [ $n_queries -gt 0 ]; then
    #      tabix -h -R query.bed.gz $uri \
    #      | bgzip -c > remote_slices/$local_name
    #    else
    #      tabix --only-header $uri \
    #      | bgzip -c > remote_slices/$local_name
    #    fi
    #    tabix -s 1 -b 2 -e 2 -f remote_slices/$local_name
    #  done < ~{write_lines(select_all(vep_remote_files))}
    #fi

  # New code to Slice VCF
  if [ ~{n_remote_files} -gt 0 ]; then
    echo -e "\nSLICING REMOTE ANNOTATION FILES:\n"
    mkdir remote_slices
    mv ~{sep=" " select_all(vep_remote_file_indexes)} ./

    while read uri; do
      local_name=$(basename "$uri")
      echo -e "$local_name"

      # detect if filename contains "spliceai" (case-insensitive)
      if echo "$local_name" | grep -qi "spliceai"; then
        query_bed="splice_ai.query.bed.gz"
      else
        query_bed="query.bed.gz"
      fi

      # slice VCF based on query BED
      if [ $n_queries -gt 0 ]; then
        tabix -h -R "$query_bed" "$uri" | bgzip -c > "remote_slices/$local_name"
      else
        tabix --only-header "$uri" | bgzip -c > "remote_slices/$local_name"
      fi

      # if this is a SpliceAI file, prepend "chr" to non-header lines before indexing
      if echo "$local_name" | grep -qi "spliceai"; then
        zcat "remote_slices/$local_name" \
          | awk 'BEGIN{OFS="\t"} /^#/ {print; next} { $1="chr"$1; print }' \
          | bgzip -c > "remote_slices/${local_name%.gz}.tmp.gz"
        mv "remote_slices/${local_name%.gz}.tmp.gz" "remote_slices/$local_name"
      fi

      # index final file
      tabix -s 1 -b 2 -e 2 -f "remote_slices/$local_name"

    done < ~{write_lines(select_all(vep_remote_files))}
  fi

    if [ ~{n_gnomad_files} -gt 0 ]; then
      echo -e "\nSLICING GNOMAD VCFs:\n"
      mkdir gnomad_slices/
      mv ~{sep=" " select_all(gnomad_vcf_indexes)} ./

      query_fmt="%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER"
      query_keys=$( cat ~{write_lines(select_all(gnomad_infos))} | awk '{ print $1"=%INFO/"$1 }' | paste -s -d\; )
      query_fmt="$query_fmt\t$query_keys\n"
      echo "Interpreted gnomAD query format as $query_fmt"

      if [ $n_queries -gt 0 ]; then
        while read uri; do
          local_name=$( basename $uri )
          echo -e "$local_name"
          tabix -H $uri | cut -f1-8 > "$local_name".header
          bcftools query \
            -R query.bed.gz \
            -f "$query_fmt" \
            $uri \
          | cat "$local_name".header - | bgzip -c \
          > gnomad_slices/"$local_name"
        done < ~{write_lines(select_all(gnomad_vcf_uris))}

        echo -e "\nMERGING GNOMAD VCFs:\n"
        bcftools concat --naive \
          -O z -o ~{gnomad_outfile_prefix}.vcf.gz \
          gnomad_slices/*.vcf.*z
      else
        tabix -H ~{select_all(gnomad_vcf_uris)[0]} \
        | bgzip -c > ~{gnomad_outfile_prefix}.vcf.gz
      fi
      tabix -p vcf -f ~{gnomad_outfile_prefix}.vcf.gz
    fi

  >>>

  output {
    Array[File?] remote_slices = glob("remote_slices/*gz")
    Array[File?] remote_slice_idxs = glob("remote_slices/*gz.tbi")
    File? gnomad_vcf = "~{gnomad_outfile_prefix}.vcf.gz"
    File? gnomad_vcf_idx = "~{gnomad_outfile_prefix}.vcf.gz.tbi"
    File? bed_file = "query.bed.gz"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    bootDiskSizeGb: 30
  }
}

task GrabLinesInGnomadShard {
  input {
    File? gnomad_vcf
  }
  command <<<
  set -euxo pipefail
  bcftools view -H ~{gnomad_vcf} | wc -l > numlines.txt 
  >>>
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
  }
  output {
    Int numlines = read_int("numlines.txt")
  }
}

task RunVep {
  input {
    File vcf
    File vcf_idx

    File reference_fasta
    File vep_cache_tarball
    Array[File?] other_vep_files
    String vep_assembly

    Array[String] vep_options
    Array[String] gnomad_infos = []
    Int vep_max_sv_size = 50
    Int vep_version = 110

    Int num_lines
    Float mem_gb = 7.5
    Int n_cpu = 4
    Int? disk_gb

    String docker
  }

  String out_filename = basename(vcf, ".vcf.bgz") + ".vep.vcf.bgz"
  Int default_disk_gb = ceil(10 * size([vcf, vep_cache_tarball, reference_fasta], "GB")) + 50

  command <<<
    set -euxo pipefail

    # Unpack contents of cache into $VEP_CACHE/
    # Note that $VEP_CACHE is a default ENV variable set in VEP docker
    tar -xzvf ~{vep_cache_tarball} -C $VEP_CACHE/

    # Relocate other_vep_files to execution directory
    if [ ~{defined(other_vep_files)} == "true" ]; then
      while read file; do
        mv $file ./
      done < ~{write_lines(select_all(other_vep_files))}
    fi

    # Unzip the LOFTEE directory
    #tar -xzf loftee.tar.gz

    # Build gnomad annotation command based on gnomad_infos
    gnomad_option=""
    if [ ~{length(gnomad_infos)} -gt 0 ]; then
      gnomad_option="--custom gnomad.vcf.gz,gnomAD,vcf,exact,0,~{sep=',' gnomad_infos}"
    fi
    clinvar_option="--custom file=clinvar.oct22_2025.chr_prefixed.vcf.gz,short_name=clinvar,fields=CLNSIG%CLNSIGCONF,format=vcf,type=exact"
    if [ ~{num_lines} -gt 0 ]; then
      # --species homo_sapiens \
      vep \
        --input_file ~{vcf} \
        --format vcf \
        --output_file ~{out_filename} \
        --vcf \
        --verbose \
        --compress_output bgzip \
        --force_overwrite \
        --assembly ~{vep_assembly} \
        --max_sv_size ~{vep_max_sv_size} \
        --offline \
        --cache \
        --dir_cache $VEP_CACHE/ \
        --cache_version ~{vep_version} \
        --dir_plugins $VEP_PLUGINS/ \
        --fasta ~{reference_fasta} \
        --minimal \
        --nearest gene \
        --distance 10000 \
        --numbers \
        --hgvs \
        --no_escape \
        --symbol \
        --canonical \
        --domains \
        ~{sep=" " vep_options} \
        $clinvar_option \
        $gnomad_option
    else
      vep \
        --input_file ~{vcf} \
        --format vcf \
        --output_file ~{out_filename} \
        --vcf \
        --verbose \
        --compress_output bgzip \
        --force_overwrite \
        --assembly ~{vep_assembly} \
        --max_sv_size ~{vep_max_sv_size} \
        --offline \
        --cache \
        --dir_cache $VEP_CACHE/ \
        --cache_version ~{vep_version} \
        --dir_plugins $VEP_PLUGINS/ \
        --fasta ~{reference_fasta} \
        --minimal \
        --nearest gene \
        --distance 10000 \
        --numbers \
        --hgvs \
        --no_escape \
        --symbol \
        --canonical \
        --domains \
        $clinvar_option \
        ~{sep=" " vep_options}
    fi
 
     tabix -f ~{out_filename}

  >>>

  output {
    File annotated_vcf = "~{out_filename}"
    File annotated_vcf_idx = "~{out_filename}.tbi"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    bootDiskSizeGb: 25
    preemptible: 3
  }
}
