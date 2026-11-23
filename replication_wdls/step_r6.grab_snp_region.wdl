# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0

workflow STEP_R6_GRAB_SNP_REGION {
  input {
    Int start_sequence = 20945
    Int end_sequence = 20965
    String gene = "ZNF17"
    File relatedness = "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/relatedness_flagged_samples.tsv"
  }

  Array[Int] sequence_numbers = range(end_sequence - start_sequence + 1)

  scatter (offset in sequence_numbers) {
    Int sequence = start_sequence + offset

    call T1_process_vcf_shard {
      input:
        sequence = sequence,
        gene = gene,
        relatedness = relatedness
    }
  }

  call T2_concat_regional_vcfs {
    input:
      gene = gene,
      vcf_files = T1_process_vcf_shard.processed_vcf
  }

  output {
    File regional_vcf = T2_concat_regional_vcfs.regional_vcf
    File regional_vcf_index = T2_concat_regional_vcfs.regional_vcf_index
    Array[File] processed_shards = T1_process_vcf_shard.processed_vcf
  }
}


task T1_process_vcf_shard {
  input {
    Int sequence
    String gene
    File relatedness
  }

  String output_file_name = gene + "_region." + sequence + ".vcf.bgz"

  command <<<
    set -euo pipefail

    BASE="gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/vcf"
    RELATEDNESS="~{relatedness}"

    FILE=$(printf "%010d" ~{sequence})
    INPUT="${FILE}.vcf.bgz"
    OUTPUT_VCF="~{output_file_name}"

    echo "Downloading ${FILE}..."
    gsutil -q -u wb-cordial-diamond-9893 cp "$BASE/${FILE}.vcf.bgz" "$INPUT"

    echo "Processing ${FILE}..."
    echo "[$(TZ='America/New_York' date '+%Y-%m-%d %H:%M:%S %Z')] Processing ${FILE}..."

    bcftools view -S "^${RELATEDNESS}" --force-samples -Ou "$INPUT" |
    bcftools annotate -x FORMAT,^FORMAT/GT -Ou |
    bcftools norm -m -any -Ou |
    bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Oz --threads 1 -o "$OUTPUT_VCF"

    rm "$INPUT"
  >>>

  output {
    File processed_vcf = "~{output_file_name}"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bioinformatics_tools"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 32 HDD"
    preemptible: 3
  }
}


task T2_concat_regional_vcfs {
  input {
    String gene
    Array[File] vcf_files
  }

  String output_file_name = gene + "_region.vcf.bgz"

  command <<<
    set -euo pipefail

    OUTPUT_VCF="~{output_file_name}"

    bcftools concat -Oz -o "$OUTPUT_VCF" ~{sep=' ' vcf_files}
    bcftools index -t "$OUTPUT_VCF"
  >>>

  output {
    File regional_vcf = "~{output_file_name}"
    File regional_vcf_index = "~{output_file_name}.tbi"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bioinformatics_tools"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 50 HDD"
    preemptible: 3
  }
}

