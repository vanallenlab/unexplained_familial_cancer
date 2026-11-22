# The Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow STEP_R6_GRAB_SNP_REGION {
  input {
    Int start_sequence = 7702
    Int end_sequence = 7720
    String gene = "ZNF346"
  }

  call T1_make_regional_vcf {
    input:
      start_sequence = start_sequence,
      end_sequence = end_sequence,
      gene = gene
  }
}
task T1_make_regional_vcf {
  input {
    Int start_sequence = 9299
    Int end_sequence = 9314
    String gene
  }

  String outdir = gene + "_region"
  String output_file_name = gene + "_region.vcf.bgz"

  command <<<
    set -euo pipefail

    BASE="gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/vcf"
    OUTDIR="~{outdir}"
    OUTPUT_VCF="~{output_file_name}"

    if (( ~{start_sequence} > ~{end_sequence} )); then
      echo "ERROR: start_sequence must be less than or equal to end_sequence" >&2
      exit 1
    fi

    mkdir -p "$OUTDIR"

    # Download each VCF and its tabix index.
    vcf_files=()
    for i in $(seq ~{start_sequence} ~{end_sequence}); do
      FILE=$(printf "%010d" "$i")

      echo "Downloading ${FILE}..."

      gsutil -q -u "$GOOGLE_CLOUD_PROJECT" cp \
        "$BASE/${FILE}.vcf.bgz" \
        "$OUTDIR/"

      gsutil -q -u "$GOOGLE_CLOUD_PROJECT" cp \
        "$BASE/${FILE}.vcf.bgz.tbi" \
        "$OUTDIR/"

      vcf_files+=("$OUTDIR/${FILE}.vcf.bgz")
    done

    # Concatenate VCFs in numeric sequence order.
    bcftools concat -Oz \
      -o "$OUTPUT_VCF" \
      "${vcf_files[@]}"

    # Index the concatenated VCF.
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
    preemptible: 3
  }
}

