# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0

workflow STEP_6_VARIANT_SENSITIVITY {
  input {
    File rf_predictions	# File output from hail RF in step 5

    # VCF files for HQ sites as designated by 1000G. Used for Sensitivity Analysis
    File ThousandGenomesSNPsVCF
    File ThousandGenomesSNPsVCF_idx
    File ThousandGenomesINDELsVCF
    File ThousandGenomesINDELsVCF_idx

  }
  call SliceRemoteFiles as SliceRemoteFiles_snps{
    input:
      vcf = ThousandGenomesSNPsVCF,
      vcf_idx = ThousandGenomesSNPsVCF_idx,
      callset_name = "1000G_snps"
  }
  call SliceRemoteFiles as SliceRemoteFiles_indels{
    input:
      vcf = ThousandGenomesINDELsVCF,
      vcf_idx = ThousandGenomesINDELsVCF_idx,
      callset_name = "1000G_indels"
  }
  scatter (i in range(length(SliceRemoteFiles_snps.remote_slices)-0)){
    call Extract_Variants_From_VCF as Extract_SNPs_From_VCF{
      input:
        vcf = SliceRemoteFiles_snps.remote_slices[i]
    }
  }
  scatter (idx in range(length(SliceRemoteFiles_indels.remote_slices)-0)){
    call Extract_Variants_From_VCF as Extract_INDELs_From_VCF{
      input:
        vcf = SliceRemoteFiles_indels.remote_slices[idx]
    }
  }
  call split_tsv_chunks{
    input:
      tsv_file = rf_predictions
  }
  
  # 1000Genome Data
  call concatenateFiles as concatenateFiles_snps{
    input:
      files = Extract_SNPs_From_VCF.variants,
      callset_name = "1000G_snps"
  }
  scatter( i in range(length(split_tsv_chunks.chunk_files))){
  #scatter( i in range(length(split_tsv_chunks.out1))){
    call Probabilities as SNP_Probabilities{
      input:
        variants = concatenateFiles_snps.out,
        rf_predictions = split_tsv_chunks.chunk_files[i],
        snp_or_indel = "snp"
    }
  }
  call concat_and_weight {
    input:
      input_files = SNP_Probabilities.tp_prob_output,
  }
  call concatenateFiles as concatenateFiles_indels{
    input:
      files = Extract_INDELs_From_VCF.variants,
      callset_name = "1000G_indels"
  }
  call Probabilities as INDEL_Probabilities{
    input:
      variants = concatenateFiles_indels.out,
      rf_predictions = rf_predictions,
      snp_or_indel = "indel"
  }

  # Plotting Sensitivity
  call plot_sensitivity as plot_sensitivity_snps{
    input:
      tsv_file = concat_and_weight.out1,
      snp_or_indel = "snp",
      max_sensitivity = "0.81"
  }
  call plot_sensitivity as plot_sensitivity_indels{
    input:
      tsv_file = INDEL_Probabilities.tp_prob_output,
      snp_or_indel = "indel",
      max_sensitivity = INDEL_Probabilities.max_sensitivity
  }
}

task plot_sensitivity {
  input {
    File tsv_file
    String snp_or_indel
    Float max_sensitivity
  }

  command <<<
  # Exit immediately if a command exits with a non-zero status
  set -eu -o pipefail

  python3 <<CODE
  import pandas as pd
  import matplotlib.pyplot as plt
  import numpy as np

  # Load the TSV file
  df = pd.read_csv("~{tsv_file}", sep='\t')

  # Sort by tp_prob in ascending order
  df = df.sort_values(by='tp_prob', ascending=True)

  # Define probability thresholds from 0 to 1 in increments of 0.01
  thresholds = np.arange(0, 1.01, 0.01)

  # Initialize list to store sensitivity at each threshold
  sensitivity = []

  # Calculate sensitivity at each threshold
  total_positives = len(df)

  for threshold in thresholds:
      # Count how many rows have tp_prob >= current threshold
      true_positives_at_threshold = df[df['tp_prob'] >= threshold].shape[0]
      # Sensitivity is TP / Total positives
      sensitivity_at_threshold = (true_positives_at_threshold / total_positives) * ~{max_sensitivity}
      sensitivity.append(sensitivity_at_threshold)

  # Plot the sensitivity curve
  plt.figure(figsize=(8, 6))
  plt.plot(thresholds, sensitivity, marker='o', linestyle='-', color='b')
  plt.xlabel('tp_prob')
  plt.ylabel('Sensitivity')
  plt.title('Sensitivity Curve for 1000G ~{snp_or_indel}')
  plt.grid(True)

  # Save the plot
  plt.savefig("sensitivity_plot_~{snp_or_indel}.png")

  # Create a DataFrame with thresholds and sensitivity
  df_sensitivity = pd.DataFrame({
      'tp_prob': thresholds,
      'sensitivity': sensitivity
  })

  # Save the DataFrame to a TSV file
  df_sensitivity.to_csv("sensitivity_output_~{snp_or_indel}.tsv", sep='\t', index=False)
  CODE
  >>>

  output {
    File sensitivity_plot = "sensitivity_plot_~{snp_or_indel}.png"
    File sensitivity_tsv = "sensitivity_output_~{snp_or_indel}.tsv"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "16G"
    disks: "local-disk 50 HDD"
  }
}


task Probabilities {
    input {
      File variants
      File rf_predictions
      String snp_or_indel
    }

  command <<<
  set -eu -o pipefail

  # Filter rf_predictions to be just snps or indels
  head -n 1 ~{rf_predictions} > filtered_predictions.tsv
  grep ~{snp_or_indel} ~{rf_predictions} > tmp.tsv
  echo "Number of Variants:"
  echo "$(wc -l < tmp.tsv)"
  echo "Number of testing variants"
  echo "$(wc -l ~{variants})"

  # Filter tmp.tsv so we just have snps or indels in the HQ vcf sites
  split -l 50000 ~{variants} snp_chunk_
  for f in snp_chunk_*; do
      grep -Fwf "$f" tmp.tsv >> filtered_predictions.tsv || true
  done

  echo "Number of Matched Variants:"
  echo "$(wc -l < filtered_predictions.tsv)"

  # Count the total number of variants in variants.txt
  total_variants=$(wc -l < ~{variants})

  # Count how many variants were found in tmp.tsv
  matched_variants=$( tail -n +2 filtered_predictions.tsv | wc -l)

  # Calculate the percentage of matched variants
  percentage=$(awk "BEGIN {printf \"%.2f\", ($matched_variants / $total_variants) * 100}")

  # Write the percentage to a file
  echo "$percentage" > matched_percentage.txt

  python3 <<CODE
  import pandas as pd
  import re

  # Load the two TSV files
  df = pd.read_csv("filtered_predictions.tsv", sep='\t')

  # Function to extract fp_prob and tp_prob from the string
  def extract_probs(prob_str):
      # Use regular expressions to extract the numeric values for 'False' and 'True'
      fp_match = re.search(r'"False","value":([0-9.]+)', prob_str)
      tp_match = re.search(r'"True","value":([0-9.]+)', prob_str)

      # Extract the values, default to None if not found
      fp_prob = float(fp_match.group(1)) if fp_match else None
      tp_prob = float(tp_match.group(1)) if tp_match else None

      return fp_prob, tp_prob

  # Apply the function to each row in the 'rf_probability' column
  df['fp_prob'], df['tp_prob'] = zip(*df['rf_probability'].apply(extract_probs))

  # Ensure that 'tp_prob' is a float (it should already be a float, but just to be safe)
  df['tp_prob'] = df['tp_prob'].astype(float)

  # Sort the DataFrame by 'tp_prob' in ascending order
  df = df.sort_values(by='tp_prob', ascending=True)

  # Save as a TSV file
  df.to_csv('tp_prob.tsv', sep='\t', index=False)
  CODE
  >>>

  output {
      File tp_prob_output = "tp_prob.tsv"
      Float max_sensitivity = read_float("matched_percentage.txt")
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    disks: "local-disk 100 HDD"
    memory: "8G"
  }
}

task concat_and_weight {
  input {
    Array[File] input_files   # TSVs with header
  }

  command <<<
    set -euo pipefail

    out_file="concatenated.tsv"
    header_written=0

    i=0
    for f in ~{sep=" " input_files}; do

      # Count non-header lines
      n_lines=$(($(wc -l < "$f") - 1))

      # Concatenate (keep header only once)
      if [ $header_written -eq 0 ]; then
        cat "$f" >> "$out_file"
        header_written=1
      else
        tail -n +2 "$f" >> "$out_file"
      fi

      # Save pair (weight, value)
      i=$((i+1))
    done

  >>>

  output {
    File out1 = "concatenated.tsv"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "4G"
    disks: "local-disk 100 HDD"
    preemptible: 3
  }
}

task Extract_Variants_From_VCF {
  input{
    File vcf
  }
  command <<<
   bcftools query -f '%CHROM:%POS\t\["%REF","%ALT"\]\n' ~{vcf} > variants.tsv 
  >>>
  output {
    File variants = "variants.tsv"
  }
  runtime {
    docker: "vanallenlab/bcftools"
  }
}

task SliceRemoteFiles {
  input {
    File vcf
    File vcf_idx
    String callset_name

    Int disk_gb = 50
    Float mem_gb = 12
    Int n_cpu = 4

    String docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base:2023-07-28-v0.28.1-beta-e70dfbd7"
  }

  command <<<
    set -euxo pipefail

    echo -e "\nSLICING QUERY REGIONS FROM VCF, REMOTELY\n"
    #mv ~{vcf_idx} ./
    # Generate a BED file with intervals from 1 to 300,000,000 with steps of 1,000,000
    (
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 280000000; i += 1000000) { print "chr1", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 280000000; i += 1000000) { print "chr2", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 230000000; i += 1000000) { print "chr3", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 220000000; i += 1000000) { print "chr4", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 210000000; i += 1000000) { print "chr5", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 200000000; i += 1000000) { print "chr6", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 190000000; i += 1000000) { print "chr7", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 180000000; i += 1000000) { print "chr8", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 170000000; i += 1000000) { print "chr9", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 160000000; i += 1000000) { print "chr10", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 160000000; i += 1000000) { print "chr11", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 160000000; i += 1000000) { print "chr12", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 140000000; i += 1000000) { print "chr13", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 140000000; i += 1000000) { print "chr14", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 130000000; i += 1000000) { print "chr15", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 120000000; i += 1000000) { print "chr16", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 110000000; i += 1000000) { print "chr17", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 1100000000; i += 1000000) { print "chr18", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 90000000; i += 1000000) { print "chr19", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 90000000; i += 1000000) { print "chr20", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 70000000; i += 1000000) { print "chr21", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 70000000; i += 1000000) { print "chr22", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 180000000; i += 1000000) { print "chrX", i, (i + 1000000 - 1) } }' ;
    awk -v OFS="\t" 'BEGIN { for (i = 1; i <= 90000000; i += 1000000) { print "chrY", i, (i + 1000000 - 1) } }' ;
    ) \
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
      export GCS_OAUTH_TOKEN='gcloud auth application-default print-access-token'
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

task concatenateFiles {
  input {
    Array[File] files
    String callset_name
  }
  command <<<
    cat ~{sep=" " files} >> ~{callset_name}.tsv
  >>>
  runtime{
    docker:"ubuntu:latest"
    preemptible: 3
  }
  output{
    File out = "~{callset_name}.tsv"
  }
}

task split_tsv_chunks {
  input {
    File tsv_file
    Int n_chunks = 100
  }

  command <<<
    set -euxo pipefail

    # Count total lines including header
    total_lines=$(wc -l < ~{tsv_file})

    # Calculate lines per chunk (subtract 1 for header)
    lines_per_chunk=$(( (total_lines - 1 + ~{n_chunks} - 1) / ~{n_chunks} ))

    # Extract header
    head -n 1 ~{tsv_file} > header.tsv
    awk -F'\t' 'NR==1 || $4=="snp"' ~{tsv_file} > filtered.tsv
    rm ~{tsv_file}

    # Split remaining lines
    tail -n +2 filtered.tsv | split -l $lines_per_chunk - chunk_

    rm filtered.tsv

    # Prepend header to each chunk
    for f in chunk_*; do
        cat header.tsv "$f" > "${f}.tsv"
    done

    # Move chunk files to final output names
    mkdir -p chunks
    mv chunk_*.tsv chunks/
  >>>

  output {
    Array[File] chunk_files = glob("chunks/*.tsv")
  }

  runtime {
    docker: "ubuntu:latest"
    memory: "4G"
    disks: "local-disk 100 HDD"
  }
}

