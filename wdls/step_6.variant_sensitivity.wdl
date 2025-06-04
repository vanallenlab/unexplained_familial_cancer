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
  
  # 1000Genome Data
  call concatenateFiles as concatenateFiles_snps{
    input:
      files = Extract_SNPs_From_VCF.variants,
      callset_name = "1000G_snps"
  }
  call Probabilities as SNP_Probabilities{
    input:
      variants = concatenateFiles_snps.out,
      rf_predictions = rf_predictions,
      snp_or_indel = "snp"
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

  # Doubleton, Tripleton, False Positives
  #call Probabilities as Doubletons_Probabilities{
  #  input:
  #    variants = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/Extract_AF/01a73628-a1d3-4efc-a3cf-2ce29872ccdb/call-concatenateDoubletonTPs/cacheCopy/doubleton_tps_list.tsv",
  #    rf_predictions = rf_predictions,
  #    snp_or_indel = "All"
  #}
  #call Probabilities as Tripletons_Probabilities{
  #  input:
  #    variants = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/Extract_AF/2816eb47-faa1-46af-9631-bd4dbb2d1a16/call-concatenateTripletons/tripletons_list.tsv",
  #    rf_predictions = rf_predictions,
  #    snp_or_indel = "All"
  #}
  #call Probabilities as False_Positives_Probabilities{
  #  input:
  #    variants = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/Extract_AF/01a73628-a1d3-4efc-a3cf-2ce29872ccdb/call-concatenateFPs/cacheCopy/false_positive_list.tsv",
  #    rf_predictions = rf_predictions,
  #    snp_or_indel = "indel"
  #}
  # Plotting Sensitivity
  call plot_sensitivity as plot_sensitivity_snps{
    input:
      tsv_file = SNP_Probabilities.tp_prob_output,
      snp_or_indel = "snp",
      max_sensitivity = SNP_Probabilities.max_sensitivity
  }
  call plot_sensitivity as plot_sensitivity_indels{
    input:
      tsv_file = INDEL_Probabilities.tp_prob_output,
      snp_or_indel = "indel",
      max_sensitivity = INDEL_Probabilities.max_sensitivity
  }
  #call plot_sensitivity_all {
  #  input:
  #    doubletons = Tripletons_Probabilities.tp_prob_output,
  #    tripletons = Tripletons_Probabilities.tp_prob_output,
  #    #false_positives = Tripletons_Probabilities.tp_prob_output,
  #    snp_or_indel = "all",
  #    doubletons_max_sensitivity = Tripletons_Probabilities.max_sensitivity,
  #    tripletons_max_sensitivity = Tripletons_Probabilities.max_sensitivity
      #false_positives_max_sensitivity = Tripletons_Probabilities.max_sensitivity
  #}
}

task plot_sensitivity_all {
  input {
    File doubletons
    File tripletons
    #File false_positives
    String snp_or_indel
    Float doubletons_max_sensitivity
    Float tripletons_max_sensitivity
    #Float false_positives_max_sensitivity
  }

  command <<<
  # Exit immediately if a command exits with a non-zero status
  set -eu -o pipefail
  cp ~{doubletons} doubletons.tsv
  cp ~{tripletons} tripletons.tsv
  #cp {false_positives} false_positives.tsv

  python3 <<CODE
  import pandas as pd
  import matplotlib.pyplot as plt
  import numpy as np

  # Load the TSV files
  df_doubletons = pd.read_csv("doubletons.tsv", sep='\t')
  df_tripletons = pd.read_csv("tripletons.tsv", sep='\t')
  #df_false_positives = pd.read_csv("false_positives.tsv", sep='\t')

  # Sort by tp_prob in ascending order
  df_doubletons = df_doubletons.sort_values(by='tp_prob', ascending=True)
  df_tripletons = df_tripletons.sort_values(by='tp_prob', ascending=True)
  #df_false_positives = df_false_positives.sort_values(by='tp_prob', ascending=True)
  
  # Define probability thresholds from 0 to 1 in increments of 0.01
  thresholds = np.arange(0, 1.01, 0.01)

  # Initialize list to store sensitivity at each threshold
  doubletons_sensitivity = []
  tripletons_sensitivity = []
  #false_positives_sensitivity = []

  # Calculate sensitivity at each threshold
  doubleton_total_positives = len(df_doubletons)
  tripleton_total_positives = len(df_tripletons)
  #false_total_positives = len(df_false_positives)

  for threshold in thresholds:
      # Count how many rows have tp_prob >= current threshold
      doubletons_true_positives_at_threshold = df_doubletons[df_doubletons['tp_prob'] >= threshold].shape[0]
      tripletons_true_positives_at_threshold = df_tripletons[df_tripletons['tp_prob'] >= threshold].shape[0]
      #false_positives_at_threshold = df_false_positives[df_false_positives['tp_prob'] >= threshold].shape[0]

      # Sensitivity is TP / Total positives
      doubletons_sensitivity_at_threshold = (doubletons_true_positives_at_threshold / doubleton_total_positives) * ~{doubletons_max_sensitivity}
      tripletons_sensitivity_at_threshold = (tripletons_true_positives_at_threshold / tripleton_total_positives) * ~{tripletons_max_sensitivity}
      #false_positives_sensitivity_at_threshold = (false_positives_at_threshold / false_total_positives) * {false_positives_max_sensitivity}


      doubletons_sensitivity.append(doubletons_sensitivity_at_threshold)
      tripletons_sensitivity.append(tripletons_sensitivity_at_threshold)
      #false_positives_sensitivity.append(false_positives_sensitivity_at_threshold)

  # Plot the sensitivity curve
  plt.figure(figsize=(8, 6))
  plt.plot(thresholds, doubletons_sensitivity, linestyle='-', color='b',label='Doubletons')
  plt.plot(thresholds, tripletons_sensitivity, linestyle='-', color='r', label = 'Tripletons')
  #plt.plot(thresholds, false_positives_sensitivity, linestyle='-', color='g', label = 'False Positives')

  plt.xlabel('tp_prob')
  plt.ylabel('Sensitivity')  
  plt.title('Sensitivity Curve for Doubleton and Tripleton Variants')
  plt.grid(True)
  plt.legend()

  # Save the plot
  plt.savefig("sensitivity_plot_all_metrics_~{snp_or_indel}.png")
  plt.close()

  # Create a DataFrame with thresholds and sensitivity
  df_sensitivity = pd.DataFrame({
      'tp_prob': thresholds,
      'doubletons_sensitivity': doubletons_sensitivity,
      'tripletons_sensitivity': tripletons_sensitivity,
      #'false_positives_sensitivity': false_positives_sensitivity,
  })

  # Save the DataFrame to a TSV file
  df_sensitivity.to_csv("sensitivity_output_all_metrics_~{snp_or_indel}.tsv", sep='\t', index=False)
  CODE
  >>>

  output {
    File sensitivity_plot = "sensitivity_plot_all_metrics_~{snp_or_indel}.png"
    File sensitivity_tsv = "sensitivity_output_all_metrics_~{snp_or_indel}.tsv"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    #memory: "4G"
    #disks: "local-disk 10 HDD"
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
    memory: "8G"
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
  grep -Fwf ~{variants} tmp.tsv  >> filtered_predictions.tsv
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
    memory: "16G"
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
