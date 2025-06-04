# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0

workflow STEP_5_COUNT_VARIANTS {
  input {
    String step_4_output_dir       # Directory to STEP_4 Output VCFs

    Array[File] fam_files
    File rf_predictions
    File segmental_duplications_bed

    File subjects_list
    Int num_subjects = 4638 
    String storage_directory = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
  }
  call filter_fam_file {
    input:
      ufc_subject_list = subjects_list,
      fam_files = fam_files
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
    call make_annotation_file {
      input:
        rf_predictions = gather_tsvs.out[i] #rf_predictions
    }
  }

  Array[File] shards = sort_vcf_list.vcf_arr

  # Loop through the vcf shards and annotate each Variant with tp_prob
  scatter(i in range(length(shards)-negative_shards)){
    call annotate_vcf_with_tp_prob as annotate_vcf_with_tp_prob_snps {
      input:
        vcf = shards[i],
        tp_prob_map = make_annotation_file.snp_annotations[i]
    }
    call annotate_vcf_with_tp_prob as annotate_vcf_with_tp_prob_indels {
      input:
        vcf = shards[i],
        tp_prob_map = make_annotation_file.indel_annotations[i],
        variant_type = "indels"
    }
  }

  Array[File] snp_shards = select_all(annotate_vcf_with_tp_prob_snps.vcf_arr)
  Array[File] indel_shards = select_all(annotate_vcf_with_tp_prob_indels.vcf_arr)

  # Count Variants per Person
  scatter (i in range(length(snp_shards)-0)){
    call count_variants as count_snps {
      input:
        vcf = snp_shards[i],
        subject_list = subjects_list
    }
    call count_denovo_variants as count_denovo_snps {
      input:
        vcf = snp_shards[i],
        fam_file = filter_fam_file.output_fam_file,
        segmental_duplications_bed = segmental_duplications_bed
    }
  }

  # Count Variants per Person
  scatter (i in range(length(indel_shards))){
    call count_variants as count_indels {
      input:
        vcf = indel_shards[i],
        subject_list = subjects_list
    }
    call count_denovo_variants as count_denovo_indels {
      input:
        vcf = indel_shards[i],
        fam_file = filter_fam_file.output_fam_file,
        segmental_duplications_bed = segmental_duplications_bed
    }
  }

  # Count summary statistics for variants of the whole cohort
  #call aggregate_counts as aggregate_snps_count {
  #  input:
  #    sharded_counts = count_snps.shard_variant_count,
  #    variant_type = "snp"
  #}
  call aggregate_counts as aggregate_denovo_snps_count {
    input:
      sharded_counts = count_denovo_snps.shard_variant_count,
      variant_type = "denovo_snp"
  }

  #call aggregate_counts as aggregate_indels_count {
  #  input:
  #    sharded_counts = count_indels.shard_variant_count,
  #    variant_type = "indel"
  #}
  call aggregate_counts as aggregate_denovo_indels_count {
    input:
      sharded_counts = count_denovo_indels.shard_variant_count,
      variant_type = "denovo_indel"
  }

  # Graph the Stats of Interest
  #call GraphThresholdStats as GraphTotalVariants{
  #  input:
  #    snps_count = aggregate_snps_count.out,
  #    indels_count = aggregate_indels_count.out,
  #    file_output = "total_variants.png"
  #}

  call GraphThresholdStats as GraphTotalDenovos{
    input:
      snps_count = aggregate_denovo_snps_count.out,
      indels_count = aggregate_denovo_indels_count.out,
      file_output = "denovo_variants.png"
  }

  #call copy_to_storage {
  #  input:
  #    total_variants = GraphTotalVariants.graph,
  #    denovo_variants = GraphTotalVariants.graph, ## Edit this once DeNovos are Done!
  #    variants_tsv = aggregate_snps_count.out,
  #    denovo_variants_tsv = aggregate_snps_count.out,
  #    storage_directory = storage_directory
  #}
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
    Array[File] out = read_lines("rf_predictions.sharded.list")
  }
}

task aggregate_counts {
  input {
    Array[File] sharded_counts
    String variant_type
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  import numpy as np

  # Predefine TP_PROB values from 0 to 1 in increments of 0.01
  tp_prob_values = np.arange(0, 1.01, 0.01)

  # Initialize an aggregated DataFrame with TP_PROB as index and placeholders for subjects
  aggregated_df = pd.DataFrame({"TP_PROB": tp_prob_values})

  # List of input files
  files_list = "~{sep=' ' sharded_counts}".split()

  # Initialize an empty DataFrame for final aggregation
  final_df = None


  # Initialize variables to track the highest counts and corresponding file
  highest_tp0_count = 0
  highest_tp02_count = 0
  best_tp0_file = None
  best_tp02_file = None

  # Process each file and aggregate subject columns
  for file_path in files_list:
    # Read the file into a DataFrame
    df = pd.read_csv(file_path, sep="\t")

    # Ensure TP_PROB values are rounded to 2 decimal places and sorted
    df["TP_PROB"] = df["TP_PROB"].apply(lambda x: round(x, 2))
    df = df.sort_values("TP_PROB").set_index("TP_PROB")


    # Get sum of values in the row at TP_PROB[0] and TP_PROB[0.2]
    tp0_sum = df.loc[0].sum() if 0 in df.index else 0
    tp02_sum = df.loc[0.2].sum() if 0.2 in df.index else 0

    # Check if this file has the highest count at TP_PROB[0]
    if tp0_sum > highest_tp0_count:
        highest_tp0_count = tp0_sum
        best_tp0_file = file_path

    # Check if this file has the highest count at TP_PROB[0.2]
    if tp02_sum > highest_tp02_count:
        highest_tp02_count = tp02_sum
        best_tp02_file = file_path

    # Initialize final_df if it's the first file
    if final_df is None:
      final_df = df
    else:
      # Add the data from this file to the final aggregation
      final_df = final_df.add(df, fill_value=0)

  # Report the file with the highest counts for TP_PROB[0] and TP_PROB[0.2]
  print(f"File with highest TP_PROB[0] count: {best_tp0_file} ({highest_tp0_count})")
  print(f"File with highest TP_PROB[0.2] count: {best_tp02_file} ({highest_tp02_count})")

  # Calculate statistics for each TP_PROB across all subjects
  stats_df = pd.DataFrame(index=final_df.index)
  stats_df["MIN"] = final_df.min(axis=1).astype(int)
  stats_df["Q1"] = final_df.quantile(0.25, axis=1).astype(int)
  stats_df["Q2"] = final_df.median(axis=1).astype(int)
  stats_df["Q3"] = final_df.quantile(0.75, axis=1).astype(int)
  stats_df["MAX"] = final_df.max(axis=1).astype(int)

  # Reset the index for clean output
  stats_df.reset_index(inplace=True)
  stats_df.rename(columns={"index": "TP_PROB"}, inplace=True)

  # Save the aggregated DataFrame to a file
  output_file = "~{variant_type}.counts.tsv"
  stats_df.to_csv(output_file, sep="\t", index=False)
  final_df.to_csv("~{variant_type}.final_df.tsv",sep='\t',index=False)
  print(f"Aggregated data saved to {output_file}")

  CODE
  >>>
  output {
    File out = "~{variant_type}.counts.tsv"
    File out2 = "~{variant_type}.final_df.tsv"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "4GB"
    disks: "local-disk 16 HDD"
    preemptible: 3
  }
}

task make_annotation_file {
  input {
    File rf_predictions
  }
  Int default_mem_gb = ceil(size(rf_predictions,"GB") * 3) + 4
  command <<<
  set -euxo pipefail
  zcat ~{rf_predictions} > rf_predictions.tsv
  awk -F'\t' '!($2 ~ /\*/) {print}' rf_predictions.tsv > rf_predictions.nostar.tsv
  if [ $(wc -l < rf_predictions.nostar.tsv) -eq 0 ]; then
    touch snp_annotation_file.tsv
    touch indel_annotation_file.tsv
    exit 0
  fi
  python3 <<CODE
  import pandas as pd
  import re

  print("Starting Python3 Code")

  # Function to extract fp_prob and tp_prob from the string
  def extract_probs(prob_str):
        # Use regular expressions to extract the numeric values for 'False' and 'True'
        fp_match = re.search(r'"False","value":([0-9.]+)', prob_str)
        tp_match = re.search(r'"True","value":([0-9.]+)', prob_str)

        # Extract the values, default to None if not found
        fp_prob = float(fp_match.group(1)) if fp_match else None
        tp_prob = float(tp_match.group(1)) if tp_match else None

        return fp_prob, tp_prob


  # Load the file into a DataFrame
  df = pd.read_csv("rf_predictions.nostar.tsv", sep="\t")
  df['fp_prob'], df['tp_prob'] = zip(*df['rf_probability'].apply(extract_probs))

  print("Loaded Initial Data Frame")

  # Extract CHROM, POS, REF, and ALT from 'locus' and 'alleles'
  df[['CHROM', 'POS']] = df['locus'].str.split(':', expand=True)
  df[['REF', 'ALT']] = pd.DataFrame(df['alleles'].str.strip('[]"').str.split('","').tolist(), index=df.index)

  print("2")
  print(df.head())
  df_snp = df[df['allele_type'] == "snp"]
  df_indel = df[df['allele_type'] != "snp"]
  
  print("3")
  # Reorder and save to file
  df_snp_out = df_snp[['CHROM', 'POS', 'REF', 'ALT', 'tp_prob']]
  df_indel_out = df_indel[['CHROM', 'POS', 'REF', 'ALT', 'tp_prob']]
  df_snp_out.to_csv("snp_annotation_file.tsv", sep="\t", index=False, header=["CHROM", "POS", "REF", "ALT", "TP_PROB"])
  df_indel_out.to_csv("indel_annotation_file.tsv", sep="\t", index=False, header=["CHROM", "POS", "REF", "ALT", "TP_PROB"])

  CODE
  >>>
  runtime {
    docker:"vanallenlab/pydata_stack"
    preemptible: 3
    #memory: "120GB"
    #disks: "local-disk 250 HDD"
  }
  output {
    File snp_annotations = "snp_annotation_file.tsv"
    File indel_annotations = "indel_annotation_file.tsv"
  }
}

task annotate_vcf_with_tp_prob {
  input {
    File vcf
    File tp_prob_map
    String variant_type = "snp"
  }
  command <<<
  set -euxo pipefail
 
  if [ $(wc -l < ~{tp_prob_map} ) -eq 0 ]; then
    #mv ~{vcf} annotated.vcf.gz
    exit 0
  fi
 
  bgzip -c ~{tp_prob_map} > annotation.tsv.gz
  tabix -s 1 -b 2 -e 2 -S 1 annotation.tsv.gz
  
  # Step 2: Filter the VCF by VARIANT_TYPE and process SNPs or Indels
  if [[ "~{variant_type}" == "snp" ]]; then
    echo "Filtering for SNPs..."
    bcftools view -i 'INFO/VARIANT_TYPE == "snp"' ~{vcf} -O z -o filtered.vcf.gz
  else
    echo "Filtering for Indels..."
    bcftools view -i 'INFO/VARIANT_TYPE != "snp"' ~{vcf} -O z -o filtered.vcf.gz
  fi

  bcftools annotate \
    -a annotation.tsv.gz \
    -c CHROM,POS,REF,ALT,TP_PROB \
    -h <(echo "##INFO=<ID=TP_PROB,Number=1,Type=Float,Description=\"Probability that this variant allele is a true variant as determined by a Random Forest model\">") \
    -o annotated.vcf.gz \
    -Oz filtered.vcf.gz
  >>>
  runtime{
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
    memory: "4GB"
    disks: "local-disk 16 HDD"
  }
  output {
    File? vcf_arr = "annotated.vcf.gz"
  }
}

task count_denovo_variants {
  input {
    File vcf
    File fam_file
    File segmental_duplications_bed
  }
  #Int default_mem_gb = ceil(size(vcf,"GB") * 3) + 4

  command <<<
  set -euxo pipefail
  tail -n +2 < ~{fam_file} | cut -f2-4 | tr '\t' '\n' | sort -u > family_people.list
  bcftools view -S family_people.list --force-samples ~{vcf} -Oz -o filtered1.vcf.gz
  bcftools view -e 'ALT="*"' -Oz -o nostar.vcf.gz filtered1.vcf.gz
  bcftools view -T ^~{segmental_duplications_bed} nostar.vcf.gz -Oz -o nostar.segdup.vcf.gz 
  bcftools query -i 'AF<0.001 & GT="0/0" & FORMAT/DP>10 & FORMAT/GQ>20' \
     -f '%CHROM:%POS:%REF:%ALT\t[%SAMPLE,]\n' nostar.segdup.vcf.gz > parents.txt \
      || { echo "Error filtering variants"; exit 1; }

  bcftools query -i 'AF<0.001 & GT="alt" & FORMAT/DP>10 & FORMAT/GQ>20' \
     -f '%CHROM:%POS:%REF:%ALT\t%INFO/TP_PROB\t[%SAMPLE,]\n' nostar.segdup.vcf.gz > children.txt

  # Clear Disk
  rm ~{vcf} nostar.vcf.gz nostar.segdup.vcf.gz


  # Initialize output file
  output_file="subject_variant_counts.txt"
  echo -e "Subject\tNum_Variants" > $output_file

  python3 <<CODE
  print("Things are working")
  import pandas as pd
  import numpy as np

  # File paths
  variant_samples_file = "children.txt"
  variant_samples_and_missing_gt = "parents.txt"
  fam_file = "~{fam_file}"
  output_file = "tp_prob_variant_denovo_counts.tsv"

  family_dict = {}

  print("Checkpoint 1")

  # Read the subject list
  with open(fam_file, "r") as f:
    # Skip the header
    next(f)

    for line in f:
      columns = line.strip().split()

      # Extract IID, PID, and MID (indexing based on your .fam format)
      FID, IID, PID, MID, SEX, PHENOTYPE = columns

      # Add to dictionary: IID as the key, and (PID, MID) as the value
      family_dict[IID] = (PID, MID)

  # Initialize an empty DataFrame with TP_PROB as the index
  subjects = list(family_dict.keys())
  tp_prob_values = [round(x, 2) for x in list(np.arange(0, 1.01, 0.01))]
  tp_prob_df = pd.DataFrame(index=tp_prob_values, columns=subjects).fillna(0)

  # Initialize an empty dictionary
  alt_and_missing_gt_dict = {}

  print("Checkpoint 2")

  # Open the TSV file and process each line
  
  with open(variant_samples_and_missing_gt, "r") as f:
      for line in f:
          # Split the line by tab characters
          parts = line.strip().split("\t")
        
          chrom_pos = parts[0]  # First column is the key
          samples = parts[1]  # Third column contains the sample list
        
          # Convert sample list to a set (ignoring empty entries)
          sample_set = set(samples.split(",")) if samples else set()
        
          # Store in the dictionary
          alt_and_missing_gt_dict[chrom_pos] = sample_set


  print("Checkpoint 3")
  # Process the variant_samples file
  with open(variant_samples_file, "r") as f:
      for line in f:
          print(line)
          # Split the line into components
          chrom_pos, tp_prob, samples = line.strip().split("\t")
          if tp_prob == ".":
            continue
          tp_prob = round(float(tp_prob), 2)  # Ensure TP_PROB is rounded
          sample_list = [s for s in samples.split(",") if s]   # Extract individual samples
          potential_sample_list = alt_and_missing_gt_dict.get(chrom_pos, set())

          # Increment the count for each sample at the given TP_PROB
          for sample in sample_list:
              # Make sure that the sample is in a trio in our dataset
              if sample not in family_dict.keys():
                continue
              # Make sure that neither of the parents also have that variant
              father, mother = family_dict[sample]
              if father in potential_sample_list and mother in potential_sample_list:
                tp_prob_df.loc[tp_prob_df.index <= tp_prob, sample] += 1
                print(f"Found trio info for a sample: {sample}")
              else:
                print(f"Could not find trio info for a sample: {sample}")

  print("Checkpoint 4")
  ## Reset index and save the output file
  tp_prob_df.index.name = "TP_PROB"
  tp_prob_df.reset_index(inplace=True)
  tp_prob_df.to_csv(output_file, sep="\t", index=False)

  print(f"Variant counts by TP_PROB saved to {output_file}")
  CODE
  >>>
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
    memory: "6G"
    #disks: "local-disk 32 HDD"
  }
  output {
    File shard_variant_count = "tp_prob_variant_denovo_counts.tsv"
  }
}

task count_variants {
  input {
    File vcf
    File subject_list
  }
  Int default_mem_gb = ceil(size(vcf,"GB") * 1.5) + 2

  command <<<
  set -euxo pipefail
  bcftools view -e 'ALT="*"' -Oz -o nostar.vcf.gz ~{vcf}
  bcftools query -i 'GT="alt"' \
      -f '%CHROM:%POS\t%INFO/TP_PROB\t[%SAMPLE,]\n' nostar.vcf.gz > variant_samples.txt \
      || { echo "Error filtering variants"; exit 1; }

  # Clear Disk
  rm ~{vcf}

  echo "There are this many lines in the vcf"
  wc -l variant_samples.txt

  echo "This is the head of variant_samples.txt"
  head variant_samples.txt

  # Initialize output file
  output_file="subject_variant_counts.txt"
  echo -e "Subject\tNum_Variants" > $output_file
  
  python3 <<CODE
  import pandas as pd
  import numpy as np

  # File paths
  variant_samples_file = "variant_samples.txt"
  subject_list_file = "~{subject_list}"
  output_file = "tp_prob_variant_counts.tsv"

  # Read the subject list
  with open(subject_list_file, "r") as f:
      subjects = [line.strip() for line in f.readlines()]

  print("Subjects in tp_prob_df:", len(subjects))

  # Initialize an empty DataFrame with TP_PROB as the index
  tp_prob_values = [round(x, 2) for x in list(np.arange(0, 1.01, 0.01))]
  tp_prob_df = pd.DataFrame(index=tp_prob_values, columns=subjects).fillna(0)

  # Process the variant_samples file
  with open(variant_samples_file, "r") as f:
      for line in f:
          # Split the line into components
          chrom_pos, tp_prob, samples = line.strip().split("\t")
          if tp_prob == ".":
            continue
          tp_prob = round(float(tp_prob), 2)  # Ensure TP_PROB is rounded
          sample_list = [s for s in samples.split(",") if s]   # Extract individual samples
        
          print(f"Processing TP_PROB: {tp_prob}, Sample Length: {sample_list}")
          # Increment the count for each sample at the given TP_PROB
          for sample in sample_list:
              if sample in tp_prob_df.columns:
                  tp_prob_df.loc[tp_prob_df.index <= tp_prob, sample] += 1
              else:
                print(f"Sample not found in subjects: {sample}")

  # Reset index and save the output file
  tp_prob_df.index.name = "TP_PROB"
  tp_prob_df.reset_index(inplace=True)
  tp_prob_df.to_csv(output_file, sep="\t", index=False)

  print(f"Variant counts by TP_PROB saved to {output_file}")
  CODE
  >>>
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
    #memory: "16GB"
    disks: "local-disk 32 HDD"
  }
  output {
    File shard_variant_count = "tp_prob_variant_counts.tsv"
  }
}

task copy_to_storage {
  input {
    File total_variants
    File denovo_variants
    File variants_tsv
    File denovo_variants_tsv

    String storage_directory
  }
  command <<<
  set -eu -o pipefail
  gsutil -m cp ~{total_variants} gs://~{storage_directory}/STEP_5_COUNT_VARIANTS/output/graphs
  gsutil -m cp ~{denovo_variants} gs://~{storage_directory}/STEP_5_COUNT_VARIANTS/output/graphs

  gsutil -m cp ~{variants_tsv} gs://~{storage_directory}/STEP_5_COUNT_VARIANTS/output/tsvs
  gsutil -m cp ~{denovo_variants_tsv} gs://~{storage_directory}/STEP_5_COUNT_VARIANTS/output/tsvs
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
    disks: "local-disk 32 HDD"
  }
}

task gather_vcfs {
  input {
    String dir
  }
  command <<<
  set -euxo pipefail
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

task Probabilities {
    input {
      File denovo
      File rf_predictions
    }

  command <<<
  head -n 1 ~{rf_predictions} > denovo_rf_predictions.tsv
  grep -Fwf ~{denovo} ~{rf_predictions} >> denovo_rf_predictions.tsv

  head -n 1 ~{rf_predictions} > denovo_snp_predictions.tsv
  grep 'snp' denovo_rf_predictions.tsv >> denovo_snp_predictions.tsv

  head -n 1 ~{rf_predictions} > denovo_indel_predictions.tsv
  grep -v 'snp' denovo_rf_predictions.tsv >> denovo_indel_predictions.tsv
  python3 <<CODE
  import pandas as pd
  import re

  # Load the two TSV files
  df_snp = pd.read_csv("denovo_snp_predictions.tsv", sep='\t')
  df_indel = pd.read_csv("denovo_indel_predictions.tsv", sep='\t')

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
  df_snp['fp_prob'], df_snp['tp_prob'] = zip(*df_snp['rf_probability'].apply(extract_probs))
  df_indel['fp_prob'], df_indel['tp_prob'] = zip(*df_indel['rf_probability'].apply(extract_probs))

  # Ensure that 'tp_prob' is a float (it should already be a float, but just to be safe)
  df_snp['tp_prob'] = df_snp['tp_prob'].astype(float)
  df_indel['tp_prob'] = df_indel['tp_prob'].astype(float)

  # Sort the DataFrame by 'tp_prob' in ascending order
  df_snp = df_snp.sort_values(by='tp_prob', ascending=True)
  df_indel = df_indel.sort_values(by='tp_prob', ascending=True)

  # Extract the 'tp_prob' column and convert it to a comma-separated string
  tp_prob_string_snp = ','.join(df_snp['tp_prob'].astype(str))
  tp_prob_string_indel = ','.join(df_indel['tp_prob'].astype(str))

  # Save the result to a file
  with open('tp_prob_output_snp.txt', 'w') as f:
    f.write(tp_prob_string_snp + '\n')
  with open('tp_prob_output_indel.txt', 'w') as f:
    f.write(tp_prob_string_indel + '\n')
  CODE
  >>>

  output {
      File tp_prob_output_snp = "tp_prob_output_snp.txt"
      File tp_prob_output_indel = "tp_prob_output_indel.txt"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    disks: "local-disk 20 HDD"
  }
}

task GraphThresholdStats {
    input {
      File snps_count
      File indels_count
      String file_output
    }

    command <<<    
    set -eu -o pipefail
    python3 <<CODE
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    # Define thresholds from 0 to 1 with 0.01 increments
    thresholds = np.arange(0, 1.01, 0.01)

    snps_df = pd.read_csv("~{snps_count}", sep='\t')
    snps_df["TP_PROB"] = snps_df["TP_PROB"].astype(float)
    snps_df.iloc[:, 1:] = snps_df.iloc[:, 1:].astype(int)
    indels_df = pd.read_csv("~{indels_count}", sep='\t')
    indels_df["TP_PROB"] = indels_df["TP_PROB"].astype(float)
    indels_df.iloc[:, 1:] = indels_df.iloc[:, 1:].astype(int)

    # Plot both SNPs and Indels
    plt.plot(snps_df["TP_PROB"], snps_df["Q2"], label='Median (SNPs)', color='b')
    plt.fill_between(snps_df["TP_PROB"], snps_df["Q1"], snps_df["Q3"], color='b', alpha=0.2, label='IQR (SNPs)')
    plt.fill_between(snps_df["TP_PROB"], snps_df["MIN"], snps_df["MAX"], color='b', alpha=0.03, label='Range (SNPs)')

    plt.plot(indels_df["TP_PROB"], indels_df["Q2"], label='Median (Indels)', color='r')
    plt.fill_between(indels_df["TP_PROB"], indels_df["Q1"], indels_df["Q3"], color='r', alpha=0.2, label='IQR (Indels)')
    plt.fill_between(indels_df["TP_PROB"], indels_df["MIN"], indels_df["MAX"], color='r', alpha=0.03, label='Range (Indels)')

    # Add labels and title    
    plt.xlabel('Threshold', fontsize=14)
    plt.ylabel('Variant Count per Patient', fontsize=14)
    plt.title('True Positive Variants per Threshold', fontsize=16)
    plt.legend(loc='best', fontsize=12)

    # Increase tick label font size
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.tight_layout()
    # Save the plot
    plt.savefig('~{file_output}',dpi=300)

    CODE
    >>>

    output {
        File graph = "~{file_output}"
    }

    runtime {
        docker: "vanallenlab/pydata_stack"
        memory: "8G"
        preemptible: 3
    }
}

# The goal of this code is to finalize the fam file for downstream use
task filter_fam_file {
  input {
    File ufc_subject_list
    Array[File] fam_files
  }

  command <<<
    # Get list of sample ids in vcf file
    cp ~{ufc_subject_list} sample_ids.txt

    # Concatenate all Fam Files
    for fam_file in ~{sep=" " fam_files}; do
      tail -n +2 $fam_file >> all_families.fam
    done

    # Read sample IDs into a set
    sample_ids=$(awk '{print $1}' sample_ids.txt)

    # Filter the fam file based on subjid, mother, and father
    awk 'NR==FNR {ids[$1]; next} ($2 in ids && $3 in ids && $4 in ids)' sample_ids.txt all_families.fam > trios.tmp.fam

    # Count the number of lines in filtered.fam and write to num_trios.txt
    wc -l < trios.tmp.fam > num_trios.txt

    awk 'BEGIN {print "FID IID PID MID SEX PHENOTYPE"} {print $0, 0}' trios.tmp.fam > trios.fam 
  >>>

  output {
    File output_fam_file = "trios.fam"
    Int num_trios = read_int("num_trios.txt")
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
  }
}

