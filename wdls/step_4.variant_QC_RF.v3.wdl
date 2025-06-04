# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0

workflow STEP_4_RANDOM_FOREST {
  input {
    String hard_filter_output_dir	# Directory to STEP_3 Output VCFs

    # VCF files for HQ sites as designated by 1000G. Used for HWE analysis
    Array[File] HQ_sites_vcfs
    Array[File] HQ_sites_vcf_idxs
    
    Array[File] fam_files
    File subjects_list

    String study_cohort_callset_name = "ufc_whole_genome"
    String reference_database_callset_name = "gnomad_whole_genome"

    String storage_directory = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
  }

  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call gather_vcfs {
    input:
      dir = hard_filter_output_dir
  }

  call sort_vcf_list {
    input:
      unsorted_vcf_list = gather_vcfs.vcf_list
  }

  Array[Pair[File, File]] shards = zip(sort_vcf_list.vcf_arr,sort_vcf_list.vcf_idx_arr)

  scatter (idx in range(length(shards) - 0)){
    # Annotate the VCF w/ features to load into random forest
    call annotate_vcf_with_gatk {
      input:
        vcf = shards[idx].left,
        vcf_idx = shards[idx].right,
        fam_files = fam_files
    }
    call annotate_vcf_with_gatk_part2 {
      input:
        vcf = annotate_vcf_with_gatk.output_vcf,
        vcf_idx = annotate_vcf_with_gatk.output_vcf_idx,
        subjects_list = subjects_list,
        fam_files = fam_files
    } 
    call annotate_vcf_with_pab_max_expr {
      input:
        vcf = annotate_vcf_with_gatk_part2.output_vcf
    }

    call copy_vcfs_to_storage {
      input:
        vcf = annotate_vcf_with_pab_max_expr.output_vcf,
        storage_directory = storage_directory,
        shard_num = idx
    }

    # Grab training variants
    call grab_rf_training_variants {
      input:
        vcf = annotate_vcf_with_pab_max_expr.output_vcf,
        TP_reference_vcfs = HQ_sites_vcfs,
        TP_reference_vcf_idxs = HQ_sites_vcf_idxs
    }

    call create_fam_tp_training_set {
      input:
        vcf = shards[idx].left,
        fam_files = fam_files
    }
  }

  # Create Full Lists of True Positives and False Positives based on HQ sites
  call concatenateFiles as concatenateTPs{
    input:
      files = grab_rf_training_variants.true_positive_list,
      callset_name = "true_positive_list"
  }
  call concatenateFiles as concatenateFPs{
    input:
      files = grab_rf_training_variants.false_positive_list,
      callset_name = "false_positive_list"
  }

  # Create Full Lists of True Positives based on Transmitted Doubletons
  call concatenateFiles as concatenateDoubletonTPs{
    input:
      files = create_fam_tp_training_set.tp_list,
      callset_name = "doubleton_tps_list"
  }

  # Create Full Lists of True Positives based on Transmitted Doubletons
  call concatenateFiles as concatenateDoubletonTPs_for_RF_results{
    input:
      files = create_fam_tp_training_set.doubletons_list,
      callset_name = "doubleton_tps_list_for_RF_results"
  }

  call create_training_sets_python {
    input:
      tp_list = concatenateTPs.out,
      tp_fam = concatenateDoubletonTPs.out,
      fp_list = concatenateFPs.out
  }

  # Start Imputing Values
  scatter (i in range(length(annotate_vcf_with_pab_max_expr.output_vcf) - 0)){
    call gather_mean_feature_values as gather_mean_pab_max_expr { 
      input:
        vcf = annotate_vcf_with_pab_max_expr.output_vcf[i],
        feature = "pab_max_expr"
    }
    call gather_mean_feature_values as gather_mean_AS_InbreedingCoeff {
      input:
        vcf = annotate_vcf_with_pab_max_expr.output_vcf[i],
        feature = "AS_InbreedingCoeff"
    }
    call gather_mean_feature_values as gather_mean_QD {
      input:
        vcf = annotate_vcf_with_pab_max_expr.output_vcf[i],
        feature = "QD"
    }
    call gather_mean_feature_values as gather_mean_MQRankSum {
      input:
        vcf = annotate_vcf_with_pab_max_expr.output_vcf[i],
        feature = "MQRankSum"
    }
    call gather_mean_feature_values as gather_mean_ReadPosRankSum {
      input:
        vcf = annotate_vcf_with_pab_max_expr.output_vcf[i],
        feature = "ReadPosRankSum"
    }
    call gather_mean_feature_values as gather_mean_SOR {
      input:
        vcf = annotate_vcf_with_pab_max_expr.output_vcf[i],
        feature = "SOR"
    }
  }

  # Calculate weighted means
  call calculate_weighted_mean as calculate_weighted_mean_pab_max_expr {
    input:
      weighted_means = gather_mean_pab_max_expr.weighted_mean,
      weights = gather_mean_pab_max_expr.weight,
  }

  call calculate_weighted_mean as calculate_weighted_mean_AS_InbreedingCoeff {
    input:
      weighted_means = gather_mean_AS_InbreedingCoeff.weighted_mean,
      weights = gather_mean_AS_InbreedingCoeff.weight,
  }

  call calculate_weighted_mean as calculate_weighted_mean_QD {
    input:
      weighted_means = gather_mean_QD.weighted_mean,
      weights = gather_mean_QD.weight,
  }

  call calculate_weighted_mean as calculate_weighted_mean_MQRankSum {
    input:
      weighted_means = gather_mean_MQRankSum.weighted_mean,
      weights = gather_mean_MQRankSum.weight,
  }

  call calculate_weighted_mean as calculate_weighted_mean_ReadPosRankSum {
    input:
      weighted_means = gather_mean_ReadPosRankSum.weighted_mean,
      weights = gather_mean_ReadPosRankSum.weight,
  }
  call calculate_weighted_mean as calculate_weighted_mean_SOR {
    input:
      weighted_means = gather_mean_SOR.weighted_mean,
      weights = gather_mean_SOR.weight,
  }

  # Create a Training Set by going through each Slice of VCF
  scatter (i in range(length(annotate_vcf_with_pab_max_expr.output_vcf)-0)){
    call prepare_training_vcfs {
      input:
        vcf = annotate_vcf_with_pab_max_expr.output_vcf[i],
        tp_list = create_training_sets_python.final_tp_list,
        fp_list = create_training_sets_python.final_fp_list
    }

    # Convert VCF to Hail Table; Annotate w/ n_alleles, mixed_site, spanning_deletion, allele_type
    call prepare_training_vcfs_part2 {
      input:
        vcf = prepare_training_vcfs.labeled_vcf,
        mean_pab_max_expr = calculate_weighted_mean_pab_max_expr.weighted_mean,
        mean_AS_InbreedingCoeff = calculate_weighted_mean_AS_InbreedingCoeff.weighted_mean,
        mean_QD = calculate_weighted_mean_QD.weighted_mean,
        mean_MQRankSum = calculate_weighted_mean_MQRankSum.weighted_mean,
        mean_ReadPosRankSum = calculate_weighted_mean_ReadPosRankSum.weighted_mean,
        mean_SOR = calculate_weighted_mean_SOR.weighted_mean
    }

  }

  #call concatenate_hail_train_tables_v2 {
  #  input:
  #    hail_tables = prepare_training_vcfs_part2.hail_train_table
  #}

  #call train_hail_RF {
  #  input:
  #    hail_table = concatenate_hail_train_tables_v2.concatenated_table
  #}

  # Apply the Random Forest Model 
  scatter (i in range(length(prepare_training_vcfs.labeled_vcf) - 0)){
    call apply_hail_RF {
      input:
        tsv = prepare_training_vcfs_part2.hail_table[i],
        rf_model = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/STEP_4_RANDOM_FOREST/8b65beea-3118-48b6-bb7e-6db46849e2b0/call-train_hail_RF/rf.model.tar.gz",
        shard_num = i
    }

    call copy_tsvs_to_storage {
      input:
        tsv = apply_hail_RF.predictions,
        storage_directory = storage_directory,
        shard_num = i 
    }
  }
 

  #call concatenateFiles as concatenateRFpredictions{
  #  input:
  #    files = apply_hail_RF.predictions,
  #    callset_name = "ufc_rf_predictions.whole_genome"
  #}

  #call copy_to_storage {
  #  input:
  #    rf_model = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/cromwell-execution/STEP_4_RANDOM_FOREST/8b65beea-3118-48b6-bb7e-6db46849e2b0/call-train_hail_RF/rf.model.tar.gz",
  #    rf_predictions = concatenateRFpredictions.out,
  #    training_tp_list = create_training_sets_python.final_tp_list,
  #    training_fp_list = create_training_sets_python.final_fp_list,
  #    chr20_validation_tp_list = create_training_sets_python.chr20_tp_list,
  #    chr20_validation_fp_list = create_training_sets_python.chr20_fp_list,
  #    doubletons_list = concatenateDoubletonTPs.out,
  #    storage_directory = storage_directory
  #}
  output {
  }
}


task copy_tsvs_to_storage {
  input {
    File tsv
    String storage_directory
    Int shard_num
  }
  command <<<
  set -euxo pipefail
  tsv_bgz="$(basename ~{tsv}).bgz"
  echo "Copying over tsv: ~{shard_num}"
  cut -f1,2,3,4,15 ~{tsv} | bgzip -c > "$tsv_bgz"
  gsutil -m cp "$tsv_bgz" "gs://~{storage_directory}/STEP_4_RANDOM_FOREST/sharded_tsvs/"
  >>>
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
    #disks: "local-disk 50 HDD"
  }
}

task copy_to_storage {
  input {
    File rf_model
    File rf_predictions
    File training_tp_list
    File training_fp_list
    File doubletons_list
    File chr20_validation_tp_list
    File chr20_validation_fp_list

    String storage_directory
  }
  command <<<
  set -eu -o pipefail
  gsutil -m cp ~{rf_model} gs://~{storage_directory}/STEP_4_RANDOM_FOREST/
  gsutil -m cp ~{rf_predictions} gs://~{storage_directory}/STEP_4_RANDOM_FOREST/
  gsutil -m cp ~{training_tp_list} gs://~{storage_directory}/STEP_4_RANDOM_FOREST/
  gsutil -m cp ~{training_fp_list} gs://~{storage_directory}/STEP_4_RANDOM_FOREST/
  gsutil -m cp ~{doubletons_list} gs://~{storage_directory}/STEP_4_RANDOM_FOREST/
  gsutil -m cp ~{chr20_validation_tp_list} gs://~{storage_directory}/STEP_4_RANDOM_FOREST/
  gsutil -m cp ~{chr20_validation_fp_list} gs://~{storage_directory}/STEP_4_RANDOM_FOREST/
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
    disks: "local-disk 100 HDD"
  }
}


task gather_vcfs {
  input {
    String dir
  }
  command <<<
  set -eu -o pipefail
  gsutil ls ~{dir}/*chr*.vcf.gz > vcf.list
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


task copy_vcfs_to_storage {
  input {
    File vcf
    Int shard_num
    String storage_directory
  }
  command <<<
  set -euxo pipefail
  echo "Copying over vcf: ~{shard_num}"
  gsutil -m cp ~{vcf} gs://~{storage_directory}/STEP_4_RANDOM_FOREST/sharded_vcfs/
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
    disks: "local-disk 50 HDD"
  }
}


# Task to Make sure we are incorporating doubletons into TP for Hail RF
task create_fam_tp_training_set {
  input {
    Array[File] fam_files
    File vcf
  }
  command <<<
  # Make one big fam file
  for fam_file in ~{sep=" " fam_files}; do
    tail -n +2 $fam_file >> all_families.fam
  done

  # Initial Filtering for familial analysis
  bcftools view -i 'AC=2' ~{vcf} -o ac_filtered.vcf

  python3 <<CODE
  import pandas as pd
  import pysam

  # Load FAM file into a DataFrame
  fam_df = pd.read_csv('all_families.fam', sep='\t', names=['FAMILY_ID', 'SUBJECT_ID', 'MOTHER', 'FATHER', 'SEX'])

  # Create a dictionary to map child to parent(s)
  parent_child_map = {}
  for _, row in fam_df.iterrows():
    if row['MOTHER'] != '0':  # If a mother exists
        parent_child_map[(row['MOTHER'], row['SUBJECT_ID'])] = True
    if row['FATHER'] != '0':  # If a father exists
        parent_child_map[(row['FATHER'], row['SUBJECT_ID'])] = True

  # Open the VCF file using pysam
  vcf_in = pysam.VariantFile('ac_filtered.vcf')

  # Create a file to write the parent-child variants
  with open('fam_tp.list', 'w') as tp_list, open('doubletons.list', 'w') as doubletons_list:
    # Iterate over each record in the VCF file
    for record in vcf_in:
        # Filter genotypes to only those that are variants
        samples_with_variants = [sample for sample in record.samples if record.samples[sample]['GT'] != (0, 0)]
        
        # Check if exactly two samples have variants at this locus
        if len(samples_with_variants) == 2:
            sample1, sample2 = samples_with_variants
            
            # Check if the two samples are in a parent-child relationship
            if (sample1, sample2) in parent_child_map or (sample2, sample1) in parent_child_map:
                # Format the variant information for tp.list
                chrom = record.chrom
                pos = record.pos
                alt = record.alts[0]
                ref = record.ref
                
                # Write the variant to tp.list
                tp_list.write(f'{chrom}_{pos}_{ref}_{alt}\n')

                # Write to fam_tp_detailed.list with format chrom:pos\t["ref","alt"]
                doubletons_list.write(f'{chrom}:{pos}\t["{ref}","{alt}"]\n')

  CODE
  >>>
  runtime{
    docker:"vanallenlab/g2c_pipeline:latest"
    preemptible: 3
  }
  output{
    File tp_list = "fam_tp.list"
    File doubletons_list = "doubletons.list"
  }
}

task gather_mean_feature_values {
  input {
    File vcf
    String feature
  }
  command <<<
  set -euxo pipefail

  # Exclude positions where ALT contains '*'
  bcftools view -e 'ALT=="*"' ~{vcf} -Oz -o no_star.vcf.gz

  bcftools query -f '%INFO/~{feature}\n' no_star.vcf.gz | awk '$0 != "."' > values.list
  echo "Checkpoint 1"

  python3 <<CODE
  # Read the list of numbers from a file
  with open("values.list", "r") as file:
      # Convert each line to a float or integer, stripping whitespace
      numbers = [float(line.strip()) for line in file if line.strip()]

  # Calculate the total
  total = sum(numbers)

  # Count the numbers
  count = len(numbers)

  # Write the total to a new file
  with open("total.txt", "w") as file:
      file.write(f"{total}")

  # Write the count to another file
  with open("count.txt", "w") as file:
      file.write(f"{count}")  
  CODE
   
  >>>
  output {
    Float weighted_mean = read_float("total.txt")
    Int weight = read_int("count.txt")
  }
  runtime {
    memory: "8 GB"
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}

task calculate_weighted_mean {
  input {
    Array[Float] weighted_means
    Array[Float] weights
  }
  command <<<
  set -eu -o pipefail
  python3 <<CODE

  num = ~{sep="+" weighted_means}
  denom = ~{sep="+" weights}
  weighted_mean = num / denom

  print("Checkpoint 1")
  with open("weighted_mean.txt", "w") as file:
    file.write(f"{weighted_mean}")
  CODE
  >>>
  output {
    Float weighted_mean = read_float("weighted_mean.txt")
  }
  runtime {
    memory: "8 GB"
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task apply_hail_RF {
  input {
    File rf_model
    File tsv
    Int shard_num
  }
  command <<<
  set -euxo pipefail

  touch rf_predictions.~{shard_num}.tsv
  tar -xzvf ~{rf_model}

  python3 <<CODE
  import sys
  sys.path.append('/opt')
  import hail as hl
  import apply_rf_model as rf

  hl.init()

  # Define the desired types for columns
  desired_types = {
    'allele_type': hl.tstr,
    'qd': hl.tfloat64,
    'MQRankSum': hl.tfloat64,
    'ReadPosRankSum': hl.tfloat64,
    'SOR': hl.tfloat64,
    'pab_max_expr': hl.tfloat64,
    'Inbreeding_Coeff': hl.tfloat64,
    'has_star': hl.tbool,
    'n_alt_alleles': hl.tint64,
    'variant_type': hl.tstr,
    'was_mixed': hl.tbool,
    'TRAIN_LABEL': hl.tstr  # Import directly as boolean
  }

  #Load the Model
  print("Importing RF Model")
  rf_model = rf.load_model("rf.model")
  print("Imported RF Model")

  # Import the VCF
  ht = hl.import_table("~{tsv}",force_bgz=True, types=desired_types)
  #mt = hl.import_vcf("{vcf}", reference_genome='GRCh38', array_elements_required=False, force_bgz=True)
  variant_count = ht.count()

  if variant_count == 0:
    print("No variants found in this vcf shard")
    exit()
  else:
    print(f"The vcf shard contains {variant_count} variants.")

  features = ["allele_type", "pab_max_expr","Inbreeding_Coeff", "qd","MQRankSum","ReadPosRankSum","SOR","variant_type","n_alt_alleles","has_star","was_mixed"]

  print("Applying RF Model")
  out_ht = rf.apply_rf_model(ht, rf_model, features)
  
  # Assuming 'ht' is your Hail Table
  out_ht.export('rf_predictions.~{shard_num}.tsv')

  CODE
  >>>
  runtime{
    cpu: "4"
    memory: "8G"
    docker:"vanallenlab/hail_rf"
    preemptible: 3
  }
  output {
    File predictions = "rf_predictions.~{shard_num}.tsv"
  }
}

task create_training_sets_python {
  input {
    File tp_list
    File tp_fam
    File fp_list
  }
  command <<<
  python3 <<CODE
  import random
  import os

  # Define file paths
  tp_list_path = "tp_list.txt"
  tp_fam_path = "tp_fam.txt"
  fp_list_path = "fp_list.txt"

  # Copy input files to ensure consistent filenames
  os.system(f"cp ~{tp_list} tp_list.txt")
  os.system(f"cp ~{fp_list} fp_list.txt")
  os.system(f"cp ~{tp_fam} tp_fam.txt")

  # Function to read file lines into a list
  def read_file(file_path):
    with open(file_path, 'r') as f:
        with_chr20 = []
        without_chr20 = []
        for line in f:
            line = line.strip()
            if 'chr20' in line:
                with_chr20.append(line)
            else:
                without_chr20.append(line)
        return without_chr20, with_chr20

  # Function to write list to a file
  def write_file(file_path, data):
    with open(file_path, 'w') as f:
      f.write("\n".join(data) + "\n")

  # Function to randomly truncate a list to a specified length
  def truncate_list(data, target_length):
    if len(data) > target_length:
      return random.sample(data, target_length)
    return data

  # Read the files
  tp_list, tp_list_chr20 = read_file("tp_list.txt")
  tp_fam, tp_fam_chr20 = read_file("tp_fam.txt")
  fp_list, fp_list_chr20 = read_file("fp_list.txt")

  # Make Chr20 validation list:
  tp_chr20 = list(set(tp_list_chr20 + tp_fam_chr20))
  fp_chr20 = list(set(fp_list_chr20) - set(tp_chr20))
  write_file("chr20_tp.list", tp_chr20)
  write_file("chr20_fp.list", fp_chr20)

  # Make Whole Genome List:
  tp_wgs = list(set(tp_list + tp_fam))
  fp_wgs = list(set(fp_list) - set(tp_wgs))
  tp_high_quality = list(set(tp_wgs) - set(tp_fam))

  # Count initial lines
  total_tp_count = len(tp_wgs)
  tp_high_quality_count = len(tp_high_quality)
  tp_fam_count = len(tp_fam)
  fp_count = len(fp_wgs)

  print(f"Initial number of total TP sites: {total_tp_count}")
  print(f"Initial number of TP HQ sites: {tp_high_quality_count}")
  print(f"Initial number of TP doubletons: {tp_fam_count}")
  print(f"Initial number of FP lines: {fp_count} after pruning TP sites")

  # Truncate lists to equalize sizes
  if total_tp_count > fp_count:
    print("Randomly truncating TP High Quality list to match FP count...")
    
    target_length = fp_count - tp_fam_count
    tp_list_truncated = truncate_list(tp_high_quality, target_length)
    tp_list_combined = tp_list_truncated + tp_fam
    write_file("fp.list", fp_wgs)
    write_file("tp.list", tp_list_combined)
  elif fp_count > total_tp_count:
    print("Randomly truncating FP list to match TP count...")
    fp_list_truncated = truncate_list(fp_wgs, total_tp_count)
    write_file("tp.list", tp_wgs)
    write_file("fp.list", fp_list_truncated)
  else:
    print("TP and FP lists are already the same length.")
    write_file("fp.list", fp_wgs)
    write_file("tp.list", tp_wgs)
  print("Processing complete. Final lists written to tp.list and fp.list.")
  CODE
  >>>   
  output {
     File final_tp_list = "tp.list"
     File final_fp_list = "fp.list"
     File chr20_tp_list = "chr20_tp.list"
     File chr20_fp_list = "chr20_fp.list"
   }
   runtime {
     docker: "vanallenlab/pydata_stack"
     memory: "8 GB"
     preemptible: 3
   }
}

task prepare_training_vcfs {
  input {
    File vcf
    File tp_list
    File fp_list
  }
  Int default_disk_gb = ceil(size(vcf,"GB")) + 25

  command <<<
    set -euxo pipefail
    # Step 1: Annotate VCF with CHROM_POS_ALT in the ID field
    echo "Annotating VCF with CHROM_POS_ALT in the ID field..."
    awk -F'_' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tTrue"}' ~{tp_list} > tp_pretabix.list
    awk -F'_' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\tFalse"}' ~{fp_list} > fp_pretabix.list
    cat tp_pretabix.list fp_pretabix.list | sort -k1,1V -k2,2n > tmp.tsv

    bgzip -c tmp.tsv > annotation.tsv.gz
    tabix -s 1 -b 2 -e 2 annotation.tsv.gz

  bcftools index -t ~{vcf}
  bcftools annotate \
    -a annotation.tsv.gz \
    -c CHROM,POS,REF,ALT,TRAIN_LABEL \
    -h <(echo "##INFO=<ID=TRAIN_LABEL,Number=1,Type=String,Description=\"The Label that the RF model will use for training.\">") \
    -o annotated.vcf.gz \
    -Oz ~{vcf}

  >>>
  output {
    File labeled_vcf = "annotated.vcf.gz"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline:latest"
    disks: "local-disk ~{default_disk_gb} HDD"
    memory: "8 GB"
    preemptible: 3
  }
}

task annotate_vcf_with_gatk {
  input {
    File vcf
    File vcf_idx
    File Reference_Fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File Reference_Fasta_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    File Reference_Fasta_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    Array[File] fam_files
    Int? disk_gb
  }
  String vcf_basename=basename("~{vcf}")
  Int default_disk_gb = ceil(4 * size(vcf,"GB")) + 12
  command <<<
  set -eu -o pipefail

    # Make one big fam file
    for fam_file in ~{sep=" " fam_files}; do
      tail -n +2 $fam_file >> all_families.fam
    done
    awk '{print $0 "\t" 0}' all_families.fam > trios.fam


    gatk VariantAnnotator \
      -R ~{Reference_Fasta} \
      -V ~{vcf} \
      -O ~{vcf_basename} \
      -A VariantType \
      -A IndelClassify

    # Remove the original VCF to save space
    rm ~{vcf}
    
  >>>
  runtime {
    docker: "broadinstitute/gatk:4.6.0.0"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
  output {
    File output_vcf = "~{vcf_basename}"
    File output_vcf_idx = "~{vcf_basename}.tbi"
  }
}

task annotate_vcf_with_gatk_part2 {
  input {
    File vcf
    File vcf_idx
    File subjects_list
    File Reference_Fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File Reference_Fasta_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    File Reference_Fasta_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    Array[File] fam_files
    Int? disk_gb
  }
  String vcf_basename=basename("~{vcf}")
  Int default_disk_gb = ceil(4 * size(vcf,"GB")) + 12
  command <<<
    set -eu -o pipefail

    # Make one big fam file
    for fam_file in ~{sep=" " fam_files}; do
      tail -n +2 $fam_file >> all_families.fam
    done
    awk '{print $0 "\t" 0}' all_families.fam > trios.fam

    echo "How many lines in original trios.fam/ how many probands"
    wc -l trios.fam

    cut -f2 trios.fam > probands.list

    # Loop through each sample in all_samples.list
    while read -r sample; do
        # Check if the sample exists in the second column of trios.fam
        if grep -q "$sample" probands.list; then
            :
        else
            # Append the new row if the sample is not found
            echo -e "0\t$sample\t0\t0\t0\t0" >> trios.fam
            #echo "$sample was added to trios.2.fam"
        fi
    done < ~{subjects_list}
    # Find the directory where the VCF file is located
    vcf_dir=$(dirname ~{vcf})

    
    java -Xmx4G -jar /usr/GenomeAnalysisTK.jar \
      -R ~{Reference_Fasta} \
      -T VariantAnnotator \
      -V ~{vcf} \
      -o ~{vcf_basename} \
      -A AS_InbreedingCoeff \
      --pedigree trios.fam

  >>>
  runtime {
    docker: "broadinstitute/gatk3:3.8-1"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    memory: "4 GB"
    preemptible: 3
  }
  output {
    File output_vcf = "~{vcf_basename}"
  }
}


task annotate_vcf_with_pab_max_expr {
  input {
    File vcf
    Int? disk_gb
  }
  String vcf_basename=basename("~{vcf}")
  String vcf_basename_bgz=basename("~{vcf}",".gz") + ".bgz"
  Int default_disk_gb = ceil(6 * size(vcf,"GB")) + 15

  command <<<
  set -eu -o pipefail

  # Load the necessary Python3 environment if required, or set up Hail
  python3 <<CODE
  import hail as hl

  # Step 1: Initialize Hail and import the VCF as a Hail Matrix Table
  hl.init()
  
  # Load the normalized VCF as a Hail Matrix Table
  mt = hl.import_vcf('~{vcf}',reference_genome='GRCh38', array_elements_required=False, force_bgz=True)

  # Step 2: Define the pab_max_expr function (as provided earlier)
  # Function is taken from https://broadinstitute.github.io/gnomad_methods/_modules/gnomad/utils/annotations.html#pab_max_expr
  def pab_max_expr(gt_expr: hl.expr.CallExpression,ad_expr: hl.expr.ArrayExpression) -> hl.expr.ArrayExpression:
    expr = hl.agg.array_agg(lambda x: hl.agg.filter(gt_expr.is_het(),hl.agg.max(hl.binom_test(x, hl.sum(ad_expr), 0.5, "two-sided"))),ad_expr[1:])
    return expr

  # Step 3: Apply the pab_max_expr function to the Hail Matrix Table
  mt = mt.annotate_rows(info=mt.info.annotate(pab_max_expr=pab_max_expr(mt.GT, mt.AD)))
  mt.describe()
  # Step 4: Add the new pab_max_expr to the VCF header
  # Step 4.1: When being read in by hail, a lot of fields loose their header.
  metadata = {
    'info': {
        'pab_max_expr': {
            'Description': 'P-value from binomial test for the maximum alternate allele read depth',
            'Number': 'A',
            'Type': 'Float'
        }
    },
    'filter': {
    'VQSRTrancheSNP99.00to99.90+': {
        'Description': 'Truth sensitivity tranche level for SNP model at VQS Lod < -10.0'
    },
    'VQSRTrancheINDEL99.95to100.00': {
        'Description': 'Truth sensitivity tranche level for INDEL model at VQS Lod: -120.4376 <= x < -5.206'
    },
    'VQSRTrancheINDEL99.00to99.50': {
        'Description': 'Truth sensitivity tranche level for INDEL model at VQS Lod: -1.242 <= x < -0.7187'
    },
    'VQSRTrancheINDEL99.50to99.90': {
        'Description': 'Truth sensitivity tranche level for INDEL model at VQS Lod: -4.0735 <= x < -1.242'
    },
    'VQSRTrancheINDEL99.90to99.95': {
        'Description': 'Truth sensitivity tranche level for INDEL model at VQS Lod: -5.206 <= x < -4.0735'
    },
    'VQSRTrancheINDEL99.95to100.00+': {
        'Description': 'Truth sensitivity tranche level for INDEL model at VQS Lod < -120.4376'
    }
    },
    'format': {
    'PL': {
        'Number': 'G',
        'Type': 'Integer',
        'Description': 'Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification'
    }
    }
  }

  # Step 5: Write the annotated VCF to a file
  hl.export_vcf(mt, '~{vcf_basename_bgz}', metadata=metadata)

  # Stop the Hail session
  hl.stop()
  CODE
  >>>

  output {
    File output_vcf = "~{vcf_basename_bgz}"
  }

  runtime {
    docker: "vanallenlab/hail_rf"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    memory: "4 GB"
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


task GnomadHWE_part1 {
  input {
    File vcf_slice
  }
 
  command <<<
    # Extract Essential Columns from gnomad
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%nhomalt\t%AC\t%AN\t%AF\n' ~{vcf_slice} > tmp1.tsv
    
    # Add the nhet column
    awk 'BEGIN {OFS="\t"} \
     {
         AC = $6;
         AN = $7;
         obs_nhomalt = $5;
         AF = $8;
         obs_nhet = AC - 2 * obs_nhomalt;
         obs_nhomref = AN - obs_nhomalt - obs_nhet;
         exp_nhomalt = AF * AF * AN;
         exp_nhomref = (1 - AF) * (1 - AF) * AN;
         exp_nhet = 2 * (1 - AF) * AF * AN;
         #t1 = ((obs_nhet - exp_nhet ) * (obs_nhet - exp_nhet)) / exp_nhet;
         #t2 = ((obs_nhomref - exp_nhomref) * (obs_nhomref - exp_nhomref)) / exp_nhomref;
         #t3 = ((obs_nhomalt - exp_nhomalt) * (obs_nhomalt - exp_nhomalt)) / exp_nhomalt;
         #P_HWE = t1 + t2 + t3;
         print $1,$2, obs_nhomref, obs_nhet, obs_nhomalt, exp_nhomref, exp_nhet, exp_nhomalt;
     }' tmp1.tsv > hwe.tsv

  >>>
  output {
    File hwe_part1_file = "hwe.tsv"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
  }
}

task GnomadHWE_part2 {
  input {
    File hwe_file
  }

  command <<<
  python3 <<CODE
  import sys
  import numpy as np

  # Adopted from https://csg.sph.umich.edu/abecasis/Exact/snp_hwe.c
  def snp_hwe(obs_homref, obs_het, obs_homalt):
    if obs_homref < 0 or obs_homalt < 0 or obs_het < 0:
        raise ValueError(f"Invalid input: obs_homref={obs_homref}, obs_het={obs_het}, obs_homalt={obs_homalt}")

    obs_homc = max(obs_homref, obs_homalt)
    obs_homr = min(obs_homref, obs_homalt)

    rare_copies = 2 * obs_homr + obs_het
    genotypes = obs_het + obs_homc + obs_homr

    # Check to prevent division by zero
    if genotypes == 0:
        return "No Variants Present"
        #raise ValueError("Genotypes sum to zero, invalid input.")

    # Allocate array for heterozygote probabilities
    het_probs = np.zeros(int(rare_copies) + 1)
    
    # Start at midpoint
    mid = rare_copies * (2 * genotypes - rare_copies) // (2 * genotypes)

    # Check to ensure that midpoint and rare alleles have same parity
    if (rare_copies % 2) != (mid % 2):
        mid += 1

    mid = int(mid)
    curr_hets = mid
    curr_homr = int((rare_copies - mid)) // 2
    curr_homc = int(genotypes - curr_hets - curr_homr)

    het_probs[mid] = 1.0
    total_sum = het_probs[mid]
    
    for curr_hets in range(mid, 1, -2):
        het_probs[curr_hets - 2] = (het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) /
                                    (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0)))
        total_sum += het_probs[curr_hets - 2]

        # Update counts for next iteration
        curr_homr += 1
        curr_homc += 1

    curr_hets = mid
    curr_homr = int((rare_copies - mid)) // 2
    curr_homc = int(genotypes - curr_hets - curr_homr)
    
    for curr_hets in range(mid, int(rare_copies) - 1, 2):
        if curr_hets + 2 <= rare_copies:
            het_probs[curr_hets + 2] = (het_probs[curr_hets] * 4.0 * curr_homr * curr_homc /
                                        ((curr_hets + 2.0) * (curr_hets + 1.0)))
            total_sum += het_probs[curr_hets + 2]

            # Update counts for next iteration
            curr_homr -= 1
            curr_homc -= 1

    # Normalize probabilities
    het_probs /= total_sum

    # p-value calculation for p_hwe
    p_hwe = np.sum(het_probs[:int(obs_het) + 1])

    return min(p_hwe, 1.0)

  def main(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write('CHR\tPOS\tOBS(HOM1/HET/HOM2)\tE(HOM1/HET/HOM2)\tP_HWE\n')
        for line in infile:
            parts = line.strip().split()
            if len(parts) < 8:
                continue  # Skip lines that don't have enough data

            chrom = parts[0]
            pos = parts[1]
            obs_nhomref = int(parts[2])
            obs_nhet = int(parts[3])
            obs_nhomalt = int(parts[4])
            exp_nhomref = float(parts[5])
            exp_nhet = float(parts[6])
            exp_nhomalt = float(parts[7])

            p_value = snp_hwe(obs_nhomref, obs_nhet, obs_nhomalt)

            # At times Gnomad has locations where there is no Variant Called at a locus but it still exists
            if p_value == "No Variants Present":
              continue   
            # Write to output file with p-value appended
            outfile.write('{}\t{}\t{}/{}/{}\t{}/{}/{}\t{:.4f}\n'.format(chrom, pos, obs_nhomref, obs_nhet, obs_nhomalt, exp_nhomref, exp_nhet, exp_nhomalt, p_value))

  if __name__ == "__main__":
    # Replace ~{hwe_file} with the actual path to your input file
    input_file = '~{hwe_file}'
    output_file = 'hwe.part2.tsv'
    main(input_file, output_file)
  CODE
  >>>
  output {
    File hwe_part2_file = "hwe.part2.tsv"
  }
  runtime {
    docker: "zatonovo/numpy:3.7"
    preemptible: 3
  }
}

task RunBcfTools {
  input {
    File vcf_slice
    String callset_name
  }

  command <<<
    # Ectract Allele Frequency
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' ~{vcf_slice} > ~{callset_name}.frq
  >>>

  runtime{
    disks: "local-disk 10 HDD"
    docker: "vanallenlab/bcftools"
    preemptible: 3
  }
  output{
    File af_file = "~{callset_name}.frq"
  }
}


task grab_rf_training_variants {
  input{
    File vcf
    Array[File] TP_reference_vcfs
    Array[File] TP_reference_vcf_idxs
    String FP_filter_criteria = "AS_QD < 2.0 || AS_FS > 60.0 || AS_MQ < 30.0"
  }
  Int default_disk_size = ceil(size(vcf, "GB") * 10) + ceil(size(TP_reference_vcfs, "GB")) + 25
  command <<<
  set -euxo pipefail

  bcftools index -t ~{vcf}

  # Define the file names
  fp_file="unique_fp_positions.txt"
  tp_file="unique_tp_positions.txt"

  # Check if the VCF has any variants
  variant_count=$(bcftools view -H ~{vcf} -G | wc -l)

  # If there are no variants, create the files and exit successfully
  if [[ "$variant_count" -eq 0 ]]; then
      touch "$fp_file" "$tp_file"
      echo "VCF file contains no variants. Created empty files: $fp_file, $tp_file"
      exit 0
  fi

  # Continue with further processing if variants exist
  echo "VCF file contains variants. Proceeding with analysis..."

  # Find range that vcf spans for efficient querying
  range=$(bcftools query -f '%CHROM\t%POS\n' ~{vcf} | awk 'NR==1{start=$2; chrom=$1} {end=$2} END{print chrom":"start"-"end}')

  # Step 1: Create VCF with false positives
  bcftools filter -i "~{FP_filter_criteria}" ~{vcf} -Oz -o false_positives.vcf.gz

  # Exclude positions where ALT contains '*'
  bcftools view -e 'ALT=="*"' false_positives.vcf.gz -Oz -o false_positives_no_star.vcf.gz

  bcftools query -f '%CHROM\_%POS\_%REF\_%ALT\n' false_positives_no_star.vcf.gz  | sort -u >> unique_fp_positions.txt
  echo "Unique false positive positions saved: $(wc -l < unique_fp_positions.txt) lines"
  rm false_positives.vcf.gz false_positives_no_star.vcf.gz
  echo "False positives VCF removed."

  # Step 2: Extract TP positions (chr, pos, ref, alt) from each TP reference VCF and create a union
  for tp_vcf in ~{sep=" " TP_reference_vcfs}; do
    # Step 2.1: Filter the HQ reference sites to the range we are interested in
    bcftools view -r "$range" "$tp_vcf" -Oz -o region.vcf.gz
    # Step 2.2: Remove annotation info just for the sake of having a faster runtime
    bcftools annotate -x INFO region.vcf.gz -Oz -o region.noINFO.vcf.gz
    # Step 2.3: Normalize the vcf
    bcftools norm -m -any region.noINFO.vcf.gz -Oz -o region.noINFO.norm.vcf.gz
    # Step 2.4: Filter to only variants we want
    bcftools index -t region.noINFO.norm.vcf.gz
    bcftools isec -p isec_dir region.noINFO.norm.vcf.gz ~{vcf}

    rm region.vcf.gz
    ls isec_dir/
    bcftools query -f '%CHROM\_%POS\_%REF\_%ALT\n' isec_dir/0002.vcf >> all_tp_positions.txt
    rm region.noINFO.vcf.gz
    rm region.noINFO.norm.vcf.gz
    rm -r isec_dir
  done
  echo "All TP positions extracted and saved: $(wc -l < all_tp_positions.txt) lines"

  # Sort and remove duplicates from the combined positions
  sort -u all_tp_positions.txt > unique_tp_positions.txt
  echo "Unique TP positions saved: $(wc -l < unique_tp_positions.txt) lines"
  rm all_tp_positions.txt

  # Step 3: Remove overlapping lines from unique_fp_positions.txt
  # Create a temporary file to store the non-overlapping lines
  comm -23 "$fp_file" "$tp_file" > "${fp_file}.tmp"

  
  # Replace the original unique_fp_positions.txt with the non-overlapping lines
  mv "${fp_file}.tmp" "$fp_file"
  echo "Overlapping lines removed from FP positions. Remaining FP positions: $(wc -l < "$fp_file") lines"
  
  # Step 2: Count the number of lines in both files
  fp_count=$(wc -l < "$fp_file")
  tp_count=$(wc -l < "$tp_file")
  echo "FP count: $fp_count, TP count: $tp_count"

  >>>

  output {
    File false_positive_list = "unique_fp_positions.txt"
    File true_positive_list = "unique_tp_positions.txt"
  }
  runtime{
    docker: "vanallenlab/bcftools"
    memory: "4GB"
    disks: "local-disk ~{default_disk_size} HDD"
    preemptible: 3
  }
}

task train_hail_RF {
  input {
    File hail_table
  }
  
  command <<<
  set -eu -o pipefail
  python3 <<CODE
  import sys
  sys.path.append('/opt')
  import hail as hl
  import apply_rf_model as rf
  from datetime import datetime

  hl.init()

  # Define the desired types for columns
  desired_types = {
    'allele_type': hl.tstr,
    'qd': hl.tfloat64,
    'MQRankSum': hl.tfloat64,
    'ReadPosRankSum': hl.tfloat64,
    'SOR': hl.tfloat64,
    'pab_max_expr': hl.tfloat64,
    'Inbreeding_Coeff': hl.tfloat64,
    'has_star': hl.tbool,
    'n_alt_alleles': hl.tint64,
    'variant_type': hl.tstr,
    'was_mixed': hl.tbool,
    'TRAIN_LABEL': hl.tstr  # Import directly as boolean
  }

  # Import the false positives VCF
  ht = hl.import_table("~{hail_table}", types=desired_types)

  # Get the current date and time
  now = datetime.now()
  formatted_datetime = now.strftime("%Y-%m-%d %H:%M:%S")
  print("Current date and time for creating Training Hail Matrix:", formatted_datetime)

  #ht = mt.row()
  #ht = mt.to_table()
  ht.describe()
  print("Converted to Hail Table!!") 
  
  label = "TRAIN_LABEL"
  features = ["allele_type", "pab_max_expr","Inbreeding_Coeff", "qd","MQRankSum","ReadPosRankSum","SOR", "has_star","n_alt_alleles","variant_type","was_mixed"]
  rf_model = rf.train_rf(ht,features,label)
  rf.save_model(rf_model, out_path="rf.model", overwrite=True)
  CODE
  tar -czvf rf.model.tar.gz rf.model
  >>>
  runtime{
    docker:"vanallenlab/hail_rf:latest"
    disks: "local-disk 100 HDD"
    cpu: 4
    memory: "8G"
  }
  output{
    File rf_model = "rf.model.tar.gz"
  }
}

task prepare_training_vcfs_part2 {
  input {
    File vcf
    Float mean_pab_max_expr
    Float mean_AS_InbreedingCoeff
    Float mean_QD
    Float mean_MQRankSum
    Float mean_ReadPosRankSum
    Float mean_SOR
  }
  
  command <<<
  set -eu -o pipefail
  python3 <<CODE
  import sys
  sys.path.append('/opt')
  import hail as hl
  from datetime import datetime

  hl.init()

  # Import the false positives VCF
  print("Importing VCF")
  mt = hl.import_vcf("~{vcf}", reference_genome='GRCh38', array_elements_required=False, force_bgz=True)
  print("Successfully imported vcf into a matrix table")

  # Group by locus to aggregate information across multiallelic sites
  grouped_mt = mt.group_rows_by(mt.locus).aggregate(
    num_snvs=hl.agg.count_where(hl.is_snp(mt.alleles[0], mt.alleles[1])),
    num_indels = hl.agg.count_where(hl.is_indel(mt.alleles[0], mt.alleles[1])),
    was_mixed= (
        hl.agg.any(hl.is_indel(mt.alleles[0], mt.alleles[1])) & 
        hl.agg.any(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    ),
    variant_type = hl.case()
        .when(
            hl.agg.any(hl.is_indel(mt.alleles[0], mt.alleles[1]))
            & hl.agg.any(hl.is_snp(mt.alleles[0], mt.alleles[1])), 
            'mixed'
        ) 
        .when(hl.agg.count_where(hl.is_snp(mt.alleles[0], mt.alleles[1])) > 1, 'multi-SNV')
        .when(hl.agg.count_where(
            hl.is_indel(mt.alleles[0], mt.alleles[1])) > 1, 'multi-Indel')
        .when(hl.agg.count_where(hl.is_snp(mt.alleles[0], mt.alleles[1])) == 1, 'SNV')
        .when(
            hl.agg.count_where(
                hl.is_indel(mt.alleles[0], mt.alleles[1])) == 1, 'Indel'
        )
        .or_missing(),
    n_alt_alleles=hl.agg.count(),
    has_star=hl.agg.any(mt.alleles.contains('*'))
  )

  grouped_ht = grouped_mt.entries()
  grouped_ht = grouped_ht.key_by('locus')

  grouped_ht.describe()
  # Annotate with allele-specific features
  mt = mt.annotate_rows(
    allele_type=mt.info.VARIANT_TYPE,
    qd=mt.info.AS_QD[0],
    MQRankSum=mt.info.AS_MQRankSum[0],
    ReadPosRankSum=mt.info.AS_ReadPosRankSum[0],
    SOR=mt.info.AS_SOR[0],
    pab_max_expr=mt.info.pab_max_expr[0],
    Inbreeding_Coeff=mt.info.AS_InbreedingCoeff[0],
    has_star=hl.bool(grouped_ht[mt.locus].has_star),
    n_alt_alleles=grouped_ht[mt.locus].n_alt_alleles,
    variant_type=grouped_ht[mt.locus].variant_type,
    was_mixed=hl.bool(grouped_ht[mt.locus].was_mixed)
  )
  ht = mt.rows()

  ht = ht.select(TRAIN_LABEL=ht.info.TRAIN_LABEL, 
    allele_type=ht.allele_type,
    pab_max_expr=ht.pab_max_expr, 
    Inbreeding_Coeff=ht.Inbreeding_Coeff,
    qd=ht.qd, 
    MQRankSum=ht.MQRankSum,
    ReadPosRankSum=ht.ReadPosRankSum, 
    SOR=ht.SOR,
    has_star=ht.has_star, 
    n_alt_alleles=ht.n_alt_alleles,
    variant_type=ht.variant_type, 
    was_mixed=ht.was_mixed
  )

  # Fill Imputed Values for the Hail Table: ht   
  ## List of features and their respective imputation values
  features_to_impute = {
      "pab_max_expr": ~{mean_pab_max_expr},
      "Inbreeding_Coeff": ~{mean_AS_InbreedingCoeff},
      "qd": ~{mean_QD},
      "MQRankSum": ~{mean_MQRankSum},
      "ReadPosRankSum": ~{mean_ReadPosRankSum},
      "SOR": ~{mean_SOR}
  }

  # Annotating HAIL TABLE (THE WHOLE DATASET - NOT JUST TRAINING VARIANTS)
  # Annotate table with imputed values for missing data
  print("Checkpoint: 2")
  ht = ht.annotate(**{f: hl.or_else(ht[f], impute_val) for f, impute_val in features_to_impute.items()})

  # Annotate global values for reference (optional)
  print("Checkpoint: 3")
  ht = ht.annotate_globals(imputation_values=features_to_impute)

  print("Checkpoint: 4")
  train_ht = ht.filter(hl.is_defined(ht.TRAIN_LABEL) & ((ht.TRAIN_LABEL == "True") | (ht.TRAIN_LABEL == "False")))

  # Annotate global values for reference (optional)
  print("Checkpoint: 5")
  train_ht = train_ht.annotate_globals(imputation_values=features_to_impute)

  print("Checkpoint: 6")
  ht.export("ufc_wgs.tsv.bgz")
  train_ht.export("ufc_wgs.train.tsv.bgz")
  CODE
  >>>
  runtime{
    docker:"vanallenlab/hail_rf:latest"
    disks: "local-disk 20 HDD"
    memory: "8 GB"
    preemptible: 3
  }
  output{
    File hail_table = "ufc_wgs.tsv.bgz"
    File hail_train_table = "ufc_wgs.train.tsv.bgz"
    #File test_file = "table1.tsv.bgz"
  }
}

task concatenate_hail_train_tables {
  input {
    Array[File] hail_tables   # Array of Hail Tables to concatenate
  }

  command <<<
  set -euxo pipefail

  python3 <<CODE
  import hail as hl

  # Input: Array of Hail Table file paths
  hail_tables = "~{sep=' ' hail_tables}".split()

  # Load a Hail Table
  concatenated_table = hl.import_table(hail_tables[0])

  # Concatenate the tables
  for table_path in hail_tables[1:]:
    table = hl.import_table(table_path)
    concatenated_table = concatenated_table.union(table)
    # Optionally, you can delete 'table' explicitly to free up memory after each iteration
    del table

  # Write the concatenated table to the specified output file
  concatenated_table.export("ufc_wgs.tsv.bgz")
  CODE
  >>>

  output {
    File concatenated_table = "ufc_wgs.tsv.bgz"
  }

  runtime {
    docker: "vanallenlab/hail_rf"
    memory: "16G"
    preemptible: 3
  }
}

task concatenate_hail_train_tables_v2 {
  input {
    Array[File] hail_tables   # Array of Hail Tables to concatenate
  }

  command <<<
  set -euxo pipefail

  # Flag to track if the header has been added
  header_added=false

   # Iterate over the array of files
   for i in ~{sep=' ' hail_tables}; do
     if [ "$header_added" = false ]; then
       # For the first file, include the header
       bgzip -dc "$i" > concatenated_output.tsv
       header_added=true
     else
       # For subsequent files, skip the header
       bgzip -dc "$i" | tail -n +2 >> concatenated_output.tsv
     fi

     # Remove the processed file to save disk space
     rm "$i"
   done

  # Compress the concatenated file
  bgzip -c concatenated_output.tsv > ufc_wgs.train.tsv.bgz

  >>>

  output {
    File concatenated_table = "ufc_wgs.train.tsv.bgz"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    #memory: "8G"
    disks: "local-disk 50 HDD"
    preemptible: 3
  }
}


task concatenateFiles {
  input {
    Array[File] files
    String callset_name
  }
  command <<<
    cat ~{sep=" " files} > ~{callset_name}.tsv
  >>>
  runtime{
    disks: "local-disk 400 HDD"
    docker:"ubuntu:latest"
    preemptible: 0
  }
  output{
    File out = "~{callset_name}.tsv"
  }
}
