# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


version 1.0

workflow STEP_1_Make_Batches {
  input {
    File intake_qc_file
    File make_batches_script
    File global_qc_thresholds
    File batch_qc_thresholds
    File ids_to_keep_file
    File ufc_samples
  }
  call filter_discordant_sex {
    input:
      input_file = intake_qc_file,
      ids_to_keep_file = ids_to_keep_file,
      ufc_samples = ufc_samples
  }    
  call make_batches {
    input:
      filtered_intake_qc = filter_discordant_sex.filtered_file,
      make_batches_script = make_batches_script,
      global_qc_thresholds = global_qc_thresholds,
      batch_qc_thresholds = batch_qc_thresholds,
      ids_to_keep_file = ids_to_keep_file
  }
  call generate_sample_map {
    input:
      post_batching_intake_file = make_batches.post_qc
  }
  call graph_cancer_frequencies {
    input:
      post_batching_intake_file = make_batches.post_qc
  }
}

task plot_metrics {
  input{
    File post_batching_tsv
  }
  command <<<
  python3 /opt/STEP_1_Make_Batches.plot_qc.py --pre-filtered ~{post_batching_tsv} --post-filtered ~{post_batching_tsv}
  >>>
  output {
    File graphs = glob("*.png")
  }
  runtime {
    docker: "vanallenlab/g2c_ufc"
  }
}

task generate_sample_map {
  input {
    File post_batching_intake_file
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  df = pd.read_csv("~{post_batching_intake_file}",sep='\t')
  print(f"Length of Intake File: {len(df)}")
  df = df[df['global_qc_pass'] == True]
  print(f"Length of Intake File (passed): {len(df)}")

  # Define base paths for each cohort
  base_paths = {
      'aou': 'gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45/dfci-g2c-inputs/aou/gatk-hc/reblocked',
      'ceph': 'gs://dfci-g2c-inputs/ceph/gatk-hc/reblocked',
      'mesa': 'gs://dfci-g2c-inputs/mesa/gatk-hc/reblocked',
      'ufc': 'gs://dfci-g2c-inputs/ufc/gatk-hc/reblocked'
  }

  # Initialize a list to collect the rows for the final concatenated file
  all_rows = []

  # Iterate through the rows of the dataframe
  for index, row in df.iterrows():
      cohort = row['cohort'].lower()  # Assuming 'Cohort' column exists
      sample_name = row['original_id']  # Assuming 'original_id' column exists

      # Generate the full path only for valid cohorts in the base_paths dictionary
      if cohort in base_paths:
          gs_path = f"{base_paths[cohort]}/{sample_name}.reblocked.g.vcf.gz"
          # Add to the final list for concatenation
          all_rows.append(f"{sample_name}\t{gs_path}")

  # Write the concatenated output to a single file
  with open("ufc_sample_map.tsv", "w") as f:
      for row in all_rows:
          f.write(row + "\n")
  CODE 
  >>>

  output{
    File sample_map = "ufc_sample_map.tsv"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task graph_cancer_frequencies {
  input {
    File post_batching_intake_file
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  import matplotlib.pyplot as plt
  from collections import Counter
  import seaborn as sns

  # Function to split cancer types and count frequencies by the chosen feature
  def count_cancers_by_feature(df, feature):
      # Split 'cancer' column by semicolon and flatten the list
      df = df.assign(cancer=df['cancer'].str.split(';')).explode('cancer')  # Explode the list into separate rows
    
      # Count frequencies of each cancer type, grouped by the feature (e.g., sex, ancestry, cohort)
      grouped_counts = df.groupby(['cancer', feature]).size().unstack(fill_value=0)
      return grouped_counts

  # Plot function for bar charts with cancer types on x-axis, colored by another feature (sex, ancestry, etc.)
  def plot_cancer_frequencies(df, feature, output_file):
      cancer_counts = count_cancers_by_feature(df, feature)

      # Get distinct colors based on the number of unique values in the feature
      num_categories = len(cancer_counts.columns)
      color_palette = sns.color_palette("bright", num_categories)  # Use a diverse color palette

      # Plot stacked bar chart with cancer types on x-axis and feature (e.g., sex, ancestry) as colors
      fig, ax = plt.subplots(figsize=(10, 6))
      cancer_counts.plot(kind='bar', stacked=True, ax=ax, color=color_palette)

      plt.title(f'Cancer Type Frequency Colored by {feature}')
      plt.ylabel('Frequency')
      plt.xlabel('Cancer Type')
      plt.legend(title=feature, bbox_to_anchor=(1.05, 1), loc='upper left')
      plt.tight_layout()

      # Save the figure
      plt.savefig(output_file)
      plt.show()  # Close the figure after saving to avoid memory issues

  # Example usage (replace 'your_data.csv' with your actual data)
  df = pd.read_csv("~{post_batching_intake_file}", sep='\t')
  df = df[(df['cancer'] != "control") & (df['cancer'] != "unknown")]

  # Generate and save plots for 'cohort', 'inferred_sex', and 'grafpop_ancestry'
  plot_cancer_frequencies(df, 'cohort', "cancer_cohort_occurence.png")
  plot_cancer_frequencies(df, 'inferred_sex', "cancer_sex_occurence.png")
  plot_cancer_frequencies(df, 'grafpop_ancestry', "cancer_ancestry_occurence.png")
  CODE
  >>>
  runtime {
    docker: "vanallenlab/pydata_stack"
  }
  output {
    Array[File] graphs = glob("*.png")
  }
} 

task make_batches {
  input {
    File filtered_intake_qc
    File make_batches_script
    File global_qc_thresholds
    File batch_qc_thresholds
    File ids_to_keep_file
  }
  command <<<
  set -eu -o pipefail
  python3 ~{make_batches_script} ~{filtered_intake_qc} \
          --match-on batching_sex --match-on batching_pheno \
          --batch-by wgd_score --batch-by median_coverage \
          --global-qc-cutoffs ~{global_qc_thresholds} \
          --batch-qc-cutoffs ~{batch_qc_thresholds} \
          --custom-qc-pass-samples ~{ids_to_keep_file} \
          --prefix dfci-ufc.gatk-sv \
          --logfile dfci-g2c.intake_qc.all.post_qc_batching.log \
          -o dfci-ufc.post_qc.all.tsv
  >>>
  output {
    File post_qc = "dfci-ufc.post_qc.all.tsv"
    File logfile = "dfci-g2c.intake_qc.all.post_qc_batching.log"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    read_from_cache: false
  }
}

task filter_discordant_sex {
    input {
        File input_file         # The input TSV file
        File ids_to_keep_file   # File containing G2C_ids to keep (one per line)
        File ufc_samples
    }

    command <<<
    head -n 1 ~{input_file} > ufc_input.tsv
    grep -E 'ceph|mesa|aou|ufc' ~{input_file} >> ufc_input.tsv
    python3 <<CODE
    import pandas as pd

    # Read the input data
    df = pd.read_csv("ufc_input.tsv", sep="\t")

    # Read the list of IDs from the file and ensure they are treated as strings
    with open('~{ufc_samples}') as f:
        ufc_list = [line.strip() for line in f]
    print(f"Length of UFC list: {len(ufc_list)}")

    # Ensure that both the 'original_id' column and the IDs in the list are strings
    df['original_id'] = df['original_id'].astype(str).str.strip()
    ufc_list = [str(id).strip() for id in ufc_list]

    # Filter the DataFrame in place
    filtered_df = df[df['original_id'].isin(ufc_list)]
    df = filtered_df
    print(f"Length of UFC DF: {len(df)}")

    # Read the IDs to keep
    with open("~{ids_to_keep_file}", 'r') as f:
        ids_to_keep = [line.strip() for line in f.readlines()]
    print(f"Length of Sex Chromosome Exceptions list: {len(ids_to_keep)}")

    # Filter rows where inferred_sex and sex_karyotype are concordant
    filtered_df = df[((df['inferred_sex'] == 'male') & (df['sex_karyotype'] == 'XY')) | 
                     ((df['inferred_sex'] == 'female') & (df['sex_karyotype'] == 'XX'))]

    # Keep discordant rows for specific G2C_ids
    to_keep_df = df[df['G2C_id'].isin(ids_to_keep)]

    # Combine the two DataFrames
    final_df = pd.concat([filtered_df, to_keep_df])
    print(f"Length of Final df: {len(final_df)}")

    # Drop duplicates to ensure each G2C_id appears only once
    final_df = final_df.drop_duplicates(subset='G2C_id')

    # Write the result to an output file
    final_df.to_csv("out.tsv", sep="\t", index=False)
    CODE
    >>>

  runtime {
    docker: "vanallenlab/pydata_stack"  # Or any image that has Python and pandas
  }

  output {
    File filtered_file = "out.tsv"
  }
}

