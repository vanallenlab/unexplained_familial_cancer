# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_5_GSEA {
  input {
    String step_10_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_8_VISUALIZE_VEP/sharded_vcfs"
    String cancer_type
    String analysis_5_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_5_GSEA/"
    String workspace_bucket = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
    String biological_pathway = "Mismatch_Repair"
    Array[String] allowed_consequences = ["frameshift_variant","stop_gained"]
  }
  File sample_data = "gs://~{workspace_bucket}/UFC_REFERENCE_FILES/analysis/~{cancer_type}/~{cancer_type}.metadata"
  Array[String] genes_of_interest = read_lines("gs://~{workspace_bucket}/UFC_REFERENCE_FILES/gene_lists/" + biological_pathway + ".list")
  Int negative_shards = 21

  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call Tasks.list_files_from_directory {
    input:
      dir = step_10_output_dir,
      suffix = ".tsv.gz" 
  }

  scatter (i in range(length(list_files_from_directory.out1) - negative_shards)){
    call T1_get_rows {
      input:
        genes_of_interest = genes_of_interest,
        allowed_consequences = allowed_consequences,
        variant_tsv = list_files_from_directory.out1[i]

    }
  }
  call Tasks.concatenateFiles {
    input:
      files = T1_get_rows.out1,
      output_name = biological_pathway
  }
}

task T1_get_rows {
  input {
    File variant_tsv
    Array[String] genes_of_interest
    Array[String] allowed_consequences
  }
  command <<<
  python3 <<CODE
  import pandas as pd

  df = pd.read_csv("~{variant_tsv}",sep='\t',index_col=False)
  genes = "~{sep=' ' genes_of_interest}".split()
  consequences = "~{sep=' ' allowed_consequences}".split()

  # Split 'gene_consequence' into gene and consequence
  df[['gene', 'consequence']] = df['gene_consequence'].str.split('_', n=1, expand=True)

  # Filter rows where both gene and consequence match
  df = df[df['gene'].isin(genes) & df['consequence'].isin(consequences)]
  df.to_csv("out.tsv",sep='\t',index=False)
  CODE
  gzip out.tsv
  >>>
  output{
    File out1 = "out.tsv.gz"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    disk: "local-disk 4 HDD"
    preemptible: 3
    cpu: 1
  }
}

task T2_gsea {
  input {
    File variant_counts
    Array[String] gene_list
    Array[String] damaging_terms
  }

  command <<<
  set -euo pipefail

  # Save gene list to file
  printf "%s\n" ~{sep='\n' gene_list} > genes.txt

  # Filter variant count matrix to just rows with matching genes
  gunzip -c ~{variant_counts} > all_variants.tsv
  grep -Ff genes.txt all_variants.tsv > filtered.tsv || true

  # Prepend header line (first line of original file)
  head -n 1 all_variants.tsv > filtered_with_header.tsv
  cat filtered.tsv >> filtered_with_header.tsv

  # Now parse using Python
  python3 <<CODE
  import pandas as pd

  # Load the filtered data
  df = pd.read_csv("filtered_with_header.tsv", sep="\t")

  # Damage-related consequences
  damaging = set(~{sep=', ' damaging_terms})

  # Extract gene and consequence from the first column
  df[['gene', 'consequence']] = df.iloc[:, 0].str.split('_', n=1, expand=True)

  # Keep only rows with damaging consequences
  df = df[df['consequence'].isin(damaging)]

  # Drop metadata columns, sum across rows for each patient
  patient_counts = (df.iloc[:, 1:-2] > 0).astype(int).sum()

  # Output to file
  patient_counts.reset_index().to_csv("damaging_counts.tsv", sep="\t", header=False, index=False)
  CODE
  >>>

  output {
    File out1 = "damaging_counts.tsv"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "4G"
    disk: "local-disks 10 HDD"
  }
}
