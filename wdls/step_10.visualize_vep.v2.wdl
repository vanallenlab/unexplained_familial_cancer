# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow STEP_10_VISUALIZE_VEP {
  input{
    String step_9_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/sharded_vcfs" # Directory to STEP_9 Output VCFs
    File aou_subjects
    File gene_list
    String output_name
    String gene_list_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/gene_lists/"
    String step_10_storage_dir = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
    String step_10_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_10_VISUALIZE_VEP/stats/" 
  }

  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call Tasks.gather_vcfs as gather_vcfs{
    input:
      dir = step_9_output_dir
  }

  call Tasks.sort_vcf_list {
    input:
      unsorted_vcf_list = gather_vcfs.vcf_list
  }

  call Tasks.list_files_from_directory {
    input:
      dir = gene_list_dir
  }

  Int negative_shards = 0

  scatter (i in range(length(sort_vcf_list.vcf_arr)-negative_shards)){
    call Convert_To_TSV {
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        aou_subjects = aou_subjects,
        gene_list = gene_list
    }

    call Filter_Vep_TSV{
      input:
        input_tsv = Convert_To_TSV.out1
    }

    call Count_Shard_Mutations{
      input:
        zipped_tsvs = Filter_Vep_TSV.out1,
        shard = i,
        aou_subjects = aou_subjects
    }
  }
  
  call Tasks.sum_tables_by_sample{
    input:
      input_tables = Count_Shard_Mutations.out1
  }

  call Tasks.append_covariates{
    input:
      data = sum_tables_by_sample.out1,
      output_name = output_name
  }

  call Tasks.copy_file_to_storage{
    input:
      text_file = append_covariates.out1,
      output_dir = step_10_output_dir
  }
}




task Count_Shard_Mutations {
  input {
    File zipped_tsvs
    Int shard
    File aou_subjects
  }

  Int default_mem_gb = ceil(size(zipped_tsvs,"GB")) * 10 + 2
  String zipped_tsvs_basename = basename(zipped_tsvs)
  
  command <<<
  set -euxo pipefail
  mv ~{zipped_tsvs} .
  unzip ~{zipped_tsvs_basename}
  
  python3 <<CODE
  import pandas as pd
  import os
  import gc

  # Step 1: Read the sample IDs from 'subjects.list' file
  with open("~{aou_subjects}") as f:
      samples = [line.strip() for line in f]
  
  master_df = pd.DataFrame(index=samples)
  master_df.index.name = 'Sample'
 
  def count_shards(file, master_df):
    print(f"Currently processing: {file}")
    variant_type = file.split('/')[-1].split('.')[0]
    df = pd.read_csv(f"tsv_files/{file}",sep='\t',index_col=False)

    sample_scores = df[samples].sum().to_dict()

    # Convert the scores dictionary to a DataFrame for output
    scores_df = pd.DataFrame(list(sample_scores.items()), columns=['Sample', variant_type])
 
    master_df = master_df.merge(scores_df, how='left', on='Sample')
    del df
    del scores_df
    del sample_scores
    gc.collect()
    return master_df


  # Get all files in the current directory ending with .tsv
  files = [f for f in os.listdir('tsv_files/') if f.endswith('.tsv')]

  # Process each file
  for file in files:
      # Run your code on each file
      print(f"Processing file: {file}")
      # You can call your function here with the file name
      master_df = count_shards(file, master_df)

  #master_df = master_df[['Sample'] + sorted(c for c in master_df.columns if c != 'Sample')]
  master_df.to_csv("Variant_Count.~{shard}.tsv",sep='\t',index=False) 
  CODE
  >>>
  output {
    File out1 = "Variant_Count.~{shard}.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    memory: "~{default_mem_gb}GB"
    preemptible: 3
  }
}

task SumPatientScores {
  input {
    Array[File] input_files
    String output_file = "summed_scores.tsv"
    String variable_name
  }

  command <<<
  # Run a Python script to sum the scores
  python3 <<CODE
  import pandas as pd
  import sys

  # Initialize an empty DataFrame to accumulate results
  result_df = pd.DataFrame(columns=["Sample", "~{variable_name}"])

  # Loop through each input file
  for file in ["~{sep='\",\"' input_files}"]:
      # Read the TSV file
      df = pd.read_csv(file, sep='\t')

      # If there are no columns, skip
      if df.empty:
          continue

      # Accumulate the scores by Sample
      for index, row in df.iterrows():
          sample = row['Sample']
          score = row['~{variable_name}']

          # Add or update the score for this sample
          if sample in result_df['Sample'].values:
              result_df.loc[result_df['Sample'] == sample, '~{variable_name}'] += score
          else:
              # Use pd.concat() to append new rows
              new_row = pd.DataFrame({'Sample': [sample], '~{variable_name}': [score]})
              result_df = pd.concat([result_df, new_row], ignore_index=True)

  # Write the accumulated results to an output file
  result_df.to_csv("~{output_file}", sep='\t', index=False)
  CODE
  >>>

  output {
    File out = "~{output_file}"  # Output file with summed scores
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task Filter_Vep_TSV {
    input {
      File input_tsv
    }
    Int default_mem_gb = ceil(size(input_tsv,"GB")) * 20 + 4
    String unzipped_tsv = basename(input_tsv, ".gz")
    String zipped_tsv = basename(input_tsv)
    String output_file_basename = basename(input_tsv, ".tsv.gz") 
    command <<<
    set -x
    mv ~{input_tsv} .
    gunzip ~{zipped_tsv}


    echo -e "This many lines in the file"
    wc -l ~{unzipped_tsv}
    echo -e "This file is this large"
    du -sh ~{unzipped_tsv}
    python3 <<CODE
    import pandas as pd
    import numpy as np
    from collections import defaultdict
    import gc

    dtype_map = {
        "ID": "string",
        "AF": "float32",
        "IMPACT": "category",
        "clinvar_clnsig": "category",
        "Consequence": "category",
        "REVEL_score": "string",
        "SYMBOL": "string",
    }


    # Read in the TSV file
    df = pd.read_csv("~{unzipped_tsv}", sep='\t', index_col=False, dtype=dtype_map)

    # Identify gnomAD columns and convert them to string data type
    gnomad_cols = [col for col in df.columns if col.startswith('gnomAD')]
    df[gnomad_cols] = df[gnomad_cols].astype(str)

    # Keep only columns in dtype_map or gnomad_cols
    #df = df[[col for col in df.columns if col in dtype_map or col in gnomad_cols]]

    print(f"Length of pandas df: {len(df)}")

    # Group by 'ID', and join 'Consequence' with '&'
    df = df.groupby("ID", as_index=False).agg({
        **{col: 'first' for col in df.columns if col != 'Consequence' and col != 'ID'},
        'Consequence': lambda x: '&'.join(x.dropna().astype(str).unique())
    })

    def process_revel_score(score_field: str) -> float:
        score_field = str(score_field)
        if not score_field or score_field == ".":
            return 0.0
        parts = score_field.split("&")
        numeric_values = []
        for p in parts:
            try:
                if p != ".":
                    numeric_values.append(float(p))
            except ValueError:
                continue
        return max(numeric_values) if numeric_values else 0.0

    # Parse REVEL scores once
    df["REVEL_score"] = df["REVEL_score"].astype(str)
    df["REVEL_score_numeric"] = df["REVEL_score"].apply(process_revel_score).astype("float32")

    def filter_rare_variants(df):

        low_cols = [
            'gnomAD_AF_non_cancer_afr', 'gnomAD_AF_non_cancer_amr', 'gnomAD_AF_non_cancer_eas',
            'gnomAD_AF_non_cancer_fin', 'gnomAD_AF_non_cancer_nfe', 'gnomAD_AF_non_cancer_sas'
        ]
        med_cols = ['AF']
        high_cols = [
            'gnomAD_AF_non_cancer_asj', 'gnomAD_AF_non_cancer_mid',
            'gnomAD_AF_non_cancer_ami', 'gnomAD_AF_non_cancer_oth'
        ]

        all_cols = low_cols + med_cols + high_cols
        existing_cols = [col for col in all_cols if col in df.columns]
        if not existing_cols:
            return df.iloc[0:0]

        low_mask = np.ones(len(df), dtype=bool)
        for col in low_cols:
            if col in df.columns:
                values = df[col].astype(str).str.split(',').str[0]
                values = pd.to_numeric(values, errors='coerce')
                low_mask &= (values < 0.01).fillna(False)

        med_mask = np.ones(len(df), dtype=bool)
        for col in med_cols:
            if col in df.columns:
                values = df[col].astype(str).str.split(',').str[0]
                values = pd.to_numeric(values, errors='coerce')
                low_mask &= (values < 0.05).fillna(False)

        high_mask = np.ones(len(df), dtype=bool)
        for col in high_cols:
            if col in df.columns:
                values = df[col].astype(str).str.split(',').str[0]
                values = pd.to_numeric(values, errors='coerce')
                high_mask &= (values < 0.1).fillna(False)

        return df[low_mask & high_mask]


    # Define filters
    filters = {
        "high_impact": lambda df: df["IMPACT"].str.contains("HIGH", na=False),
        "plp_clinvar": lambda df: df["clinvar_clnsig"].str.contains(
          r"^(Pathogenic|Likely_pathogenic|Pathogenic/Likely_pathogenic)$", case=False, na=False),
        "stopgained": lambda df: df["Consequence"].str.contains("stop_gained", na=False),
        "missense": lambda df: df["Consequence"].str.contains("missense", na=False),
        "frameshift": lambda df: df["Consequence"].str.contains("frameshift", na=False),
        "intergenic": lambda df: df["Consequence"].str.contains("intergenic_variant", na=False),
        "intronic": lambda df: df["Consequence"].str.contains("intron_variant", na=False),
        "synonymous": lambda df: df["Consequence"].str.contains("synonymous", na=False),
        "splice_donor": lambda df: df["Consequence"].str.contains("splice_donor", na=False),
        "splice_acceptor": lambda df: df["Consequence"].str.contains("splice_acceptor", na=False),
        "5_prime_UTR_variant": lambda df: df["Consequence"].str.contains("5_prime_UTR_variant", na=False),
        "3_prime_UTR_variant": lambda df: df["Consequence"].str.contains("3_prime_UTR_variant", na=False),
        "transcript_ablation": lambda df: df["Consequence"].str.contains("transcript_ablation", na=False),
        "TF_binding_site_variant": lambda df: df["Consequence"].str.contains("TF_binding_site_variant", na=False),
        "splice_polypyrimidine_tract_variant": lambda df: df["Consequence"].str.contains("splice_polypyrimidine_tract_variant", na=False),
        "inframe_insertion":lambda df: df["Consequence"].str.contains("inframe_insertion", na=False),
        "inframe_deletion":lambda df: df["Consequence"].str.contains("inframe_deletion", na=False),
         "splice_region_variant":lambda df: df["Consequence"].str.contains("splice_region_variant", na=False),
        "incomplete_terminal_codon_variant": lambda df: df["Consequence"].str.contains("incomplete_terminal_codon_variant", na=False),
        "synonymous_variant": lambda df: df["Consequence"].str.contains("synonymous_variant", na=False),
        "stop_lost": lambda df: df["Consequence"].str.contains("stop_lost", na=False),
        "start_lost": lambda df: df["Consequence"].str.contains("start_lost", na=False),
        "non_coding_transcript_exon_variant": lambda df: df["Consequence"].str.contains("non_coding_transcript_exon_variant", na=False),
        "transcript_amplification": lambda df: df["Consequence"].str.contains("transcript_amplification"),
        "start_retained_variant":lambda df: df["Consequence"].str.contains("start_retained_variant"),
        "stop_retained_variant": lambda df: df["Consequence"].str.contains("stop_retained_variant"),
        "mature_miRNA_variant": lambda df: df["Consequence"].str.contains("mature_miRNA_variant"),
        "non_coding_transcript_exon_variant": lambda df: df["Consequence"].str.contains("non_coding_transcript_exon_variant"),
        "non_coding_transcript_exon_variant":lambda df: df["Consequence"].str.contains("non_coding_transcript_exon_variant"),
        "MODIFIER": lambda df: df["IMPACT"].str.contains("MODIFIER", na=False),
        "LOW": lambda df: df["IMPACT"].str.contains("LOW", na=False),
        "MODERATE": lambda df: df["IMPACT"].astype(str).str.contains("MODERATE", na=False),
        "incomplete_terminal_codon_variant": lambda df: df["Consequence"].str.contains("incomplete_terminal_codon_variant"),
        "REVEL050": lambda df: df["REVEL_score_numeric"] >= 0.5,
        "REVEL075": lambda df: df["REVEL_score_numeric"] >= 0.75,
        "incomplete_terminal_codon_variant": lambda df: df["Consequence"].str.contains("incomplete_terminal_codon_variant")    }

    with open("expected_num_files.txt", "w") as f:
        f.write(f"{len(filters) * 2 + 2}")

    df_grouped = df.groupby("ID", as_index=False).first()
    df_grouped.to_csv("total_variants.tsv",sep='\t',index=False)
    df_rare = filter_rare_variants(df_grouped)
    df_rare.to_csv("total_variants_rare.tsv",sep='\t',index=False)
   
    del df_rare
    del df_grouped
    gc.collect()
    # Apply filters and export
    for name, condition in filters.items():
        print(f"Starting {name} filtering.")
        mask = condition(df)
        df_filtered = df[mask]
        df_filtered = df_filtered.groupby("ID", as_index=False).first()
 
        print(f"{name}: {len(df_filtered)} variants")

        df_filtered.to_csv(f"{name}.tsv", sep='\t', index=False)

        df_rare = filter_rare_variants(df_filtered)
        df_rare.to_csv(f"{name}_rare.tsv", sep='\t', index=False)

        print(f"{name}_rare: {len(df_rare)} rare variants")

        # Explicitly delete temporary DataFrames to release memory
        del df_filtered
        del df_rare
        del mask

        # Perform garbage collection to free up memory
        gc.collect()

    CODE

    # Step 1: Create a directory to hold the .tsv files
    mkdir -p tsv_files

    # Step 2: Move all .tsv files from the home directory into that folder
    rm ~{unzipped_tsv}
    mv *.tsv tsv_files/

    expected=$(<expected_num_files.txt)
    actual=$(find tsv_files/ -maxdepth 1 -name "*.tsv" | wc -l)

    if [ "$expected" -eq "$actual" ]; then
      echo "Correct number of files: $actual"
    else
      echo "Mismatch: expected $expected but found $actual"
      exit 1
    fi

    # Step 3: Zip the folder into a file
    zip -r tsv_files.zip tsv_files
    >>>

    output {
      File out1 = "tsv_files.zip"
    }

    runtime {
      docker: "vanallenlab/g2c_pipeline"
      memory: "~{default_mem_gb}GB"
      disks: "local-disk 15 HDD"
      preemptible: 3
    }
}


task Convert_To_TSV {
  input {
    File vcf
    File aou_subjects
    File gene_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/gencode.v47.autosomal.protein_coding.genes.list"
  }
  String output_file = basename(vcf, ".vcf.bgz") + ".tsv"
  Int default_disk_gb = ceil(size(vcf, "GB")) * 4 + 10
  command <<<
  set -x
    bcftools view -S ~{aou_subjects} ~{vcf} -Oz -o aou.tmp.vcf.gz
    bcftools view -i 'AC > 0' aou.tmp.vcf.gz -Oz -o aou.tmp2.vcf.gz
    bcftools +fill-tags aou.tmp2.vcf.gz -Oz -o aou.vcf.gz -- -t AF
    bcftools index -t aou.vcf.gz

    rm ~{vcf} aou.tmp.vcf.gz aou.tmp2.vcf.gz
 
    if bcftools view -h aou.vcf.gz | grep -q "gnomAD_AF_non_cancer"; then
        echo -e "ID\nAF\nConsequence\nIMPACT\nSYMBOL\nclinvar_clnsig\nREVEL_score\ngnomAD_AF_non_cancer\ngnomAD_AF_non_cancer_afr\ngnomAD_AF_non_cancer_ami\ngnomAD_AF_non_cancer_amr\ngnomAD_AF_non_cancer_asj\ngnomAD_AF_non_cancer_eas\ngnomAD_AF_non_cancer_fin\ngnomAD_AF_non_cancer_mid\ngnomAD_AF_non_cancer_nfe\ngnomAD_AF_non_cancer_oth\ngnomAD_AF_non_cancer_raw\ngnomAD_AF_non_cancer_sas" > vertical_header.txt
        bcftools query -l aou.vcf.gz >> vertical_header.txt
        # Convert vertical header to horizontal tab-delimited format
        tr '\n' '\t' < vertical_header.txt | sed 's/\t$/\n/' > ~{output_file}
        bcftools +split-vep aou.vcf.gz -f '%ID\t%INFO/AF\t%Consequence\t%IMPACT\t%SYMBOL\t%clinvar_clnsig\t%REVEL_score\t%gnomAD_AF_non_cancer\t%gnomAD_AF_non_cancer_afr\t%gnomAD_AF_non_cancer_ami\t%gnomAD_AF_non_cancer_amr\t%gnomAD_AF_non_cancer_asj\t%gnomAD_AF_non_cancer_eas\t%gnomAD_AF_non_cancer_fin\t%gnomAD_AF_non_cancer_mid\t%gnomAD_AF_non_cancer_nfe\t%gnomAD_AF_non_cancer_oth\t%gnomAD_AF_non_cancer_raw\t%gnomAD_AF_non_cancer_sas\t[%GT\t]' -d > tmp.txt
    else
        echo -e "CHROM\nPOS\nID\nREF\nALT\nAF\nAN\nConsequence\nIMPACT\nSYMBOL\nclinvar_clnsig\nREVEL_score" > vertical_header.txt
        bcftools query -l aou.vcf.gz >> vertical_header.txt
        # Convert vertical header to horizontal tab-delimited format
        tr '\n' '\t' < vertical_header.txt | sed 's/\t$/\n/' > ~{output_file}
        bcftools +split-vep aou.vcf.gz -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AN\t%Consequence\t%IMPACT\t%SYMBOL\t%clinvar_clnsig\%REVEL_score\t[%GT\t]\n' -d > tmp.txt
        
    fi

    # Clear Disk Space
    rm aou.vcf.gz

    # Replace All Genotype Data with counts to cut memory in half
    sed -i \
      -e 's#0/0#0#g' \
      -e 's#0/1#1#g' \
      -e 's#1/1#1#g' \
      -e 's#\./\.#0#g' \
      -e 's#0|1#1#g' \
      -e 's#1|0#1#g' \
      -e 's#1|1#1#g' \
      -e 's#0|0#0#g' \
      tmp.txt

    # Filter to Gene List
    sort -u < tmp.txt > tmp2.txt
    awk 'NR==FNR {cpg[$1]; next} $5 in cpg' ~{gene_list} tmp2.txt >> ~{output_file}
     
    gzip ~{output_file}
  >>>
  output {
    File out1 = "~{output_file}.gz"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    disks: "local-disk ~{default_disk_gb} HDD"
    preemptible: 3
  }
}



task Copy_To_Bucket {
  input {
    Array[File] files
    String google_bucket
    String subdir = "STEP_9_VEP_STATISTICS"
  }
  command <<<
    # Loop through each file and copy it to the subdir in the Google bucket
    for file in ~{sep=" " files}; do
      gsutil cp -m "$file" ~{google_bucket}/~{subdir}/
    done
  >>>
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    read_from_cache: false
  }
}
