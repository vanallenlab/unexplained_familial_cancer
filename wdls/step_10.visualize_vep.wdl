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
    #Array[File] tier_files = ["gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier0_001.tsv"]
    Array[File] tier_files = ["gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier0_001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier0_0001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier1_001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier1_0001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier2_001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier2_0001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier3_001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier3_0001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier4_001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier4_0001.tsv"]
    String step_10_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/FINAL_FILES/"
  }

  call Tasks.gather_chromosome_level_vcfs {
    input:
      dir = step_9_output_dir
  }
  #call Tasks.sort_vcf_list {
  #  input:
  #    unsorted_vcf_list = gather_a_chromosome_level_vcfs.out1[0]
  #}
  call Tasks.concatenateFiles {
    input:
      files = gather_chromosome_level_vcfs.out1,
      output_name = "out"
  }
  call Tasks.sort_vcf_list {
    input:
      unsorted_vcf_list = concatenateFiles.out2
  }

  Int negative_shards = 0
  scatter (tier_file in tier_files){
    scatter (i in range(length(sort_vcf_list.vcf_arr)-negative_shards)){
      call Convert_To_TSV {
        input:
          vcf = sort_vcf_list.vcf_arr[i],
          aou_subjects = aou_subjects,
          gene_list = gene_list,
          tier_variants = tier_file
      }

      call Filter_Vep_TSV{
          input:
            input_tsv = Convert_To_TSV.out1,
            subjects_list = aou_subjects,
            variants_file = tier_file
      }
    }
  
    call merge_variant_counts {
      input:
        tsv_files = Filter_Vep_TSV.out1
    }
  }

  call Tasks.concatenateFiles_noheader {
    input:
      files = merge_variant_counts.out1,
      callset_name = "step_10_output.dec10_2025"
  }

  call Tasks.copy_file_to_storage{
    input:
      text_file = concatenateFiles_noheader.out1,
      output_dir = step_10_output_dir
  }
}



task merge_variant_counts {
  input {
    Array[File] tsv_files
  }

  Int default_mem_gb = ceil(size(tsv_files, "GB")) * 5 + 10
  command <<<
  set -euxo pipefail

  python3 <<CODE
  import pandas as pd

  file_paths = "~{sep=' ' tsv_files}".split()
  combined_df = pd.read_csv(file_paths[0], sep="\t", index_col=0)
  combined_df = combined_df.apply(pd.to_numeric, errors="coerce").fillna(0)

  for path in file_paths[1:]:
      df = pd.read_csv(path, sep="\t", index_col=0)
      df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

      all_rows = combined_df.index.union(df.index)
      all_cols = combined_df.columns.union(df.columns)

      combined_df = combined_df.reindex(index=all_rows, columns=all_cols, fill_value=0)
      df = df.reindex(index=all_rows, columns=all_cols, fill_value=0)

      combined_df += df

  combined_df.to_csv("ufc.cpg.variant_counts.tsv", sep="\t")

  CODE
  #gzip ufc.cpg.variant_counts.tsv
  >>>
  output {
    File out1 = "ufc.cpg.variant_counts.tsv"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "~{default_mem_gb}GB"
    preemptible: 3
  }
}

task Filter_Vep_TSV {
    input {
      File input_tsv
      File subjects_list
      File variants_file
    }
    command <<<
    set -x

    python3 <<CODE
    import pandas as pd
    import numpy as np
    from collections import defaultdict
    import gc
    import re
    import os

    tier, af = re.search(r'tier(\d+)_(\d+)', os.path.basename("~{variants_file}")).groups()
    tier = int(tier)

    # --- Read (gene, variant) pairs from the tier file ---
    valid_pairs = set()
    with open("~{variants_file}", 'r') as f:
        for line in f:
            if line.strip():
                gene, var_id = line.strip().split()[:2]
                valid_pairs.add((gene, var_id))
 
    dtype_map = {
        "ID": "string",
        "AF": "float32",
        "IMPACT": "string",
        "SYMBOL": "string",
    }

    # Read in files
    with open("~{subjects_list}",'r') as f:
        patients = [line.strip() for line in f]

    # Read in the TSV file
    df = pd.read_csv("~{input_tsv}", sep='\t', index_col=False, dtype=dtype_map)
    #df = df[(df['IMPACT'] == "MODERATE")]
    df = df[~df['CLINVAR'].isin(['Pathogenic','Likely_pathogenic','Pathogenic/Likely_pathogenic'])]
    df = df.drop_duplicates(subset=["ID","SYMBOL"])

    # --- Keep only variants in the variant list ---
    df = df[df.apply(lambda r: (r['SYMBOL'], r['ID']) in valid_pairs, axis=1)]

    # Initialize the summary dictionary
    summary = defaultdict(int)

    # Loop through each row
    for index, row in df.iterrows():

        gene = row["SYMBOL"]
        for patient in patients:
          value = pd.to_numeric(row[patient],errors='coerce')
          summary[(gene, f"Tier{tier}_{af}", patient)] += value if pd.notnull(value) else 0

    # Convert summary to a DataFrame
    summary_df = pd.DataFrame.from_dict(summary, orient='index', columns=["variant_count"])

    # Create MultiIndex with names
    summary_df.index = pd.MultiIndex.from_tuples(summary_df.index, names=["gene", "impact", "patient"])

    # Reshape so rows = gene_consequence, columns = patients
    summary_df = summary_df.reset_index()
    summary_df["gene_impact"] = summary_df["gene"] + "_" + f"Tier{tier}_{af}"
    summary_df = summary_df.pivot(index="gene_impact", columns="patient", values="variant_count").fillna(0)

    # Save to TSV
    summary_df.to_csv("variant_counts_by_gene_consequence.tsv", sep="\t")

    CODE

    gzip variant_counts_by_gene_consequence.tsv
    >>>
    output {
      File out1 = "variant_counts_by_gene_consequence.tsv.gz"
    }

    runtime {
      docker: "vanallenlab/g2c_pipeline"
      preemptible: 3
    }
}


task Convert_To_TSV {
  input {
    File vcf
    File aou_subjects
    File gene_list
    File tier_variants

  }

  String output_file = basename(vcf, ".vcf.bgz") + ".tsv"
  Int default_disk_gb = ceil(size(vcf, "GB")) * 10 + 32

  command <<<
  set -euxo pipefail

  awk 'NR==FNR {genes[$1]; next} $1 in genes' ~{gene_list} ~{tier_variants} | cut -f2 > filtered_tier_variants.txt

  echo "### Step 1: Filter to All of Us cohort"
  bcftools view -S ~{aou_subjects} ~{vcf} -Oz -o aou.tmp.vcf.gz
  bcftools +fill-tags aou.tmp.vcf.gz -Oz -o aou.tmp2.vcf.gz -- -t AF
  bcftools view --include ID==@filtered_tier_variants.txt aou.tmp2.vcf.gz -O z -o aou.vcf.gz
  bcftools index -t aou.vcf.gz

  echo "### Step 3: Build Header"
  {
    echo -e "ID"
    echo -e "AF"
    echo -e "SYMBOL"
    echo -e "IMPACT"
    echo -e "CLINVAR"
    bcftools query -l aou.vcf.gz
  } > vertical_header.txt

  tr '\n' '\t' < vertical_header.txt | sed 's/\t$/\n/' > ~{output_file}

  echo "### Step 4: Extract VEP annotations + genotype matrix"
  if bcftools view -h aou.vcf.gz | grep -q "gnomAD_AF_non_cancer"; then
    bcftools +split-vep aou.vcf.gz \
      -f '%ID\t%INFO/AF\t%SYMBOL\t%IMPACT\t%clinvar_CLNSIG\t[%GT\t]' \
      -d > tmp.txt
    rm aou.vcf.gz

    echo "### Step 5: Restrict to Gene List"
    awk 'NR==FNR {cpg[$1]; next} $3 in cpg' ~{gene_list} tmp.txt > tmp1.txt
    mv tmp1.txt tmp.txt

  fi

  echo "### Step 7: Collapse genotypes to allele counts"
  sed -i \
    -e 's#0/0#0#g' \
    -e 's#0/1#1#g' \
    -e 's#1/1#1#g' \
    -e 's#\./\.#0#g' \
    -e 's#0|1#1#g' \
    -e 's#1|0#1#g' \
    -e 's#1|1#1#g' \
    -e 's#0|0#0#g' \
    tmp.txt || touch tmp.txt

  if [ -s tmp.txt ]; then
    cat tmp.txt >> ~{output_file}
  else
    echo "tmp.txt is empty — skipping append"
  fi

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
