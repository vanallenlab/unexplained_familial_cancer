# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_3B_INVESTIGATE {

  input {
    #Array[File] shard_vcfs
    #File metadata_file
    String cancer_type = "kidney"
    Array[String] keep_genes = ["ZNF346"]
    Array[String] exclude_genes = ['TTN']
  }
  String analysis_3_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_3_SAIGE_GENE/results/"
  File saige_results = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_3_SAIGE_GENE/results/" + cancer_type + ".saige.recalibrated.stats.tsv.gz"
  File groupfile_001 = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/ufc_001.groupfile"
  File groupfile_0001 = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/ufc_0001.groupfile"
  ############################
  # Task 1: Filter Genes by P-value
  ############################
  call T1_Filter_Genes {
    input:
      saige_results = saige_results,
      groupfile_001 = groupfile_001,
      groupfile_0001 = groupfile_0001,
      keep_genes = keep_genes,
      exclude_genes = exclude_genes
  }

  ############################
  # Task 3: Find patient carriers in VCFs
  ############################
  scatter (i in range(length(T1_Filter_Genes.out2_genes))) {
    call T2_Find_Carriers {
      input:
        gene = T1_Filter_Genes.out2_genes[i],
        variants = T1_Filter_Genes.out3_variants[i],
        vcfs = read_lines(T1_Filter_Genes.out4_vcfs[i])
    }
  }

  call Tasks.concatenateFiles_noheader {
    input:
      files = T2_Find_Carriers.out1
  }

  ############################
  # Task 4: Merge patient info with metadata
  ############################
  File metadata_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".metadata"
  call T3_Merge_Metadata {
    input:
      carriers_file = concatenateFiles_noheader.out1,
      metadata_file = metadata_file,
      cancer_type = cancer_type
  }
  call Tasks.copy_file_to_storage{
    input:
      text_file = T3_Merge_Metadata.out1,
      output_dir = analysis_3_output_dir
  }

  # output {
  #   File final_patient_report = Merge_Metadata.final_report
  #   File filtered_genes_report = Filter_Genes.filtered_genes
  # }

}

#####################################
# Task Definitions
#####################################

task T1_Filter_Genes {
  input {
    File saige_results
    File groupfile_001
    File groupfile_0001
    Array[String] keep_genes
    Array[String] exclude_genes
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  import os

  keep_genes = "~{sep=' ' keep_genes}".split(' ')
  exclude_genes = "~{sep=' ' exclude_genes}".split(' ')
  def read_groupfile(path):
      data = {}
      with open(path) as f:
          lines = [line.strip() for line in f if line.strip()]
      for i in range(0, len(lines), 2):
          gene = lines[i].split()[0]
          vars_ = lines[i].split()[2:]  # skip "var"
          annos = lines[i + 1].split()[2:]  # skip "anno"
          data[gene] = {"var": vars_, "anno": annos}
      return data

  # Load SAIGE results
  df = pd.read_csv("~{saige_results}", sep='\t', index_col=False)
  df = df[(df['MAC'] <= 50)]
  df_tmp1 = df[(df['recalibrated_p'] <= 0.00001)]
  df_tmp2 = df[(df['#gene'].isin(keep_genes)) & (~df['#gene'].isin(exclude_genes))]
  df_tmp2 = df_tmp2.loc[df_tmp2.groupby('#gene')['recalibrated_p'].idxmin()]
  df = pd.concat([df_tmp1,df_tmp2])
  df = df[df['criteria'].isin(['T1','T3', 'T1;T2', 'T1;T2;T3', 'T1;T2;T3;T4','T5'])]
  # For each Region, find the row with the smallest Pvalue
  df = df.loc[df.groupby('#gene')['recalibrated_p'].idxmin()].reset_index(drop=True)

  print("1: SAIGE results filtered")

  # Load groupfiles
  gf_001 = read_groupfile("~{groupfile_001}")
  gf_0001 = read_groupfile("~{groupfile_0001}")

  print("2: Groupfiles loaded")

  # Collect variants
  gene_variants = {}
  for _, row in df.iterrows():
      gene = row['#gene']
      print(gene)
      maf = row['max_AF']
      group = row['criteria']

      # Choose groupfile
      gf = gf_001 if abs(maf - 0.01) < 1e-6 else gf_0001

      # Skip if gene not in groupfile
      if gene not in gf:
          continue

      variants = gf[gene]['var']
      annotations = gf[gene]['anno']

      annotations_to_keep = group.split(';')
      filtered_vars = [v for v, a in zip(variants, annotations) if a in annotations_to_keep]

      gene_variants[gene] = filtered_vars


  output_file = "variants_summary.txt"
  genes_list_file = "genes_list.txt"

  # Write summary (gene + comma-separated variants)
  with open(output_file, "w") as out, open(genes_list_file, "w") as genes_out:
      for gene in sorted(gene_variants.keys()):
          variants = gene_variants[gene]

          # Write gene name to master list
          genes_out.write(f"{gene}\n")

          # Write summary line to combined file
          joined = ",".join(variants)
          out.write(f"{gene}\t{joined.replace(':', '_')}\n")

          # Write one variant per line to its own file
          per_gene_file = f"{gene}.variants"
          with open(per_gene_file, "w") as gf:
              for v in variants:
                  gf.write(f"{v.replace(':', '_')}\n")
  # Base GCS path
  gcs_base = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_8_FILTER_TO_TP_VARIANTS/sharded_vcfs/"

  # Collect global shards in a set to avoid duplicates
  shards_needed_global = set()

  for gene, variants in gene_variants.items():
      # Collect per-gene shards
      shards_needed_gene = set()

      for v in variants:
          # Example v: chr22:28695232:A:G
          chrom, pos, *_ = v.split(":")
          pos = int(pos)
          mb = pos // 1_000_000  # calculate megabase shard

          # Construct the VCF and index paths
          vcf_path = f"{gcs_base}dfci-ufc-{chrom}-finalrun.{mb}.variant_filtered.vcf.bgz"
          tbi_path = f"{vcf_path}.tbi"

          shards_needed_gene.add(vcf_path)
          shards_needed_global.add(vcf_path)

      # Write per-gene shard file
      per_gene_shards_file = f"{gene}.shards"
      with open(per_gene_shards_file, "w") as gf:
          for shard in sorted(shards_needed_gene):
              gf.write(f"{shard}\n")
              gf.write(f"{shard}.tbi\n")

  CODE
  >>>
  output {
    File out1 = "variants_summary.txt"
    Array[String] out2_genes = read_lines("genes_list.txt")
    Array[File] out3_variants = glob("*.variants")
    Array[String] out4_vcfs = glob("*.shards")
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task Extract_Variants {
  input {
    File filtered_genes
    File groupfile_1pct
    File groupfile_01pct
  }
  command {
    python3 <<CODE
import pandas as pd

filtered_genes = pd.read_csv("filtered_genes.tsv", sep="\t", header=None, names=["Gene","Group","MAF"])
group1 = pd.read_csv("~{groupfile_1pct}", sep="\t")
group01 = pd.read_csv("~{groupfile_01pct}", sep="\t")

# select variants matching filtered genes and thresholds
variants = []
for _, row in filtered_genes.iterrows():
    gene = row['Gene']
    maf = row['MAF']
    group = row['Group']
    df = group1 if maf=="0.01" else group01
    variants.extend(df[(df['Gene']==gene) & (df['Group']==group)]['Variant'].tolist())

with open("filtered_variants.tsv","w") as out:
    for v in variants:
        out.write(v + "\n")
CODE
  }
  output {
    File filtered_variants = "filtered_variants.tsv"
  }
}

task T2_Find_Carriers {
  input {
    String gene
    File variants
    Array[File] vcfs
  }
  command <<<
    set -euo pipefail

    # Filter out .tbi index files
    vcfs_to_use=()
    for vcf in ~{sep=" " vcfs}; do
      if [[ "$vcf" != *.tbi ]]; then
        vcfs_to_use+=("$vcf")
      fi
    done

    # Create output directory
    mkdir -p filtered_vcfs

    # Filter each VCF for variants matching the gene list
    for vcf in "${vcfs_to_use[@]}"; do
      shard_name=$(basename "$vcf" .vcf.bgz)
      bcftools view --include "ID==@~{variants}" "$vcf" -O z -o "filtered_vcfs/${shard_name}.filtered.vcf.gz"
      bcftools index -t "filtered_vcfs/${shard_name}.filtered.vcf.gz"
    done

    # Merge all filtered outputs if there’s more than one
    num_filtered=$(ls filtered_vcfs/*.filtered.vcf.gz | wc -l)
    if [[ "$num_filtered" -gt 1 ]]; then
      echo "Merging $num_filtered filtered VCFs..."
      bcftools concat -a -O z -o "~{gene}.filtered.vcf.gz" filtered_vcfs/*.filtered.vcf.gz
    else
      mv filtered_vcfs/*.filtered.vcf.gz "~{gene}.filtered.vcf.gz"
    fi

    bcftools query -i 'GT="alt"' -f '%ID\t[%SAMPLE\t]\n' "~{gene}.filtered.vcf.gz" > "carrier_file_raw.tsv"

  python3 <<CODE
  # Read the bcftools output where samples are separated by tabs
  with open("carrier_file_raw.tsv") as f:
      lines = [line.strip().split('\t') for line in f if line.strip()]

  # Write exploded output
  with open("carrier_file.tsv", "w") as out:
    out.write("variant\tgene\toriginal_id\n")
    for parts in lines:
      variant_id = parts[0]
      samples = parts[1:]  # remaining fields
      for s in samples:
        out.write(f"{variant_id}\t~{gene}\t{s}\n")
  CODE

  >>>
  output {
    File out1 = "carrier_file.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task T3_Merge_Metadata {
  input {
    File carriers_file
    File metadata_file
    String cancer_type
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  meta = pd.read_csv("~{metadata_file}", sep="\t",index_col=False)
  carriers_df = pd.read_csv("~{carriers_file}", sep='\t',index_col=False)
  # Convert both original_id columns to string
  meta['original_id'] = meta['original_id'].astype(str)
  carriers_df['original_id'] = carriers_df['original_id'].astype(str)
  merged = carriers_df.merge(meta,on='original_id', how="inner")
  columns_to_keep = [
    'variant',
    'gene',
    'original_id',
    'inferred_sex',
    'intake_qc_pop',
    'age',
    'original_dx',
    'maternal_family_dx_counts',
    'paternal_family_dx_counts'
  ]

  merged_subset = merged[columns_to_keep]
  # Define a sorting key for 'original_dx'
  merged_subset['dx_sort_key'] = merged_subset['original_dx'].apply(
      lambda x: (1, '') if str(x).lower() == 'control' else (0, str(x).lower())
  )

  # Sort first by gene, then by dx_sort_key
  merged_subset = merged_subset.sort_values(by=['gene', 'dx_sort_key'])

  # Drop helper column
  merged_subset = merged_subset.drop(columns=['dx_sort_key'])

  merged_subset.to_csv("~{cancer_type}.patient_report", sep="\t", index=False)
  CODE
  >>>
  output {
    File out1 = "~{cancer_type}.patient_report"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible:3
  }
}
