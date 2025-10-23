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
    String step_10_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_10_VISUALIZE_VEP/v2/"
  }

  call Tasks.gather_chromosome_level_vcfs {
    input:
      dir = step_9_output_dir
  }

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
  scatter (i in range(length(sort_vcf_list.vcf_arr)-negative_shards)){
    call Convert_To_TSV {
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        aou_subjects = aou_subjects,
        gene_list = gene_list
    }

    call Filter_Vep_TSV{
        input:
          input_tsv = Convert_To_TSV.out1,
          subjects_list = aou_subjects
    }
  }
 
  call merge_variant_counts {
    input:
      tsv_files = Filter_Vep_TSV.out1
  } 

  call Tasks.copy_file_to_storage{
    input:
      text_file = merge_variant_counts.out1,
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
  gzip ufc.cpg.variant_counts.tsv
  >>>
  output {
    File out1 = "ufc.cpg.variant_counts.tsv.gz"
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
    }
    Int default_mem_gb = ceil(size(input_tsv, "GB")) * 10 + 4
    String unzipped_tsv = basename(input_tsv, ".gz")
    String zipped_tsv = basename(input_tsv)
    String output_file_basename = basename(input_tsv, ".tsv.gz") 
    command <<<
    set -x
    mv ~{input_tsv} .
    gunzip ~{zipped_tsv}

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
        "am_class": "string",
        "primateAI_pred": "string",
        "LOFTEE": "string",
        "SpliceAI": "float32"
    }

    # Read in files
    with open("~{subjects_list}",'r') as f:
        patients = [line.strip() for line in f]

    # Read in the TSV file
    df = pd.read_csv("~{unzipped_tsv}", sep='\t', index_col=False, dtype=dtype_map)
    df = df[~df['Consequence'].str.contains('NMD_transcript_variant', case=False)]

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

    def is_rare(row,major_gnomAD_AF,minor_gnomAD_AF):
        low_cols = [
            'gnomAD_AF_non_cancer_afr', 'gnomAD_AF_non_cancer_amr', 'gnomAD_AF_non_cancer_eas',
            'gnomAD_AF_non_cancer_fin', 'gnomAD_AF_non_cancer_nfe', 'gnomAD_AF_non_cancer_sas', 'AF'
        ]
        high_cols = [
            'gnomAD_AF_non_cancer_asj', 'gnomAD_AF_non_cancer_mid',
            'gnomAD_AF_non_cancer_ami', 'gnomAD_AF_non_cancer_oth'
        ]

        # Check low-frequency columns (< 0.01)
        for col in low_cols:
            if col in row and pd.notnull(row[col]):
                try:
                    af = float(str(row[col]).split(',')[0])
                    if af >= major_gnomAD_AF:
                        return False
                except ValueError:
                    continue

        # Check high-frequency columns (< 0.1)
        for col in high_cols:
            if col in row and pd.notnull(row[col]):
                try:
                    af = float(str(row[col]).split(',')[0])
                    if af >= minor_gnomAD_AF:
                        return False
                except ValueError:
                    continue

        return True


    # Initialize the summary dictionary
    summary = defaultdict(int)

    # Loop through each row
    for index, row in df.iterrows():
       
        if not is_rare(row,major_gnomAD_AF=0.01,minor_gnomAD_AF=0.1):
           continue

        gene = row["SYMBOL"]
        consequence_col = row["Consequence"]
        impact = row['IMPACT']
        clinvar_uncertain = row['clinvar_clnsig']

        if clinvar_uncertain == "Uncertain_significance":
          for patient in patients:
            value = pd.to_numeric(row[patient],errors='coerce')
            summary[(gene, "VUS_001", patient)] += value

        if impact == "HIGH":
          for patient in patients:
            value = pd.to_numeric(row[patient],errors='coerce')
            summary[(gene, "HIGH_001", patient)] += value
          continue

        consequences = consequence_col.split('&')

        for consequence in consequences:        
          if consequence == "missense_variant" and impact == "MODERATE":
            for patient in patients:
              value = pd.to_numeric(row[patient], errors='coerce')
              if row['REVEL_score_numeric'] >= 0.5:
                summary[(gene, "REVEL050_001", patient)] += value
                if row['REVEL_score_numeric'] >= 0.75:
                    summary[(gene, "REVEL075_001", patient)] += value

        if not is_rare(row,major_gnomAD_AF=0.001,minor_gnomAD_AF=0.01):
           continue

        if clinvar_uncertain == "Uncertain_significance":
          for patient in patients:
            value = pd.to_numeric(row[patient],errors='coerce')
            summary[(gene, "VUS_0001", patient)] += value

        if impact == "HIGH":
          for patient in patients:
            value = pd.to_numeric(row[patient],errors='coerce')
            summary[(gene, "HIGH_0001", patient)] += value
          continue

        consequences = consequence_col.split('&')

        for consequence in consequences:
          if consequence == "missense_variant" and impact == "MODERATE":
            for patient in patients:
              value = pd.to_numeric(row[patient], errors='coerce')
              if row['REVEL_score_numeric'] >= 0.5:
                summary[(gene, "REVEL050_0001", patient)] += value
                if row['REVEL_score_numeric'] >= 0.75:
                    summary[(gene, "REVEL075_0001", patient)] += value 

    # Convert summary to a DataFrame
    summary_df = pd.DataFrame.from_dict(summary, orient='index', columns=["variant_count"])

    # Create MultiIndex with names
    summary_df.index = pd.MultiIndex.from_tuples(summary_df.index, names=["gene", "impact", "patient"])

    # Reshape so rows = gene_consequence, columns = patients
    summary_df = summary_df.reset_index()
    summary_df["gene_impact"] = summary_df["gene"] + "_" + summary_df["impact"]
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
      memory: "~{default_mem_gb}GB"
      preemptible: 3
    }
}


task Convert_To_TSV {
  input {
    File vcf
    File aou_subjects
    File gene_list

    # Default AF cutoff for major gnomAD pops + cohort
    Float default_af = 0.01
  }

  String output_file = basename(vcf, ".vcf.bgz") + ".tsv"
  Int default_disk_gb = ceil(size(vcf, "GB")) * 10 + 32

  command <<<
  set -euxo pipefail

  echo "### Step 1: Filter to All of Us cohort"
  bcftools view -S ~{aou_subjects} ~{vcf} -Oz -o aou.tmp.vcf.gz
  bcftools view -i 'AC > 0' aou.tmp.vcf.gz -Oz -o aou.tmp2.vcf.gz
  bcftools +fill-tags aou.tmp2.vcf.gz -Oz -o aou.tmp3.vcf.gz -- -t AF
  bcftools index -t aou.tmp3.vcf.gz

  echo "### Step 2: Filter to Rare Variants in cohort (< ~{default_af})"
  bcftools view -i "AF <= ~{default_af}" aou.tmp3.vcf.gz -Oz -o aou.vcf.gz
  bcftools index -t aou.vcf.gz
  rm ~{vcf} aou.tmp.vcf.gz aou.tmp2.vcf.gz aou.tmp3.vcf.gz*

  echo "### Step 3: Build Header"
  {
    echo -e "ID"
    echo -e "AF"
    echo -e "Consequence"
    echo -e "IMPACT"
    echo -e "SYMBOL"
    echo -e "clinvar_clnsig"
    echo -e "REVEL_score"
    echo -e "gnomAD_AF_non_cancer"
    echo -e "gnomAD_AF_non_cancer_afr"
    echo -e "gnomAD_AF_non_cancer_ami"
    echo -e "gnomAD_AF_non_cancer_amr"
    echo -e "gnomAD_AF_non_cancer_asj"
    echo -e "gnomAD_AF_non_cancer_eas"
    echo -e "gnomAD_AF_non_cancer_fin"
    echo -e "gnomAD_AF_non_cancer_mid"
    echo -e "gnomAD_AF_non_cancer_nfe"
    echo -e "gnomAD_AF_non_cancer_oth"
    echo -e "gnomAD_AF_non_cancer_raw"
    echo -e "gnomAD_AF_non_cancer_sas"
    bcftools query -l aou.vcf.gz
  } > vertical_header.txt

  tr '\n' '\t' < vertical_header.txt | sed 's/\t$/\n/' > ~{output_file}

  echo "### Step 4: Extract VEP annotations + genotype matrix"
  if bcftools view -h aou.vcf.gz | grep -q "gnomAD_AF_non_cancer"; then
    bcftools +split-vep aou.vcf.gz \
      -f '%ID\t%INFO/AF\t%Consequence\t%IMPACT\t%SYMBOL\t%clinvar_clnsig\t%REVEL_score\t%gnomAD_AF_non_cancer\t%gnomAD_AF_non_cancer_afr\t%gnomAD_AF_non_cancer_ami\t%gnomAD_AF_non_cancer_amr\t%gnomAD_AF_non_cancer_asj\t%gnomAD_AF_non_cancer_eas\t%gnomAD_AF_non_cancer_fin\t%gnomAD_AF_non_cancer_mid\t%gnomAD_AF_non_cancer_nfe\t%gnomAD_AF_non_cancer_oth\t%gnomAD_AF_non_cancer_raw\t%gnomAD_AF_non_cancer_sas\t[%GT\t]' \
      -d > tmp.txt
    rm aou.vcf.gz

    echo "### Step 5: Restrict to Gene List"
    awk 'NR==FNR {cpg[$1]; next} $5 in cpg' ~{gene_list} tmp.txt > tmp1.txt
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

  echo "### Step 8: Keep HIGH/MODERATE impact"
  if grep -q -E 'HIGH|MODERATE|Uncertain_significance' tmp.txt; then
    grep -E 'HIGH|MODERATE|Uncertain_significance' tmp.txt | grep -Ev 'Benign|Likely_benign|Benign/Likely_benign' | sort -u >> ~{output_file}
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
