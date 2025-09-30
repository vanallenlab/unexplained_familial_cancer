# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0

workflow ANALYSIS_3B_INVESTIGATE {

  input {
    Array[File] shard_vcfs
    File groupfile_1pct
    File groupfile_01pct
    File saige_results
    File metadata_file
    String cancer_type
  }

  ############################
  # Task 1: Filter Genes by P-value
  ############################
  call T1_Filter_Genes {
    input:
      saige_results = saige_results
  }

  ############################
  # Task 2: Extract variants from group files
  ############################
  # call Extract_Variants {
  #   input:
  #     filtered_genes = Filter_Genes.filtered_genes
  #     groupfile_1pct = groupfile_1pct
  #     groupfile_01pct = groupfile_01pct
  # }

  ############################
  # Task 3: Find patient carriers in VCFs
  ############################
  # scatter (vcf_file in shard_vcfs) {
  #   call Find_Carriers {
  #     input:
  #       vcf_file = vcf_file
  #       variants_of_interest = Extract_Variants.filtered_variants
  #   }
  # }

  ############################
  # Task 4: Merge patient info with metadata
  ############################
  # call Merge_Metadata {
  #   input:
  #     carriers_files = Find_Carriers.*.carrier_file
  #     metadata_file = metadata_file
  # }

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
  }
  command <<<
  python3 <<CODE
  import pandas as pd

  df = pd.read_csv("~{saige_results}",sep='\t',index_col=False)
  df = df[df['Pvalue'] <= 0.0001]
  df = df[df['Group'].str.isin(['HIGH','HIGH:REVEL050','HIGH:REVEL050:REVEL075'])]
  
  gf_001 = pd.read_csv("~{groupfile_001}",sep='\t',index_col=False)
  gf_0001 = pd.read_csv("~{groupfile_0001}",sep='\t',index_col=False)

  for _, row in df.iterrows():
    gene = row['Region']
    maf = row['max_MAF']
    group = row['Group']
    
    # Decide which groupfile
    gf = gf_001 if maf == 0.01 else gf_0001
    
    # Parse group into individual annotations
    annotations_to_keep = group.split(":")  # e.g., HIGH:REVEL050 -> ["HIGH", "REVEL050"]
    
    # Get the row for this gene in the groupfile
    gf_row = gf[gf['Gene'] == gene]
    if gf_row.empty:
        continue
    
    # Extract variant names and annotations
    variant_cols = [c for c in gf_row.columns if c not in ['Gene','var','anno']]
    
    # Assuming 'var' column has variant list and 'anno' column has annotation list
    variants = gf_row['var'].iloc[0].split()
    annotations = gf_row['anno'].iloc[0].split()
    
    # Filter variants matching any annotation
    filtered_vars = [v for v, a in zip(variants, annotations) if a in annotations_to_keep]
    
    gene_variants[gene] = filtered_vars
    
  CODE
  >>>
  output {
    File filtered_genes = "filtered_genes.tsv"
  }
  runtime {
    docker: "vanallenlab/pydatastack"
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

task Find_Carriers {
  input {
    File vcf_file
    File variants_of_interest
  }
  command {
    zcat ~{vcf_file} | grep -Ff ~{variants_of_interest} > carriers_temp.tsv
    # you could process carriers_temp.tsv to get sample IDs per variant
    mv carriers_temp.tsv carrier_file.tsv
  }
  output {
    File carrier_file = "carrier_file.tsv"
  }
}

task Merge_Metadata {
  input {
    Array[File] carriers_files
    File metadata_file
  }
  command {
    python3 <<CODE
import pandas as pd

meta = pd.read_csv("~{metadata_file}", sep="\t")
carriers = []
for f in [~{carriers_files}]:
    df = pd.read_csv(f, sep="\t", header=None, names=["Variant","Sample"])
    carriers.append(df)
carriers_df = pd.concat(carriers, ignore_index=True)

merged = carriers_df.merge(meta, left_on="Sample", right_on="sample_id", how="left")
merged.to_csv("final_patient_report.tsv", sep="\t", index=False)
CODE
  }
  output {
    File final_report = "final_patient_report.tsv"
  }
}
