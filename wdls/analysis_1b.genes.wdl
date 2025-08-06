# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_1B_GENES {
  input {

    String analysis_1_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/"
    String analysis_1b_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/1B_GENES/"
    Array[String] chroms = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
    File gene_regions_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis_1b.bed.gz"
  }

  scatter(chr in chroms){
   File chr_roh_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr" + chr + ".txt"
    String output_name = "aou.chr" + chr + ".genes.tsv"
    call T2_data_process {
      input:
        roh_file = chr_roh_file,
        chr = chr,
        gene_regions = gene_regions_file,
        output_name = output_name
    }

    call Tasks.copy_file_to_storage {
      input:
        text_file = T2_data_process.out1,
        output_dir = analysis_1b_output_dir
    }
  }
}

task T2_data_process{
  input {
    File roh_file
    String chr
    File gene_regions # New input variable: the file containing gene regions
    String output_name 
  }
  command<<<
  set -euxo pipefail

  python3 <<CODE
  import pandas as pd

  # Load gene regions from the input file
  # Assuming the format: chr<TAB>start<TAB>end<TAB>gene_name
  # Example: chr1    65419    71585    OR4F5
  gene_df = pd.read_csv("~{gene_regions}", delim_whitespace=True, header=None,
                        names=["Chr", "Start", "End", "GeneName"])

  # Filter gene regions for the current chromosome
  # Ensure the 'chr' input matches the format in gene_df (e.g., "chr1" vs "1")
  gene_df['Chr'] = gene_df['Chr'].astype(str).str.replace('chr', '')
  gene_df = gene_df[gene_df['Chr'] == "~{chr}"]

  # Create 'haplotypes' DataFrame directly from gene_df
  # The 'haplotype' ID can be a combination of chr_start_end or gene name
  haplotypes = []
  for index, row in gene_df.iterrows():
      hap_id = f"{row['Chr']}_{int(row['Start'])}_{int(row['End'])}_{row['GeneName']}"
      haplotypes.append({
          "chr": row["Chr"],
          "start": row["Start"],
          "end": row["End"],
          "haplotype": hap_id,
          "gene_name": row["GeneName"] # Keep gene name for output if desired
      })
  hap_df = pd.DataFrame(haplotypes)

  # Load ROH
  roh_df = pd.read_csv("~{roh_file}", delim_whitespace=True, comment="#", header=None,
                    names=["Type", "Sample", "Chr", "Start", "End", "Length", "Markers", "Quality"])
  roh_df = roh_df[roh_df["Type"] == "RG"]

  # Filter to ROHs longer than 100,000 bps
  roh_df = roh_df[roh_df["Length"] > 100_000]

  # Get unique samples
  samples = sorted(roh_df["Sample"].unique())

  # Prepare output file
  with open("~{output_name}", "w") as f:
      # Write header
      f.write("Genomic_Region\tGeneName\t" + "\t".join(map(str, samples)) + "\n")


      for _, hap in hap_df.iterrows():
          hap_chr = str(hap["chr"])
          hap_start = int(hap["start"])
          hap_end = int(hap["end"])
          hap_length = hap_end - hap_start
          hap_id = hap["haplotype"]
          gene_name = hap["gene_name"]

          # Filter ROHs on this chromosome
          # Ensure ROH 'Chr' column matches 'hap_chr' format (e.g., "chr1" vs "1")
          chr_roh = roh_df[roh_df["Chr"].astype(str).str.replace('chr', '') == hap_chr]

          row_vals = []

          for sample in samples:
              sample_roh = chr_roh[chr_roh["Sample"] == sample]

              overlapping = sample_roh[
                  (sample_roh["End"] > hap_start) &
                  (sample_roh["Start"] < hap_end)
              ]

              # Merge overlapping intervals
              intervals = []
              for _, roh in overlapping.iterrows():
                  start = max(hap_start, roh["Start"])
                  end = min(hap_end, roh["End"])

                  if start < end:
                      intervals.append((start, end))

              intervals.sort()
              merged = []
              for interval in intervals:
                  if not merged or interval[0] > merged[-1][1]:
                      merged.append(list(interval))
                  else:
                      merged[-1][1] = max(merged[-1][1], interval[1])

              overlap_bp = sum(end - start for start, end in merged)
              percent_homo = round((overlap_bp / hap_length), 2) if hap_length > 0 else 0.0
              row_vals.append(f"{percent_homo:.2f}")

          # Write this gene region's row
          f.write(f"{hap_id}\t{gene_name}\t" + "\t".join(row_vals) + "\n")
  CODE
  gzip ~{output_name}
  >>>
  output {
    File out1 = "~{output_name}.gz"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
    disk: "local-disks 20 HDD"
    memory: "8GB"
  }
}

