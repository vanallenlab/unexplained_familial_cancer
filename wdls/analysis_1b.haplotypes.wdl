# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_1B_SLIDING_WINDOWS {
  input {
    Int sliding_window_length = 100000
    Array[Int] ROH_length_thresholds = [0]
  }
  File haplotypes_regions_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/genome_windows." + sliding_window_length + "bp.txt"
  String analysis_1b_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/1B_ALL_BY_ALL/"
  File roh_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/ufc_roh_genome.tsv.gz"

  call T1_split_file {
    input:
      file_to_split = haplotypes_regions_file,
      max_lines_per_chunk = 100
  }
  scatter (roh_cutoff in ROH_length_thresholds){
    scatter(i in range(length(T1_split_file.out1) - 0)){
      call T2_data_process {
        input:
          roh_file = roh_file,
          gene_regions = T1_split_file.out1[i],
          roh_cutoff = roh_cutoff,
          output_name = "roh_output"
      }
    }
    call Tasks.concatenateFiles_noheader as concat1{
      input:
        files = select_all(T2_data_process.out1),
        callset_name = "analysis_1b_output.jan9.sw_" + sliding_window_length + ".roh_" + roh_cutoff + "kb.tsv.gz"
    }
    call Tasks.copy_file_to_storage {
      input:
        text_file = concat1.out1,
        output_dir = analysis_1b_output_dir
    }
  }
}

task T1_split_file {
  input {
    File file_to_split
    Int max_lines_per_chunk = 1000
  }

  command <<<
  set -euxo pipefail

    # Make output directory
    mkdir chunks
    head -n 1 ~{file_to_split} > header.tmp

    # Skip the header and split the rest
    tail -n +2 ~{file_to_split} | \
      split -l ~{max_lines_per_chunk} - chunks/part_

    # Prepend the header to each split file
    for f in chunks/part_*; do
      cat header.tmp "$f" > "$f.with_header"
      mv "$f.with_header" "$f"
    done

    # Clean up
    rm header.tmp
  >>>

  output {
    Array[File] out1 = glob("chunks/part_*")
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}

task T2_data_process{
  input {
    File roh_file
    Int roh_cutoff
    File gene_regions # New input variable: the file containing gene regions
    String output_name 
  }
  command<<<
  set -euxo pipefail

  python3 <<CODE
  import pandas as pd

  gene_df = pd.read_csv("~{gene_regions}",sep='\t',index_col=False,header=None,names=['Chr','Start','End','GeneName','Length'])
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
                    names=["Chr", "Start", "End", "Sample","Length"])

  # Filter to ROHs longer than (insert value) bps
  roh_df = roh_df[roh_df["Length"] > ~{roh_cutoff}]

  # Get unique samples
  samples = sorted(roh_df["Sample"].unique())

  # Prepare output file
  with open("~{output_name}.tsv", "w") as f:
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
          #chr_roh = roh_df[roh_df["Chr"].astype(str).str.replace('chr', '') == hap_chr]
          chr_roh = roh_df[roh_df["Chr"].astype(str) == hap_chr]

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
  >>>
  output {
    File? out1 = "~{output_name}.tsv"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
    disk: "local-disks 10 HDD"
    memory: "2GB"
  }
}

