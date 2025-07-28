# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_1A_CREATE_HAPLOTYPES {
  input {

    File genetic_map = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/filtered_genetic_map_hg38_withX.txt.gz"
    String analysis_1a_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/HAPLOTYPE_DATA/"
    Array[String] chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
  }
  
  scatter(chr in chromosomes) {
    File chr_roh_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr" + chr + ".txt"
    String output_name = "aou.chr" + chr + ".haplotypes"
    call T1_split_genetic_map {
      input:
        genetic_map = genetic_map,
        chr = chr 
    }
    scatter(i in range(length(T1_split_genetic_map.out1)-0)) {
      call T2_data_process {
        input:
          genetic_map = T1_split_genetic_map.out1[i],
          roh_file = chr_roh_file, 
          chr = chr
      }
    }
    call Tasks.concatenateFiles_noheader {
      input:
        files = T2_data_process.out1,
        callset_name = output_name
    }
    call Tasks.copy_file_to_storage {
      input:
      text_file = concatenateFiles_noheader.out1,
      output_dir = analysis_1a_output_dir    
    } 
  }
}

task T2_data_process{
  input {
    File genetic_map
    File roh_file
    String chr
  }
  command<<<
  set -euxo pipefail
  python3 <<CODE
  import pandas as pd

  # Load genetic map and build haplotypes
  genetic_map = pd.read_csv("~{genetic_map}", delim_whitespace=True, comment="#",skiprows=1, names=["chr", "pos", "rate", "genetic_cM"], dtype={"chr": int, "pos": int})
  genetic_map = genetic_map[genetic_map['chr'] == ~{chr}]

  haplotypes = []
  for i in range(len(genetic_map) - 1):
      row = genetic_map.iloc[i]
      next_row = genetic_map.iloc[i + 1]
      hap_id = f"{int(row['chr'])}_{int(row['pos'])}_{int(next_row['pos'])}"
      haplotypes.append({
          "chr": row["chr"],
          "start": row["pos"],
          "end": next_row["pos"],
          "haplotype": hap_id
      })

  hap_df = pd.DataFrame(haplotypes)

  # Load ROH
  roh_df = pd.read_csv("~{roh_file}", delim_whitespace=True, comment="#", header=None,
                     names=["Type", "Sample", "Chr", "Start", "End", "Length", "Markers", "Quality"])
  roh_df = roh_df[roh_df["Type"] == "RG"]


  # Get unique samples
  samples = sorted(roh_df["Sample"].unique())

  # Prepare output file
  with open("aou.haplotyping_homozygosity_matrix.~{chr}.tsv", "w") as f:
      # Write header
      f.write("Haplotype\t" + "\t".join(map(str, samples)) + "\n")


      for _, hap in hap_df.iterrows():
          hap_chr = str(int(hap["chr"]))
          hap_start = int(hap["start"])
          hap_end = int(hap["end"])
          hap_length = hap_end - hap_start
          hap_id = hap["haplotype"]

          # Filter ROHs on this chromosome
          chr_roh = roh_df[roh_df["Chr"] == "chr" + hap_chr]

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
              percent_homo = round((overlap_bp / hap_length) * 100, 2) if hap_length > 0 else 0.0
              row_vals.append(f"{percent_homo:.2f}")

          # Write this haplotype's row
          f.write(f"{hap_id}\t" + "\t".join(row_vals) + "\n")
  CODE
  >>>
  output {
    File out1 = "aou.haplotyping_homozygosity_matrix.~{chr}.tsv"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}


task T1_split_genetic_map {
  input {
    File genetic_map
    Int chr
    Int max_lines_per_chunk = 1000
  }

  command <<<
  set -euxo pipefail

    # Make output directory
    mkdir chunks
    
    # Normalize Text
    awk 'NR==1 || $1 == "~{chr}" {print $1,$2,$3,$4}' <(zcat ~{genetic_map}) > genetic_map.txt
    rm ~{genetic_map}

    # Extract the header (first line)
    head -n 1 genetic_map.txt > header.tmp

    # Skip the header and split the rest
    tail -n +2 genetic_map.txt | \
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
    docker: "ubuntu:latest"
    preemptible: 3
  }
}
