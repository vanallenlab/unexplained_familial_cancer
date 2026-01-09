# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_1E_ROH_STATS {
  input { 
    String output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/stats/"
    File phenotype_data = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/dfci-ufc.aou.phenos.v2.tsv.gz"
    File sample_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cohorts/ufc_subjects.aou.list"
    File centromeric_regions = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/centromeres.txt.gz"
    Array[File] bcftools_roh_outputs = ["gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr1.txt","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr2.txt","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr3.txt","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr4.txt","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr5.txt","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr6.txt","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr7.txt","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr8.txt","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/roh.aou.chr9.txt"]
  }
  Int negative_shards = 0

  call T1_Centromeric_Regions {
    input:
      centromeric_regions = centromeric_regions
  }

  scatter (roh_file in bcftools_roh_outputs){
    call T2_RoH_Regions {
      input:
        chr_roh_file = roh_file,
        centromeric_regions = T1_Centromeric_Regions.out1
    }
  }
  call Tasks.concatenateFiles {
    input:
      files = T2_RoH_Regions.out1,
      output_name = "ufc_roh_genome"
  }
}

task T2_RoH_Regions {
  input {
    File chr_roh_file
    File centromeric_regions
  }
  command <<<
  grep -v '#' ~{chr_roh_file} > tmp.tsv
  awk '{print $3"\t"$4"\t"$5"\t"$2}' tmp.tsv > tmp2.tsv 
  bedtools sort -i tmp2.tsv > roh_sorted.bed

  awk '{if($3-$2 <= 20000000) print $0}' roh_sorted.bed > roh_filtered.bed

  awk '{OFS="\t"; print $1, $2-2000000, $3+2000000}' ~{centromeric_regions} > centromere_buffered.bed
  bedtools intersect -a roh_filtered.bed -b centromere_buffered.bed -v > roh_no_centromere.bed

  # Merge **per sample**
  cut -f4 roh_no_centromere.bed | sort -u | while read sample; do
      grep -P "\t$sample\$" roh_no_centromere.bed \
      | bedtools sort -i - \
      | bedtools merge -i - -c 4 -o distinct \
      | awk '{print $0"\t"$3-$2}'
  done > roh_merged.bed
  >>>
  output {
    File out1 = "roh_merged.bed"
    File out2 = "roh_no_centromere.bed"
    File out3 = "roh_sorted.bed"
  }
  runtime {
    docker:"vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task T1_Centromeric_Regions {
  input {
    File centromeric_regions
  }

  command <<<
  # Extract relevant columns (chr, start, end)
  zcat ~{centromeric_regions} | cut -f2-4 | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - > merged_regions.bed
  >>>

  output {
    File out1 = "merged_regions.bed"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}
 
