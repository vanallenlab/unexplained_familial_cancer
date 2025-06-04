# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow ANALYSIS_1_ROH {
  input {

    String step_8_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_8_FILTER_TO_TP_VARIANTS/sharded_vcfs"
    File genetic_map = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/genetic_map_hg38_withX.txt.gz"
    String analysis_1_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_1_ROH/"
    String chr
    String filter_on_hwe
    String cohort
    File subjects_list
  }

  call Tasks.gather_a_chromosome_level_vcfs {
    input:
      dir = step_8_output_dir,
      chr_num = chr
  }

  call Tasks.sort_vcf_list {
    input:
      unsorted_vcf_list = gather_a_chromosome_level_vcfs.out1[0]
  }

  # TO DO: Filter to highest quality spots
  scatter (i in range(length(sort_vcf_list.vcf_arr))) {
    call T2_filter_vcfs {
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        filter_on_hwe = filter_on_hwe,
        subjects_list = subjects_list
    }
  }

  call Tasks.ConcatVcfs as ConcatIntraChromosomes{
    input:
      vcfs = T2_filter_vcfs.out1,
      vcf_idxs = T2_filter_vcfs.out2,
      callset_name = chr
  } 

  call T4_runs_of_homozygosity {
    input:
      vcf = ConcatIntraChromosomes.merged_vcf,
      genetic_map = genetic_map,
      chr = chr,
      cohort = cohort
  }
  call Tasks.copy_file_to_storage {
    input:
     text_file = T4_runs_of_homozygosity.out1,
     output_dir = analysis_1_output_dir
      
  } 
}

task T4_runs_of_homozygosity {
  input {
    File vcf
    File genetic_map
    String cohort
    String chr
  }
  Int default_disk_gb = ceil(size(vcf,"GB")) * 2 + 100
  Int default_mem_gb = ceil(size(vcf,"GB")) + 4
  command<<<
  set -euxo pipefail

  # Prep genetic map 
  awk 'NR==1 || $1 == "~{chr}" {print $2,$3,$4}' <(zcat ~{genetic_map}) > genetic_map.txt 

  # 2. Find excess runs of homozygosity
  # https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
  bcftools roh -G 30 --estimate-AF - --skip-indels -m genetic_map.txt ~{vcf} -o roh.chr~{chr}.tmp.txt
  grep -E '^#|^RG' roh.chr~{chr}.tmp.txt > roh.~{cohort}.chr~{chr}.txt
  >>>
  output {
    File out1 = "roh.~{cohort}.chr~{chr}.txt"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
    disks: "local-disk ~{default_disk_gb} HDD"
    memory: "50GB"
  }
}

task T2_filter_vcfs {
  input {
    File vcf
    File subjects_list
    Float min_allele_frequency = 0.01
    Float max_allele_frequency = 0.99
    Int min_QD = 10
    String filter_on_hwe = "False"
  }
  String filename = basename(vcf, ".vcf.bgz") + ".analysis1.vcf.bgz"
  command<<<
  set -euxo pipefail

  # 0. Filter to subjects of interest
  bcftools view -S ~{subjects_list} ~{vcf} -Oz -o tmp.vcf.gz
  bcftools +fill-tags tmp.vcf.gz -Oz -o tmp0.vcf.gz -- -t AN,AC,AF,HWE

  # 1. Filter on allele frequency, quality, allele number
  bcftools view -i \
    'AF > ~{min_allele_frequency} && AF < ~{max_allele_frequency} && AS_QD > ~{min_QD} && AN >= 100' \
    tmp0.vcf.gz -O z -o tmp1.vcf.gz
  bcftools index -t tmp1.vcf.gz
  rm ~{vcf}

  # 2. Filter to just biallelic SNPs
  bcftools view -m 2 -M 2 -v snps -O z -o tmp2.vcf.gz tmp1.vcf.gz
  rm tmp1.vcf.gz*

  # 3. Optional Filter on HWE?
  if [[ "~{filter_on_hwe}" == "True" ]]; then
    bcftools view -i 'INFO/HWE > 0.01' tmp2.vcf.gz -Oz -o tmp3.vcf.gz
  else
    mv tmp2.vcf.gz tmp3.vcf.gz
  fi
  rm tmp2.vcf.gz*

  # 3. Strip all INFO fields except GT from FORMAT
  bcftools annotate -x ^INFO/AF,^FORMAT/GT tmp3.vcf.gz -Oz -o ~{filename}

  # 4. Index the final file
  bcftools index -t ~{filename}
  >>>
  output {
    File out1 = "~{filename}"
    File out2 = "~{filename}.tbi"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}
