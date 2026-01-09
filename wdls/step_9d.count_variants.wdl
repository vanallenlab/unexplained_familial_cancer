# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# # Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow STEP_9D_COUNT_VARIANTS {
  input {
    String step_8_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_8_FILTER_TO_TP_VARIANTS/sharded_vcfs"
    File subjects_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cohorts/ufc_subjects.aou.list"
    File exclude_samples = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/samples_with_pvs.dec29.list" 
    Array[File] tiered_variants = ["gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier1_001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier2_001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier3_001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier4_001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier5_001.tsv","gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/tier6_001.tsv"]
  }
  
  Int positive_shards = 3000
  call Tasks.gather_positive_vcfs {
    input:
      dir = step_8_output_dir,
      positive_shards = positive_shards
  }

  call Tasks.sort_vcf_list {
    input:
      unsorted_vcf_list = gather_positive_vcfs.vcf_list
  }

  Int negative_shards = 0
  scatter (i in range(length(sort_vcf_list.vcf_arr)-negative_shards)){
    #call T1_Get_Variants as get_tier1{
    #  input:
    #    vcf = sort_vcf_list.vcf_arr[i],
    #    aou_list = subjects_list,
    #    tier_variants = tiered_variants[0]
    #}
    call T1_Get_Variants as get_tier2{
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        aou_list = subjects_list,
        tier_variants = tiered_variants[1]
    }
    #call T1_Get_Variants as get_tier3{
    #  input:
    #    vcf = sort_vcf_list.vcf_arr[i],
    #    aou_list = subjects_list,
    #    tier_variants = tiered_variants[2]
    #}
    #call T1_Get_Variants as get_tier4{
    #  input:
    #    vcf = sort_vcf_list.vcf_arr[i],
    #    aou_list = subjects_list,
    #    tier_variants = tiered_variants[3]
    #}
    #call T1_Get_Variants as get_tier5{
    #  input:
    #    vcf = sort_vcf_list.vcf_arr[i],
    #    aou_list = subjects_list,
    #    tier_variants = tiered_variants[4]
    #}
    #call T1_Get_Variants as get_tier6{
    #  input:
    #    vcf = sort_vcf_list.vcf_arr[i],
    #    aou_list = subjects_list,
    #    tier_variants = tiered_variants[5]
    #}
  }

  #call Tasks.concatenateFiles as concat1{
  #  input:
  #    files = get_tier1.out1
  #}
  call Tasks.concatenateFiles as concat2{
    input:
      files = get_tier2.out1
  }
  #call Tasks.concatenateFiles as concat3{
  #  input:
  #    files = get_tier3.out1
  #}
  #call Tasks.concatenateFiles as concat4{
  #  input:
  #    files = get_tier4.out1
  #}
  #call Tasks.concatenateFiles as concat5{
  #  input:
  #    files = get_tier5.out1
  #}
  #call Tasks.concatenateFiles as concat6{
  #  input:
  #    files = get_tier6.out1
  #}

  # Count Variants
  #call count_variants as count1 {
  #  input:
  #    input_file = concat1.out1,
  #    tier_num = "1"
  #}
  call count_variants as count2 {
    input:
      input_file = concat2.out1,
      tier_num = "2"
  }
  #call count_variants as count3 {
  #  input:
  #    input_file = concat2.out1,
  #    tier_num = "3"
  #}
  #call count_variants as count4 {
  #  input:
  #    input_file = concat2.out1,
  #    tier_num = "4"
  #}
  #call count_variants as count5 {
  #  input:
  #    input_file = concat2.out1,
  #    tier_num = "5"
  #}
  #call count_variants as count6 {
  #  input:
  #    input_file = concat2.out1,
  #    tier_num = "6"
  #}
  
  #Array[File] tier_files = [count1.out1,count2.out1,count3.out1,count4.out1,count5.out1,count6.out1]
 
  #call summarize_variant_counts {
  #  input:
  #    tier_files = tier_files
  #}
  
}
task T1_Get_Variants {
  input {
    File tier_variants
    File aou_list
    File exclude_samples = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/samples_with_pvs.dec29.list"
    File vcf
  }
  command <<<
  cut -f2 ~{tier_variants} > tiered_variants.list
  bcftools view -i 'ID=@tiered_variants.list' ~{vcf} -Oz -o tmp.vcf.gz
  grep -Fvwf ~{exclude_samples} ~{aou_list} > final_samples.list
  bcftools view -S final_samples.list tmp.vcf.gz -Oz -o tmp1.vcf.gz
  bcftools query  -i 'GT="alt"' -f '[%SAMPLE\n]\n' tmp1.vcf.gz | sort | grep -v '^$' > variant_counts.tsv || touch variant_counts.tsv
  >>>
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
  output {
    File out1 = "variant_counts.tsv"
  }
}

task count_variants {
  input {
    File input_file
    String tier_num
  }
  command <<<
  zcat ~{input_file} | sort | uniq -c > variant_counts.tier~{tier_num}.tsv
  >>>
  output{
    File out1 = "variant_counts.tier~{tier_num}.tsv"
  }
  runtime{
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}

task summarize_variant_counts {
  input {
    Array[File] tier_files  # index 0 = tier1 ... index 5 = tier6
  }

  command <<<
  python3 << CODE
  import statistics

  output_file = "variant_count_summary.log"

  with open(output_file, "w") as out:
      out.write("Tier\tMin\tMedian\tMax\n")

      for i, f in enumerate("~{sep=' ' tier_files}".split(' '), start=1):
          counts = []

          with open(f) as infile:
              for line in infile:
                  if not line.strip():
                      continue
                  # uniq -c format: count is first field
                  count = int(line.strip().split()[0])
                  counts.append(count)

          if len(counts) == 0:
              min_v = median_v = max_v = "NA"
          else:
              min_v = min(counts)
              median_v = statistics.median(counts)
              max_v = max(counts)

          out.write(f"Tier{i}\t{min_v}\t{median_v}\t{max_v}\n")

  CODE
  >>>

  output {
    File summary_log = "variant_count_summary.log"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}

