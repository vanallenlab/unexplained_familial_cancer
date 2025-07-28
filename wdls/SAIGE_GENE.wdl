# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# # Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow SAIGE_GENE {
  input {
    String step_9_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/sharded_vcfs"
    File subjects_list
    String cancer_type
    File ld_pruned_vcf = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_11_GENETIC_RELATEDNESS/ufc_1000G_snps.ld_pruned.vcf.gz"
    File autosomal_gene_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/gencode.v47.autosomal.protein_coding.genes.list"
    File sample_data
    File aou_phenotypes = "gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45/dfci-g2c-inputs/phenotypes/dfci-g2c.aou.phenos.tsv.gz"
    File non_aou_phenotypes = "gs://dfci-g2c-inputs/phenotypes/dfci-g2c.non_aou.phenos.tsv.gz"
    String analysis_3_saige_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_3_SAIGE_GENE/results/"
    String single_sex_analysis
  }

  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call gather_vcfs {
    input:
      dir = step_9_output_dir
  }

  call sort_vcf_list {
    input:
      unsorted_vcf_list = gather_vcfs.vcf_list
  }

  Int negative_shards = 0
  
  scatter (i in range(length(sort_vcf_list.vcf_arr)-negative_shards)){
    # Filter VCF to Low Frequency Variants, variants found w/in cohort, and removed INFO
    call process_vcf_part1 {
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        subjects_list = subjects_list,
    }

    # Get the important variants to keep around 
    call make_group_file_part1 as make_group_file_part1_001 {
      input:
        vcf = process_vcf_part1.out1,
        autosomal_gene_file = autosomal_gene_file
    }

    # Get the important variants to keep around
    call make_group_file_part1 as make_group_file_part1_0001 {
      input:
        vcf = process_vcf_part1.out2,
        autosomal_gene_file = autosomal_gene_file
    }
    
    # Clean the vcf to only variants of interest
    call process_vcf_part2 as process_vcf_part2_001 {
      input:
        vcf = process_vcf_part1.out1,
        important_variants = make_group_file_part1_001.important_variants
    }

    # Clean the vcf to only variants of interest
    call process_vcf_part2 as process_vcf_part2_0001 {
      input:
        vcf = process_vcf_part1.out2,
        important_variants = make_group_file_part1_0001.important_variants
    }
  }
  #Array[File] groupfiles = make_group_file_part2.group_file
  
  call concatenateFiles as concatenateFiles_001{
    input:
      files = make_group_file_part1_001.vep_output,
      callset_name = "groupfile.part1.txt"
  }
  call concatenateFiles as concatenateFiles_0001 {
    input:
      files = make_group_file_part1_0001.vep_output,
      callset_name = "groupfile.part1.txt"
  }

  call make_group_file_part2 as make_group_file_part2_001 {
    input:
      split_vep_output = concatenateFiles_001.out,
      cancer_type = cancer_type,
      freq = "001"
  }
 
  call make_group_file_part2 as make_group_file_part2_0001 {
    input:
      split_vep_output = concatenateFiles_0001.out,
      cancer_type = cancer_type,
      freq = "0001"
  }

  call Tasks.copy_file_to_storage as copy0{
    input:
      text_file = make_group_file_part2_001.group_file,
      output_dir = analysis_3_saige_dir
  }
  call Tasks.copy_file_to_storage as copy1{
    input:
      text_file = make_group_file_part2_0001.group_file,
      output_dir = analysis_3_saige_dir
  } 
  call get_covariates {
    input:
      sample_data = sample_data,
      aou_phenotypes = aou_phenotypes,
      non_aou_phenotypes = non_aou_phenotypes,
      subjects_list = subjects_list
  }

  call ConcatVcfs as ConcatVcfs_001 {
    input:
      vcfs = process_vcf_part2_001.output_vcf, 
      vcf_idxs = process_vcf_part2_001.output_vcf_idx,
      callset_name = cancer_type
  } 

  call ConcatVcfs as ConcatVcfs_0001 {
    input:
      vcfs = process_vcf_part2_0001.output_vcf,
      vcf_idxs = process_vcf_part2_0001.output_vcf_idx,
      callset_name = cancer_type
  }

  call plink as plink_LD {
    input:
      vcf = ld_pruned_vcf,
      output_prefix = "ufc_ld"
  }

  call saige_gene_step0 {
    input:
      bed = plink_LD.bed,
      bim = plink_LD.bim,
      fam = plink_LD.fam
  }

  call saige_gene_step1 {
    input:
      bed = plink_LD.bed,
      bim = plink_LD.bim,
      fam = plink_LD.fam,
      single_sex_analysis = single_sex_analysis,
      sample_data = get_covariates.data,
      sparseGRMFile = saige_gene_step0.sparseGRMFile,
      sparseGRMSampleIDFile = saige_gene_step0.sparseGRMSampleIDFile
  }

  call saige_gene_step2 as saige_gene_step2_001 {
    input:
      cancer_type = cancer_type,
      vcf = ConcatVcfs_001.merged_vcf,
      vcf_idx = ConcatVcfs_001.merged_vcf_idx,
      sampleFile = subjects_list,
      rda = saige_gene_step1.rda, 
      group_file = make_group_file_part2_001.group_file,
      varianceRatio = saige_gene_step1.varianceRatio,
      sparseGRMFile = saige_gene_step0.sparseGRMFile,
      sparseGRMSampleIDFile = saige_gene_step0.sparseGRMSampleIDFile
  }

  call saige_gene_step2_beta as saige_gene_step2_beta_001 {
    input:
      cancer_type = cancer_type,
      vcf = ConcatVcfs_001.merged_vcf,
      vcf_idx = ConcatVcfs_001.merged_vcf_idx,
      sampleFile = subjects_list,
      rda = saige_gene_step1.rda,
      group_file = make_group_file_part2_001.group_file,
      varianceRatio = saige_gene_step1.varianceRatio,
      sparseGRMFile = saige_gene_step0.sparseGRMFile,
      sparseGRMSampleIDFile = saige_gene_step0.sparseGRMSampleIDFile
  }

  call saige_gene_step2 as saige_gene_step2_0001 {
    input:
      cancer_type = cancer_type,
      vcf = ConcatVcfs_0001.merged_vcf,
      vcf_idx = ConcatVcfs_0001.merged_vcf_idx,
      sampleFile = subjects_list,
      rda = saige_gene_step1.rda,
      group_file = make_group_file_part2_0001.group_file,
      varianceRatio = saige_gene_step1.varianceRatio,
      sparseGRMFile = saige_gene_step0.sparseGRMFile,
      sparseGRMSampleIDFile = saige_gene_step0.sparseGRMSampleIDFile
  }

  call saige_gene_step2_beta as saige_gene_step2_beta_0001 {
    input:
      cancer_type = cancer_type,
      vcf = ConcatVcfs_0001.merged_vcf,
      vcf_idx = ConcatVcfs_0001.merged_vcf_idx,
      sampleFile = subjects_list,
      rda = saige_gene_step1.rda,
      group_file = make_group_file_part2_0001.group_file,
      varianceRatio = saige_gene_step1.varianceRatio,
      sparseGRMFile = saige_gene_step0.sparseGRMFile,
      sparseGRMSampleIDFile = saige_gene_step0.sparseGRMSampleIDFile
  }

  call sortSAIGE_Output {
    input:
      saige_output_001 = saige_gene_step2_001.out1,
      saige_output_beta_001 = saige_gene_step2_beta_001.out1,
      saige_output_0001 = saige_gene_step2_0001.out1,
      saige_output_beta_0001 = saige_gene_step2_beta_0001.out1,
      output_name = cancer_type
  }

  call Tasks.copy_file_to_storage as copy3{
    input:
      text_file = sortSAIGE_Output.out1,
      output_dir = analysis_3_saige_dir
  }

}

task copy_to_storage {
  input {
    String cancer_type
    File gene_burden_output_file
    File log_file
    File single_var_output_file
    File group_file
  }
  command <<<
  gsutil -m cp ~{gene_burden_output_file}\
          gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/ \
          analysis/~{cancer_type}/~{cancer_type}.saige_results.tsv
  gsutil -m cp ~{single_var_output_file}\
          gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/ \
          analysis/~{cancer_type}/~{cancer_type}.saige_results.single_var.tsv
  gsutil -m cp ~{log_file} \
          gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/ \
          analysis/~{cancer_type}/~{cancer_type}.saige.log
  gsutil -m cp ~{group_file} \
          gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/ \
          analysis/~{cancer_type}/~{cancer_type}.group_file
  >>>
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preepmtible: 3
  }
}

task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String callset_name

    String bcftools_concat_options = ""

    Float mem_gb = 3.5
    Int cpu_cores = 2
    Int? disk_gb

    String bcftools_docker = "vanallenlab/bcftools"
  }

  String out_filename = callset_name + ".vcf.gz"

  Int default_disk_gb = ceil(2.5 * size(vcfs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools concat \
      ~{bcftools_concat_options} \
      --file-list ~{write_lines(vcfs)} \
      -O z \
      -o ~{out_filename} \
      --threads ~{cpu_cores}

    bcftools index --csi ~{out_filename}
  >>>

  output {
    File merged_vcf = "~{out_filename}"
    File merged_vcf_idx = "~{out_filename}.csi"
  }

  runtime {
    docker: bcftools_docker
    memory: mem_gb + " GB"
    cpu: cpu_cores
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
}

task sort_vcf_list {
  input {
    File unsorted_vcf_list
  }
  command <<<
  set -eu -o pipefail
  python3 <<CODE
  import re

  # Input file containing the list of file paths
  file_path = '~{unsorted_vcf_list}'

  # Read and parse file paths
  with open(file_path, 'r') as f:
      paths = f.readlines()

  # Function to extract chromosome and shard number for sorting
  def extract_key(path):
      # Extract chromosome and shard numbers using regex
      match = re.search(r'chr([0-9XY]+)-finalrun\.([0-9]+)', path)
      if match:
          chrom_str, shard_str = match.groups()

          # Convert chromosome to integer, treating 'X' as 23 and 'Y' as 24
          chrom_num = 23 if chrom_str == 'X' else 24 if chrom_str == 'Y' else int(chrom_str)
          shard_num = int(shard_str)
          
          # Return a tuple with chromosome and shard for sorting
          return (chrom_num, shard_num)
      else:
          return (float('inf'), float('inf'))  # Unmatched lines go to the end
  # Sort paths using the extracted keys
  sorted_paths = sorted(paths, key=extract_key)

  # Write to vcf.sorted.list and vcf_idx.sorted.list
  with open('vcf.sorted.list', 'w') as vcf_file, open('vcf_idx.sorted.list', 'w') as vcf_idx_file:
      for path in sorted_paths:
          clean_path = path.strip()
          vcf_file.write(clean_path + '\n')
          vcf_idx_file.write(clean_path + '.tbi\n')
  CODE
  >>>
  output {
    Array[String] vcf_arr = read_lines("vcf.sorted.list")
    Array[String] vcf_idx_arr = read_lines("vcf_idx.sorted.list")
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task gather_vcfs {
  input {
    String dir
  }
  command <<<
  set -euxo pipefail
  gsutil ls ~{dir}/*.vcf.bgz > vcf.list #tmp.txt
  #grep chr17 tmp.txt > vcf.list
  >>>
  output {
    File vcf_list = "vcf.list"
  }
  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    preemptible: 3
  }
}

task process_vcf_part1 {
  input {
    File vcf
    File subjects_list
  }
  Int default_disk_gb = ceil(size(vcf,"GB") * 2) + 10
  command <<<
    set -euxo pipefail

    if ! grep -q '^##INFO=<ID=gnomAD_AF_non_cancer_afr' <(bcftools view -h ~{vcf}); then
        echo "Error: gnomAD annotations not found in VCF header." >&2

        # Add missing FILTER header
        echo '##FILTER=<ID=ExcessHet,Description="Excess heterozygosity filter">' > extra_header.txt
        bcftools annotate -h extra_header.txt -O z -o tmp1.vcf.gz ~{vcf}
        rm ~{vcf}

        bcftools view -S ~{subjects_list} -O z -o tmp2.vcf.gz tmp1.vcf.gz

        bcftools view -h tmp2.vcf.gz -O z -o output.001.vcf.gz
        bcftools view -h tmp2.vcf.gz -O z -o output.0001.vcf.gz
        exit 0
    fi

    # Step 1: Extract relevant fields using bcftools +split-vep
    # Format: Variant ID + gnomAD subpopulation AFs + Cohort AF + Impact
    bcftools +split-vep ~{vcf} \
      -f '%ID\t%gnomAD_AF_non_cancer_afr\t%gnomAD_AF_non_cancer_eas\t%gnomAD_AF_non_cancer_amr\t%gnomAD_AF_non_cancer_fin\t%gnomAD_AF_non_cancer_nfe\t%gnomAD_AF_non_cancer_sas\t%gnomAD_AF_non_cancer_mid\t%gnomAD_AF_non_cancer_ami\t%gnomAD_AF_non_cancer_asj\t%gnomAD_AF_non_cancer_oth\t%AF\t%IMPACT\n' \
      -d | uniq > extracted_variants.txt

    # Sort and deduplicate
    sort -u extracted_variants.txt > sorted_variants.txt

    # Step 2a: Filter variants with:
    # - AF < 1% in major populations (AFR, EAS, AMR, FIN, NFE, SAS)
    # - AF < 10% in minor populations (MID, AMI, ASJ, OTH)
    # - AF < 1% in cohort
    awk -F'\t' '
    BEGIN { OFS = "\t" }
    {
      # Replace missing AFs (".") with 0 and convert to numeric
      for (i = 2; i <= 12; i++) {
        if ($i == ".") $i = 0
        $i += 0
      }

      if (
        ($2  < 0.01) && ($3  < 0.01) && ($4  < 0.01) && ($5  < 0.01) &&
        ($6  < 0.01) && ($7  < 0.01) && ($8  < 0.10) && ($9  < 0.10) &&
        ($10 < 0.10) && ($11 < 0.10) && ($12 < 0.01)
      ) {
        print $1
      }
    }
    ' sorted_variants.txt > variant_ids_AF001.txt

    # Step 2b: Stricter filter:
    # - AF < 0.1% in major populations
    # - AF < 1% in minor populations
    # - AF < 0.1% in cohort
    awk -F'\t' '
    BEGIN { OFS="\t" }
    {
      for (i = 2; i <= 12; i++) {
        $i = ($i == "." ? 0 : $i) + 0
      }

      if (
        $2  < 0.001 && $3  < 0.001 && $4  < 0.001 && $5  < 0.001 &&
        $6  < 0.001 && $7  < 0.001 && $8  < 0.01 && $9  < 0.01 &&
        $10 < 0.01 && $11 < 0.01 && $12 < 0.001
      ) print $1
    }
    ' sorted_variants.txt > variant_ids_AF0001.txt

    # Cleanup intermediate files
    rm extracted_variants.txt sorted_variants.txt

    # Make a list of variant IDs to keep
    #bcftools +split-vep ~{vcf} -f '%ID\t%gnomAD_AF_non_cancer_afr\t%gnomAD_AF_non_cancer_eas\t%gnomAD_AF_non_cancer_amr\t%gnomAD_AF_non_cancer_fin\t%gnomAD_AF_non_cancer_nfe\t%gnomAD_AF_non_cancer_sas\t%gnomAD_AF_non_cancer_mid\t%gnomAD_AF_non_cancer_ami\t%gnomAD_AF_non_cancer_asj\t%gnomAD_AF_non_cancer_oth\t%AF\t%IMPACT\n' -d | uniq > tmp1.txt
    #sort -u tmp1.txt > tmp2.txt

    # Filter to Variants w/ <1% in major gnomAD pops, < 10% in minor gnomAD pops, and < 1% in cohort
    #awk 'BEGIN{OFS="\t"} {$2=($2=="."?0:$2); $3=($3=="."?0:$3); $4=($4=="."?0:$4); $5=($5=="."?0:$5);$6=($6=="."?0:$6); $7=($7=="."?0:$7); $8=($8=="."?0:$8); $9=($9=="."?0:$9); $10=($10=="."?0:$10); $11=($11=="."?0:$11); if($2+0<0.01 && $3+0<0.01 && $4+0<0.01 && $5+0<0.01 && $6+0<0.01 && $7+0<0.01 && $8+0<0.1 && $9+0<0.1 && $10+0<0.1 && $11+0<0.1 && $12+0<0.01) print}' tmp2.txt | cut -f1 > variant_ids.001.txt
    # Filter to Variants w/ <0.1% in major gnomAD pops, < 1% in minor gnomAD pops, and < 0.1% in cohort
    #awk 'BEGIN{OFS="\t"} {$2=($2=="."?0:$2); $3=($3=="."?0:$3); $4=($4=="."?0:$4); $5=($5=="."?0:$5);$6=($6=="."?0:$6); $7=($7=="."?0:$7); $8=($8=="."?0:$8); $9=($9=="."?0:$9); $10=($10=="."?0:$10); $11=($11=="."?0:$11); if($2+0<0.001 && $3+0<0.001 && $4+0<0.001 && $5+0<0.001 && $6+0<0.001 && $7+0<0.001 && $8+0<0.01 && $9+0<0.01 && $10+0<0.01 && $11+0<0.01 && $12+0<0.001) print}' tmp2.txt | cut -f1 > variant_ids.0001.txt
    #rm tmp1.txt tmp2.txt

    bcftools view --include ID==@variant_ids.001.txt ~{vcf} -O z -o tmp1.001.vcf.gz
    bcftools view --include ID==@variant_ids.0001.txt ~{vcf} -O z -o tmp1.0001.vcf.gz
    rm ~{vcf}

    # Add missing FILTER header
    echo '##FILTER=<ID=ExcessHet,Description="Excess heterozygosity filter">' > extra_header.txt
    bcftools annotate -h extra_header.txt -O z -o tmp2.001.vcf.gz tmp1.001.vcf.gz
    bcftools annotate -h extra_header.txt -O z -o tmp2.0001.vcf.gz tmp1.0001.vcf.gz
    rm tmp1.001.vcf.gz tmp1.0001.vcf.gz

    bcftools view -S ~{subjects_list} -O z -o tmp3.001.vcf.gz tmp2.001.vcf.gz
    bcftools view -S ~{subjects_list} -O z -o tmp3.0001.vcf.gz tmp2.0001.vcf.gz
    rm tmp2.001.vcf.gz tmp2.0001.vcf.gz

    bcftools view -i 'AC > 0' tmp3.001.vcf.gz -O z -o output.001.vcf.gz
    bcftools view -i 'AC > 0' tmp3.0001.vcf.gz -O z -o output.0001.vcf.gz
  >>>
  output {
    File out1 = "output.001.vcf.gz"
    File out2 = "output.0001.vcf.gz"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    disks: "local-disk ~{default_disk_gb} HDD"
    preemptible: 3
  }
}

task process_vcf_part2 {
  input {
    File vcf
    File important_variants
  }
  Int default_disk_gb = ceil(size(vcf,"GB") * 2) + 10
  command <<<
    set -euxo pipefail
    bcftools annotate -x INFO -O z -o tmp1.vcf.gz ~{vcf}
    
    if [ $(wc -l < ~{important_variants}) -eq 0 ]; then
      bcftools view -h tmp1.vcf.gz -O z -o  output.vcf.gz
      bcftools index --csi output.vcf.gz
      exit 0
    fi

    bcftools view --include ID==@~{important_variants} tmp1.vcf.gz -O z -o output.vcf.gz

    bcftools index --csi output.vcf.gz
  >>>
  output {
    File output_vcf = "output.vcf.gz"
    File output_vcf_idx = "output.vcf.gz.csi"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    disks: "local-disk ~{default_disk_gb} HDD"
    preemptible: 3
  }
}

task make_group_file_part1 {
  input {
    File vcf
    File autosomal_gene_file
  }
  command <<<
  set -euxo pipefail

  # Get the VEP HIGH and MODERATE IMPACT Variants
  echo "Checkpoint 1"
  bcftools +split-vep ~{vcf} -f '%SYMBOL\t%CHROM:%POS:%REF:%ALT\t%IMPACT\n' -d > tmp0.txt
  grep -v '^\.' tmp0.txt > tmp1.txt || touch tmp1.txt
  grep -E 'HIGH|MODERATE' tmp1.txt > tmp2.txt || touch tmp2.txt
  grep -v '^$' tmp2.txt > tmp3.txt || touch tmp3.txt
  sort -u --compress-program=gzip < tmp3.txt > vep_impact.txt || touch vep_impact.txt

  # Get the ClinVar Important Variants
  #head vep_impact.txt
  #echo "Checkpoint 2"
  #bcftools +split-vep ~{vcf} -f '%SYMBOL\t%CHROM:%POS:%REF:%ALT\t%clinvar_clnsig\t%IMPACT\n' -d > tmp1.txt
  #grep -E 'Pathogenic|Likely_pathogenic|Pathogenic/Likely_pathogenic' tmp1.txt > tmp2.txt || touch tmp2.txt
  #rm tmp1.txt 
  ## Remove instances that are MODERATE or LOW because these variants are pathogenic in other genes most likely
  #grep -Ev 'MODIFIER|LOW' tmp2.txt | cut -f1-3 > tmp3.txt || touch tmp3.txt
  #rm tmp2.txt
  #sort -u tmp3.txt > vep_clinvar.txt || touch vep_clinvar.txt
  #rm tmp3.txt

  # Keep variants if they are greater than 0.5 and are classified as missense variants
  bcftools +split-vep ~{vcf} -f '%SYMBOL\t%CHROM:%POS:%REF:%ALT\t%REVEL_score\t%Consequence\n' -d | \
  grep "missense_variant" > tmp.txt || touch tmp.txt
  cut -f1-3 < tmp.txt | uniq > revel_scores.tsv

  # CODE to parse VEP REVEL Scores
  python3 <<CODE
  def process_revel_score(score_field: str) -> float:
      parts = score_field.split("&")
      numeric_values = [float(p) for p in parts if p != "."]
      return max(numeric_values) if numeric_values else None

  # Example:
  with open("revel_scores.tsv", "r") as infile, open("vep_revel_score.tmp.txt", "w") as outfile:
      for line in infile:
          parts = line.strip().split("\t")
          max_score = process_revel_score(parts[2])
          if max_score is None or max_score < 0.5:
              continue
          elif max_score >= 0.75:
              outfile.write(f"{parts[0]}\t{parts[1]}\tREVEL075\n") 
          elif max_score >= 0.5:
              outfile.write(f"{parts[0]}\t{parts[1]}\tREVEL050\n")
  CODE
  sort -u vep_revel_score.tmp.txt > vep_revel_score.txt

  # Get the synonymous variants
  bcftools +split-vep ~{vcf} -f '%SYMBOL\t%CHROM:%POS:%REF:%ALT\t%Consequence\t%IMPACT\n' -d > tmp0.txt
  grep -v '^\.' tmp0.txt > tmp1.txt || touch tmp1.txt
  grep -E 'synonymous_variant' tmp1.txt > tmp2.txt || touch tmp2.txt
  grep 'LOW' tmp2.txt | cut -f1-3 > tmp3.txt || touch tmp3.txt
  grep -v '^$' tmp3.txt > tmp4.txt || touch tmp4.txt 
  sort -u --compress-program=gzip < tmp4.txt > synonymous.txt || touch synonymous.txt  

  echo "Checkpoint 3"
  cat vep_impact.txt vep_revel_score.txt synonymous.txt > vep_output.tmp.txt
  awk 'NR==FNR {genes[$1]; next} $1 in genes' ~{autosomal_gene_file} vep_output.tmp.txt > vep_output.txt

  echo "Make Important Variant ID List"
  cut -f2 vep_output.txt > tmp0.txt || touch tmp0.txt
  tr ':' '_' < tmp0.txt > tmp1.txt || touch tmp1.txt
  sort -u tmp1.txt > important_variants.list || touch important_variants.list
  >>>
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
    disks: "local-disk 20 HDD"
    memory: "4GB"
  }
  output {
    File vep_output = "vep_output.txt"
    File important_variants = "important_variants.list"
  }
}

task make_group_file_part2 {
  input {
    File split_vep_output
    String cancer_type
    String freq
  }
  command <<<
  set -eu -o pipefail
  python3 <<CODE
  import pandas as pd

  # Load the VEP output (tab-delimited file from +split-vep)
  vep_file = "~{split_vep_output}"
  df = pd.read_csv(vep_file, sep="\t", names=['Gene', 'Variant', 'Consequence'], header=None)
  df = df[df['Gene'] != '.']
  #df = df[~df['Consequence'].str.contains('&',na=False)]
  # Drop Instances of Nonsense Mediated Decay Variants
  df = df.loc[~df['Consequence'].str.contains('nmd_transcript_variant', na=False)]
  df.loc[df['Consequence'].str.contains('synonymous_variant', na=False), 'Consequence'] = 'synonymous_variant'
 
  df = df[~df['Gene'].str.contains('-', na=False)]
  df['Consequence'] = df['Consequence'].replace('Pathogenic/Likely_pathogenic', 'PLP')
  # Removing Pathogenic Clinvar Terms, but if that same variant is annotated otherwise it remains in analysis
  df = df[~df['Consequence'].isin(['Pathogenic', 'Likely_pathogenic', 'PLP'])]
 
  # Print the number of unique values in 'Gene' and 'Variant'
  print(f"Unique values in 'Gene': {df['Gene'].nunique()}")
  print(f"Unique values in 'Variant': {df['Variant'].nunique()}")

  # Define the hierarchy for 'Consequence'
  consequence_hierarchy = {
      'Pathogenic': 1,
      'PLP': 2,
      'Likely_pathogenic': 3,
      'HIGH': 4,
      'REVEL075': 5,
      'REVEL050': 6,
      'MODERATE': 7,
      "synonymous_variant": 8,
      'LOW': 9
  }

  # Add a rank column based on the hierarchy
  df['Consequence_Rank'] = df['Consequence'].map(consequence_hierarchy)

  # Sort by Gene, Variant, and Consequence_Rank
  df_sorted = df.sort_values(by=['Gene', 'Variant', 'Consequence_Rank'])

  # Drop duplicates, keeping the first (highest-ranked Consequence)
  df_filtered = df_sorted.drop_duplicates(subset=['Gene', 'Variant'], keep='first')

  # Drop the rank column (optional)
  df_filtered = df_filtered.drop(columns=['Consequence_Rank'])
  df = df_filtered

  # Extract and format data
  grouped = df.groupby('Gene')

  # Prepare the output
  output = []
  for gene, group in grouped:
      variants = " ".join(group['Variant'].tolist())
      annotations = " ".join(group['Consequence'].tolist())
      output.append(f"{gene} var {variants}")
      output.append(f"{gene} anno {annotations}")

  # Write to the group file
  with open("~{cancer_type}.~{freq}.groupfile", "w") as f:
      f.write("\n".join(output) + "\n")

  CODE
  >>>
  output {
    File group_file = "~{cancer_type}.~{freq}.groupfile"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
    memory: "16GB"
  }
}

task make_cohorts {
  input{
    File sample_data
  }
  command <<<
  set -eu -o pipefail
  python3 <<CODE
  import pandas as pd
  
  df = pd.read_csv("~{sample_data}",sep='\t')
  
  CODE
  >>>
  output {
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
}

task get_covariates {
  input {
    File sample_data
    File aou_phenotypes
    File non_aou_phenotypes
    File subjects_list
  }
  command <<<
  set -eu -o pipefail
  python3 <<CODE
  import pandas as pd
  import numpy as np
  from scipy.stats import zscore

  sample_data = pd.read_csv("~{sample_data}",sep='\t')
  aou_phenotypes = pd.read_csv("~{aou_phenotypes}",sep='\t')
  non_aou_phenotypes = pd.read_csv("~{non_aou_phenotypes}",sep='\t')

  # Normalize sample_data PCS
  covariates = [f"PC{i}" for i in range(1, 5)] 
  sample_data[covariates] = sample_data[covariates].apply(zscore)

  # Concatenate the two DataFrames
  phenotypes = pd.concat([aou_phenotypes, non_aou_phenotypes], ignore_index=True)

  # Identify common columns (excluding the key column for merging)
  common_columns = sample_data.columns.intersection(phenotypes.columns).difference(["Sample", "original_id"])

  # Drop common columns from `phenotypes`
  phenotypes = phenotypes.drop(columns=common_columns)

  # Making Sure that we will merge correctly!
  sample_data['original_id'] = sample_data['original_id'].astype(str).str.strip()
  phenotypes['Sample'] = phenotypes['Sample'].astype(str).str.strip()

  # Merge the DataFrames
  merged_df = pd.merge(sample_data, phenotypes, left_on="original_id", right_on="Sample")
  
  # Add in cancer as a binary trait
  merged_df['cancer_status'] = merged_df['cancer'].apply(lambda x: 1 if x not in ['unknown', 'control'] else 0)

  # List of columns to check for missing values
  columns_to_check = ['age', 'height', 'weight', 'bmi', 'smoking_history', 'birth_year']

  # Create a report on missing values for each column
  missing_report = merged_df[columns_to_check].isnull().sum()

  # Print the detailed report
  for col, missing in missing_report.items():
      print(f"{col}: {missing} missing values")

  # Group by 'cohort' and 'sex' and calculate the median for each group
  grouped_median = merged_df.groupby(['cohort', 'sex_karyotype'])[columns_to_check].median()

  # Function to impute values for a specific column
  def impute_by_group(df, column, grouped_median):
      for idx, row in df.iterrows():
          if pd.isnull(row[column]):
              cohort = row['cohort']
              sex = row['sex_karyotype']
              df.at[idx, column] = grouped_median.loc[(cohort, sex), column]
      return df

  # Impute missing values for each column (excluding 'smoking_history')
  for column in ['age', 'height', 'weight', 'bmi', 'birth_year']:
      merged_df = impute_by_group(merged_df, column, grouped_median)

  # For 'smoking_history', simply impute to 0
  merged_df['smoking_history'].fillna(0, inplace=True)

  # Replace NaN values based on the data type of each column
  for column in merged_df.columns:
      if merged_df[column].dtype == 'object':  # Check if column is of string type
          merged_df[column].fillna('NA', inplace=True)  # Replace NaN with 'NA' in string columns
      else:  # Numeric columns
          merged_df[column].fillna(-1, inplace=True)  # Replace NaN with -1 in numeric columns

  # Read the ~{subjects_list} file into a Python set for efficient lookup
  with open("~{subjects_list}", "r") as file:
      subjects = set(line.strip() for line in file)

  # Filter the DataFrame to include only rows where 'original_id' is in the subjects set
  filtered_df = merged_df[merged_df['original_id'].isin(subjects)]
  filtered_df['saige_cohort'] = np.where(filtered_df['cohort'] == "aou", 1, 0)
  filtered_df['sex'] = (filtered_df['inferred_sex'] == 'male').astype(int)
  filtered_df['age_sex'] = filtered_df['age'] * (filtered_df['inferred_sex'] == 'male').astype(int)

  # Write to file
  filtered_df.to_csv("ufc_covariates.tsv",sep='\t', index=False)
  print("File written to ufc_covariates.tsv")
  CODE
  >>>
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
  }
  output {
    File data = "ufc_covariates.tsv"
  }
}

task saige_gene_step0 {
  input {
    File bed
    File bim
    File fam
  }
  command <<<
  set -eu -o pipefail

  # Ensure all files are in the same directory
  cp ~{bed} ufc.bed
  cp ~{bim} ufc.bim
  cp ~{fam} ufc.fam

  if ! createSparseGRM.R       \
       --plinkFile=./ufc \
       --nThreads=4  \
       --outputPrefix=sparseGRM       \
       --numRandomMarkerforSparseKin=2000      \
       --relatednessCutoff=0.125; then
    echo "Error: Sparse GRM creation failed" >&2
    exit 1
  fi

  echo "Files created by SparseGRM.R"
  ls *mtx*
  >>>
  runtime {
    docker: "wzhou88/saige:1.1.9"
    preemptible: 3
    memory: "8 GB"
  }
  output {
    File sparseGRMFile = "sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
    File sparseGRMSampleIDFile = "sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
  }
}

task saige_gene_step1 {
  input {
    File bed
    File bim
    File fam
    String single_sex_analysis
    File sample_data
    File sparseGRMFile
    File sparseGRMSampleIDFile
  }
  command <<<
  set -eu -o pipefail

  # Ensure all files are in the same directory
  cp ~{bed} ufc.bed
  cp ~{bim} ufc.bim
  cp ~{fam} ufc.fam
 
  echo "Checkpoint 1"

  #run step 1
  #step1_fitNULLGLMM.R \
  #    --plinkFile=./ufc \
  #    --sparseGRMFile=~{sparseGRMFile} \
  #    --sparseGRMSampleIDFile=~{sparseGRMSampleIDFile} \
  #    --useSparseGRMtoFitNULL=TRUE \
  #    --isCateVariance=TRUE \
  #    --phenoFile=~{sample_data} \
  #    --phenoCol=cancer_status \
  #    --traitType=binary \
  #    --covarColList=PC1,PC2,PC3,PC4,age,sex\
  #    --sampleIDColinphenoFile=original_id \
  #    --traitType=binary \
  #    --useSparseGRMforVarRatio=TRUE \
  #    --outputPrefix=./test_run \
  #    --sexCol=sex_karyotype \
  #    --FemaleCode=XO,XX,XXX \
  #    --MaleCode=XY,XYY,XXY
  #
  # Shared parameters
  CMD="step1_fitNULLGLMM.R \
    --plinkFile=./ufc \
    --sparseGRMFile=~{sparseGRMFile} \
    --sparseGRMSampleIDFile=~{sparseGRMSampleIDFile} \
    --useSparseGRMtoFitNULL=TRUE \
    --isCateVariance=TRUE \
    --phenoFile=~{sample_data} \
    --phenoCol=cancer_status \
    --traitType=binary \
    --sampleIDColinphenoFile=original_id \
    --useSparseGRMforVarRatio=TRUE \
    --outputPrefix=./test_run"

  # Add covariates and sex-related arguments if this is NOT a single-sex analysis
  # Examples of when this should be True are Breast, Prostate, Ovarian, etc
  if [ "~{single_sex_analysis}" = "False" ]; then
    CMD="$CMD \
      --covarColList=PC1,PC2,PC3,PC4,age,sex \
      --sexCol=sex_karyotype \
      --FemaleCode=XO,XX,XXX \
      --MaleCode=XY,XYY,XXY"
  else
    CMD="$CMD \
      --covarColList=PC1,PC2,PC3,PC4,age"
  fi

  # Run the final command
  eval $CMD

  >>>
  output {
    File rda = "test_run.rda"
    File varianceRatio = "test_run.varianceRatio.txt"
  }
  runtime {
    docker: "wzhou88/saige:1.3.0"
    preemptible: 3
    memory: "8 GB"
  }
}

task saige_gene_step2 {
  input {
    File vcf
    File vcf_idx
    String cancer_type
    File sampleFile
    File rda
    File group_file
    File varianceRatio
    File sparseGRMFile
    File sparseGRMSampleIDFile
  }

  command <<<

  if ! step2_SPAtests.R \
    --vcfFile=~{vcf} \
    --vcfFileIndex=~{vcf_idx} \
    --vcfField=GT \
    --SAIGEOutputFile=saige_gene_output.~{cancer_type}.txt \
    --LOCO=FALSE \
    --minMAF=0 \
    --minMAC=0.5 \
    --sampleFile=~{sampleFile} \
    --GMMATmodelFile=~{rda} \
    --varianceRatioFile=~{varianceRatio} \
    --sparseGRMFile=~{sparseGRMFile} \
    --sparseGRMSampleIDFile=~{sparseGRMSampleIDFile} \
    --groupFile=~{group_file} \
    --annotation_in_groupTest=HIGH,HIGH:REVEL075,HIGH:REVEL075:REVEL050,HIGH:REVEL075:REVEL050:MODERATE,synonymous_variant \
    --maxMAF_in_groupTest=0.05 \
    2>&1 | tee saige.~{cancer_type}.log ; then
  
    echo "Error: Step2 failed" >&2
    exit 0
  fi

  >>>
  output {
    File out1 = "saige_gene_output.~{cancer_type}.txt"
    File out2 = "saige_gene_output.~{cancer_type}.txt.singleAssoc.txt"
    File out3 = "saige.~{cancer_type}.log" 
  }
  runtime{
    docker:"wzhou88/saige:1.3.0"
    disks: "local-disk 64 HDD"
    preemptible: 3
    memory: "8 GB" 
  }
}

task saige_gene_step2_beta {
  input {
    File vcf
    File vcf_idx
    String cancer_type
    File sampleFile
    File rda
    File group_file
    File varianceRatio
    File sparseGRMFile
    File sparseGRMSampleIDFile
  }

  command <<<

  if ! step2_SPAtests.R \
    --vcfFile=~{vcf} \
    --vcfFileIndex=~{vcf_idx} \
    --vcfField=GT \
    --SAIGEOutputFile=saige_gene_output.~{cancer_type}.txt \
    --LOCO=FALSE \
    --minMAF=0 \
    --minMAC=0.5 \
    --sampleFile=~{sampleFile} \
    --GMMATmodelFile=~{rda} \
    --varianceRatioFile=~{varianceRatio} \
    --sparseGRMFile=~{sparseGRMFile} \
    --sparseGRMSampleIDFile=~{sparseGRMSampleIDFile} \
    --groupFile=~{group_file} \
    --is_no_weight_in_groupTest=TRUE \
    --annotation_in_groupTest=HIGH,HIGH:REVEL075,HIGH:REVEL075:REVEL050,HIGH:REVEL075:REVEL050:MODERATE,synonymous_variant\
    --maxMAF_in_groupTest=0.05 \
    2>&1 | tee saige.~{cancer_type}.log ; then

    echo "Error: Step2 failed" >&2
    exit 0
  fi

  echo "Locations of file"
  ls saige_gene_output*
  >>>
  output {
    File out1 = "saige_gene_output.~{cancer_type}.txt"
    File out2 = "saige_gene_output.~{cancer_type}.txt.singleAssoc.txt"
    File out3 = "saige.~{cancer_type}.log"
  }
  runtime{
    docker:"wzhou88/saige:1.3.0"
    disks: "local-disk 64 HDD"
    preemptible: 3
    memory: "8 GB"
  }
}

task plink {
  input {
    File vcf
    String output_prefix
  }

  Int default_disk_gb = ceil(2 * size(vcf,"GB")) + 20
  command <<<
  set -x

  #if [ "$(bcftools view -H ~{vcf} | wc -l)" -eq 0 ]; then
  #  touch ~{output_prefix}.bed ~{output_prefix}.bim ~{output_prefix}.fam
  #  exit 0
  #fi
  # Creating Plink Files from VCF; --double-id flag makes IID and FID the same value. Needs to be changed eventually
  plink --vcf ~{vcf} --make-bed --out ~{output_prefix} --double-id

  >>>
  output {
    File bed = "~{output_prefix}.bed"
    File bim = "~{output_prefix}.bim"
    File fam = "~{output_prefix}.fam"
  }
  runtime {
    docker: "elixircloud/plink:1.9-20210614"
    disks: "local-disk ~{default_disk_gb} HDD"
    memory: "16 GB"
    preemptible: 3
  }
}

task sortSAIGE_Output {
  input {
    File saige_output_001
    File saige_output_beta_001
    File saige_output_0001
    File saige_output_beta_0001
    String output_name
  }
  String out_filename = output_name + ".saige.tsv"
  command <<<
  set -x
    # 1% AF
    head -n 1 ~{saige_output_001} > header.txt
    tail -n +2 ~{saige_output_001} | grep -v Cauchy > tmp_001.tsv
    tail -n +2 ~{saige_output_beta_001} | grep -v Cauchy > tmp.beta_001.tsv
    
    paste <(cut -f1-5 tmp_001.tsv) <(cut -f6,7 tmp.beta_001.tsv) <(cut -f8- tmp_001.tsv) > tmp.updated_001.tsv
    awk 'BEGIN{OFS="\t"} NR==1 {print; next} { $3 = 0.01; print }' tmp.updated_001.tsv > tmp.updated_again_001.tsv
    sort -k4,4g tmp.updated_again_001.tsv > tmp.final_001.tsv

    # 0.1% AF
    tail -n +2 ~{saige_output_0001} | grep -v Cauchy > tmp_0001.tsv
    tail -n +2 ~{saige_output_beta_0001} | grep -v Cauchy > tmp.beta_0001.tsv

    paste <(cut -f1-5 tmp_0001.tsv) <(cut -f6,7 tmp.beta_0001.tsv) <(cut -f8- tmp_0001.tsv) > tmp.updated_0001.tsv
    awk 'BEGIN{OFS="\t"} NR==1 {print; next} { $3 = 0.001; print }' tmp.updated_0001.tsv > tmp.updated_again_0001.tsv
    sort -k4,4g tmp.updated_again_0001.tsv > tmp.final_0001.tsv

    # Put it all together
    cat header.txt tmp.final_001.tsv tmp.final_0001.tsv > ~{out_filename}
  >>>
  runtime{
    docker:"ubuntu:latest"
    preemptible: 3
  }
  output{
    File out1 = "~{out_filename}"
  }
}

task mergePlinkFiles {
  input {
    Array[File] beds
    Array[File] bims
    Array[File] fams
    Int num_shards
  }
  command <<<
  set -euxo pipefail

  echo "Moving bed files"
  # Move BED files to home directory
  for bed_file in ~{sep=' ' beds}; do
    mv "$bed_file" .
  done
  
  echo "Moving bim files"
  # Move BIM files to home directory
  for bim_file in ~{sep=' ' bims}; do
    mv "$bim_file" .
  done

  echo "Moving fam files"
  # Move FAM files to home directory
  for fam_file in ~{sep=' ' fams}; do
    mv "$fam_file" .
  done

  seq 1 $((~{num_shards} - 1)) > merge_list.txt
  plink --bfile 0 --merge-list merge_list.txt --make-bed --out ufc
  >>>
  output {
    File bed = "ufc.bed"
    File bim = "ufc.bim"
    File fam = "ufc.fam"
  }
  runtime {
    docker: "elixircloud/plink:1.9-20210614"
    preemptible: 3
    memory: "16 GB"
    disks: "local-disk 30 HDD"
  }
}

task concatenateFiles {
  input {
    Array[File] files
    String callset_name
  }
  command <<<
    cat ~{sep=" " files} >> ~{callset_name}.tsv
  >>>
  runtime{
    docker:"ubuntu:latest"
    preemptible: 3
  }
  output{
    File out = "~{callset_name}.tsv"
  }
}
