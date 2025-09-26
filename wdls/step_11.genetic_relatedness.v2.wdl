# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g
# Adapted from Ryan Collins code: 
# https://github.com/vanallenlab/ped_germline_SV/blob/main/gatksv_scripts/get_kinship_and_pcs.py`

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks
workflow STEP_11_GENETIC_RELATEDNESS {
  input {
    File aou_list
    Float min_allele_frequency = 0.01
    Float max_allele_frequency = 0.99
    String step_8_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_8_FILTER_TO_TP_VARIANTS/sharded_vcfs"
    String storage_directory = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
    Array[File] thousand_genomes_vcfs
  }

  Int negative_shards = 0
  call Tasks.gather_chromosome_level_vcfs {
    input:
      dir = step_8_output_dir,
      sex_chromosomes = "FALSE" 
  }

  scatter (i in range(length(thousand_genomes_vcfs))) {
    call GetVariantIDs {
      input:
        vcf = thousand_genomes_vcfs[i]
    }
  }


  call Tasks.concatenateFiles as ThousandGenomesVariants{
    input:
      files = GetVariantIDs.out1
  }

  scatter (chr in range(length(gather_chromosome_level_vcfs.out1) - 0)) {
    call Tasks.sort_vcf_list {
      input:
        unsorted_vcf_list = gather_chromosome_level_vcfs.out1[chr]
    }
  

    scatter( i in range(length(sort_vcf_list.vcf_arr)-negative_shards)){
      call Tasks.Filter_VCF_USING_ID as Filter_VCF_USING_ID_1{
        input:
          vcf = sort_vcf_list.vcf_arr[i],
          variant_list = ThousandGenomesVariants.out1 
      }

      call T4_filter_vcfs {
        input:
          vcf = Filter_VCF_USING_ID_1.out1,
          vcf_idx = Filter_VCF_USING_ID_1.out2,
          min_allele_frequency = min_allele_frequency,
          max_allele_frequency = max_allele_frequency,
          min_QD = 10,
          aou_list = aou_list
      }
    }

    call Tasks.ConcatVcfs as ConcatIntraChromosomes{
      input:
        vcfs = T4_filter_vcfs.out1, 
        vcf_idxs = T4_filter_vcfs.out2,
        callset_name = "chr_specific"
    }

    call T6_prune_variants {
      input:
        vcf = ConcatIntraChromosomes.merged_vcf
    }

    call Tasks.Filter_VCF_USING_ID as Filter_VCF_USING_ID_2{
      input:
        vcf = ConcatIntraChromosomes.merged_vcf,
        variant_list = T6_prune_variants.out1
    }
  }
  
  call Tasks.ConcatVcfs as ConcatVcfsGenome {
    input:
      vcfs = Filter_VCF_USING_ID_2.out1,
      vcf_idxs = Filter_VCF_USING_ID_2.out2,
      callset_name = "ufc_1000G_snps"
  }

  call T9_calculate_plink_pcs {
    input:
      vcf = ConcatVcfsGenome.merged_vcf
  }

}

task T4_filter_vcfs {
  input {
    File aou_list
    File vcf
    File vcf_idx
    Float min_allele_frequency
    Float max_allele_frequency
    Int min_QD
  }
  String filename = basename(vcf, ".vcf.bgz") + ".step11.vcf.bgz"
  command<<<
  set -euxo pipefail

  # 0. Filter the vcf to just AoU samples
  bcftools view -S ~{aou_list} ~{vcf} -Oz -o tmp.vcf.gz
  bcftools index -t tmp.vcf.gz
  bcftools view -i 'AC > 0' tmp.vcf.gz -Oz -o tmp0.vcf.gz
  bcftools index -t tmp0.vcf.gz
  rm ~{vcf} ~{vcf_idx} tmp.vcf.gz*
 
  # 1. Filter on allele frequency, quality, allele number, and pass status 
  bcftools view -i \
    'AF > ~{min_allele_frequency} && AF < ~{max_allele_frequency} && AS_QD > ~{min_QD} && AN >= 100 && FILTER="PASS"' \
    tmp0.vcf.gz -O z -o tmp1.vcf.gz
  bcftools index -t tmp1.vcf.gz
  rm tmp0.vcf.gz*

  # 2. Filter to just SNPs
  bcftools view -v snps -O z -o tmp2.vcf.gz tmp1.vcf.gz
  bcftools index -t tmp2.vcf.gz

  # 3. Filter to just Variants in HWE
  bcftools +fill-tags tmp2.vcf.gz -Oz -o tmp3.vcf.gz -- -t HWE
  bcftools index -t tmp3.vcf.gz
  bcftools view -i 'HWE > 0.05' tmp3.vcf.gz -Oz -o tmp4.vcf.gz
  rm tmp1.vcf.gz* tmp2.vcf.gz* tmp3.vcf.gz*

  # 4. Strip all INFO fields except GT from FORMAT
  bcftools annotate -x INFO,^FORMAT/GT tmp4.vcf.gz -Oz -o ~{filename}
  rm tmp4.vcf.gz*

  # 5. Index the final file
  bcftools index -t ~{filename} 
  >>>
  output {
    File out1 = "~{filename}"
    File out2 = "~{filename}.tbi"
  }
  runtime {
    docker: "vanallenlab/bcftools"
    preemptible: 3
  }
}

task T6_prune_variants {
  input {
    File vcf
  }
  Int default_disk_gb = ceil(size(vcf,"GB")) * 2 + 20
  command <<<
  # perform linkage pruning - i.e. identify prune sites
  plink --vcf ~{vcf} --double-id --allow-extra-chr \
    --indep-pairwise 10000kb 5 0.1 \
    --threads 4 \
    --memory 6000 \
    --out ufc
  >>>
  output {
    File out1 = "ufc.prune.in"
    File out2 = "ufc.prune.out"
  }
  runtime {
    disks: "local-disk ~{default_disk_gb} HDD"
    preemptible: 3
    cpu: "4"
    docker: "elixircloud/plink:1.9-20210614"
  }
}

task T9_calculate_plink_pcs {
  input {
    File vcf
  }
  Int default_disk_gb = ceil(size(vcf,"GB")) * 2 + 20
  command <<<
  plink2 \
    --threads 8 \
    --memory 32000 \
    --vcf ~{vcf} \
    --geno 0.95 \
    --pca \
    --make-king-table \
    --king-table-filter 0.1 \
    --out ufc
  >>>
  output {
    File out1 = "ufc.eigenvec"
    File out2 = "ufc.eigenval"
    File out3 = "ufc.kin0"
  }
  runtime {
    disks: "local-disk ~{default_disk_gb} HDD"
    preemptible: 3
    cpu: "4"
    docker: "vanallenlab/g2c_pipeline"
  }
}

task GetVariantIDs {
  input {
    File vcf
  }

  command <<<
    # Use bcftools query with a filter: only variants where FILTER == "PASS" and AF between 0.01 and 0.99.
    # Note: the condition string uses double equals for FILTER comparison.
    bcftools query -f '%CHROM\_%POS\_%REF\_%ALT\n' -i 'FILTER=="PASS" && AF>=0.01 && AF<=0.99' ~{vcf} > variant_ids.txt
  >>>

  output {
    File out1 = "variant_ids.txt"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    cpu: 1
    memory: "2G"
    preemptible: 3
  }
}
