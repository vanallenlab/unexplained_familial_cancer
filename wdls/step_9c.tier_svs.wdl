# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <noah_fields@dfci.harvard.edu>
# # Distributed under the terms of the GNU GPL v2.0

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow tier_sv_variants {
  input {
    String step_9b_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_9_RUN_VEP/"
    Array[File] sv_vcfs = ["gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr1/ConcatVcfs/dfci-ufc.sv.v1.0.chr1.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr2/ConcatVcfs/dfci-ufc.sv.v1.0.chr2.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr3/ConcatVcfs/dfci-ufc.sv.v1.0.chr3.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr4/ConcatVcfs/dfci-ufc.sv.v1.0.chr4.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr5/ConcatVcfs/dfci-ufc.sv.v1.0.chr5.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr6/ConcatVcfs/dfci-ufc.sv.v1.0.chr6.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr7/ConcatVcfs/dfci-ufc.sv.v1.0.chr7.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr8/ConcatVcfs/dfci-ufc.sv.v1.0.chr8.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr9/ConcatVcfs/dfci-ufc.sv.v1.0.chr9.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr10/ConcatVcfs/dfci-ufc.sv.v1.0.chr10.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr11/ConcatVcfs/dfci-ufc.sv.v1.0.chr11.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr12/ConcatVcfs/dfci-ufc.sv.v1.0.chr12.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr13/ConcatVcfs/dfci-ufc.sv.v1.0.chr13.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr14/ConcatVcfs/dfci-ufc.sv.v1.0.chr14.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr15/ConcatVcfs/dfci-ufc.sv.v1.0.chr15.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr16/ConcatVcfs/dfci-ufc.sv.v1.0.chr16.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr17/ConcatVcfs/dfci-ufc.sv.v1.0.chr17.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr18/ConcatVcfs/dfci-ufc.sv.v1.0.chr18.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr19/ConcatVcfs/dfci-ufc.sv.v1.0.chr19.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr20/ConcatVcfs/dfci-ufc.sv.v1.0.chr20.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr21/ConcatVcfs/dfci-ufc.sv.v1.0.chr21.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chr22/ConcatVcfs/dfci-ufc.sv.v1.0.chr22.annotated.vcf.gz","gs://fc-secure-8137c263-494f-45d8-a460-ca2a49b9ebe4/data/AnnotateVcf/chrX/ConcatVcfs/dfci-ufc.sv.v1.0.chrX.annotated.vcf.gz"]
  }

  scatter (i in range(length(sv_vcfs))){
    call filter_sv{
      input:
        vcf = sv_vcfs[i]
    }
  }
  call Tasks.concatenateFiles as concat_001{
    input:
      files = filter_sv.out1,
      output_name = "rare_001"
  }
  call Tasks.concatenateFiles as concat_0001{
    input:
      files = filter_sv.out2,
      output_name = "rare_0001"
  }
  call Tasks.concatenateFiles_noheader{
    input:
      files = filter_sv.out3,
      callset_name = "ufc.rare_svs.tsv"
  }
  call filter_by_sv_type {
    input:
      rare_variants_001 = concat_001.out2,
      rare_variants_0001 = concat_0001.out2
  }
  call Tasks.copy_file_to_storage as copy0{
    input:
      text_file = concatenateFiles_noheader.out1,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy1{
    input:
      text_file = filter_by_sv_type.inversion_001,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy2{
    input:
      text_file = filter_by_sv_type.inversion_0001,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy3{
    input:
      text_file = filter_by_sv_type.duplication_001,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy4{
    input:
      text_file = filter_by_sv_type.duplication_0001,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy5{
    input:
      text_file = filter_by_sv_type.insertion_001,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy6{
    input:
      text_file = filter_by_sv_type.insertion_0001,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy7{
    input:
      text_file = filter_by_sv_type.deletion_001,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy8{
    input:
      text_file = filter_by_sv_type.deletion_0001,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy9{
    input:
      text_file = filter_by_sv_type.cnv_001,
      output_dir = step_9b_output_dir
  }
  call Tasks.copy_file_to_storage as copy10{
    input:
      text_file = filter_by_sv_type.cnv_0001,
      output_dir = step_9b_output_dir
  }
}

task filter_by_sv_type {
  input {
    File rare_variants_001
    File rare_variants_0001
  }
  command <<<
  grep CNV ~{rare_variants_001} > rare_variants.CNV_001.tsv || touch rare_variants.CNV_001.tsv
  grep CNV ~{rare_variants_0001} > rare_variants.CNV_0001.tsv || touch rare_variants.CNV_0001.tsv

  grep DEL ~{rare_variants_001} > rare_variants.DEL_001.tsv || touch rare_variants.DEL_001.tsv
  grep DEL ~{rare_variants_0001} > rare_variants.DEL_0001.tsv || touch rare_variants.DEL_0001.tsv

  grep INS ~{rare_variants_001} > rare_variants.INS_001.tsv || touch rare_variants.INS_001.tsv
  grep INS ~{rare_variants_0001} > rare_variants.INS_0001.tsv || touch rare_variants.INS_0001.tsv

  grep DUP ~{rare_variants_001} > rare_variants.DUP_001.tsv || touch rare_variants.DUP_001.tsv
  grep DUP ~{rare_variants_0001} > rare_variants.DUP_0001.tsv || touch rare_variants.DUP_0001.tsv

  grep INV ~{rare_variants_001} > rare_variants.INV_001.tsv || touch rare_variants.INV_001.tsv
  grep INV ~{rare_variants_0001} > rare_variants.INV_0001.tsv || touch rare_variants.INV_0001.tsv
  >>>
  output {
    File cnv_001 = "rare_variants.CNV_001.tsv"
    File cnv_0001 = "rare_variants.CNV_0001.tsv"
    File deletion_001 = "rare_variants.DEL_001.tsv"
    File deletion_0001 = "rare_variants.DEL_0001.tsv"
    File insertion_001 = "rare_variants.INS_001.tsv"
    File insertion_0001 = "rare_variants.INS_0001.tsv"
    File duplication_001 = "rare_variants.DUP_001.tsv"
    File duplication_0001 = "rare_variants.DUP_0001.tsv"
    File inversion_001 = "rare_variants.INV_001.tsv"
    File inversion_0001 = "rare_variants.INV_0001.tsv"
  }
  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }
}

task filter_sv {
  input {
    File vcf
    File g2c_id_map = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cohorts/ufc_subjects.aou.g2c_id.map" 
    File ufc_subjects_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cohorts/ufc_subjects.aou.list"
  }

  command <<<
    set -euxo pipefail

  bcftools reheader \
    --samples ~{g2c_id_map} \
    -o tmp.vcf.gz \
    ~{vcf}

  tabix -p vcf tmp.vcf.gz

  bcftools view -S ~{ufc_subjects_list} --force-samples tmp.vcf.gz -Oz -o aou.vcf.gz
  bcftools query -i 'AF<0.01 && gnomad_v4.1_sv_AFR_AF < 0.01 && gnomad_v4.1_sv_EUR_AF < 0.01 && gnomad_v4.1_sv_AMR_AF < 0.01 && gnomad_v4.1_sv_EAS_AF < 0.01 && gnomad_v4.1_sv_SAS_AF < 0.01 && gnomad_v4.1_sv_FIN_AF < 0.01 && gnomad_v4.1_sv_ASJ_AF < 0.1 && gnomad_v4.1_sv_AMI_AF < 0.1 && gnomad_v4.1_sv_RMI_AF < 0.1 && gnomad_v4.1_sv_MID_AF < 0.1' -f '%ID\n' aou.vcf.gz > rare_svs.001.tsv

  bcftools query -i 'AF<0.001 && gnomad_v4.1_sv_AFR_AF < 0.001 && gnomad_v4.1_sv_EUR_AF < 0.001 && gnomad_v4.1_sv_AMR_AF < 0.001 && gnomad_v4.1_sv_EAS_AF < 0.001 && gnomad_v4.1_sv_SAS_AF < 0.001 && gnomad_v4.1_sv_FIN_AF < 0.001 && gnomad_v4.1_sv_ASJ_AF< 0.01 && gnomad_v4.1_sv_AMI_AF < 0.01 && gnomad_v4.1_sv_RMI_AF < 0.01 && gnomad_v4.1_sv_MID_AF < 0.01' -f '%ID\n' aou.vcf.gz > rare_svs.0001.tsv

  bcftools view -i 'ID=@rare_svs.001.tsv' aou.vcf.gz -Oz -o tmp_001.vcf.gz
  echo -e 'ID\tAF\tSV_LENGTH\tPREDICTED_BREAKEND_EXONIC\tPREDICTED_COPY_GAIN\tPREDICTED_DUP_PARTIAL\tPREDICTED_INTERGENIC\tPREDICTED_INTRAGENIC_EXON_DUP\tPREDICTED_INTRONIC\tPREDICTED_INV_SPAN\tPREDICTED_LOF\tPREDICTED_MSV_EXON_OVERLAP\tPREDICTED_NEAREST_TSS\tPREDICTED_NONCODING_BREAKPOINT\tPREDICTED_NONCODING_SPAN\tPREDICTED_PARTIAL_DISPERSED_DUP\tPREDICTED_PARTIAL_EXON_DUP\tPREDICTED_PROMOTER\tPREDICTED_TSS_DUP\tPREDICTED_UTR\tSAMPLES\n' >> ufc.structural_variants.txt
  bcftools query -i 'GT="alt"' -f '%ID\t0.01\t%SVLEN\t%PREDICTED_BREAKEND_EXONIC\t%PREDICTED_COPY_GAIN\t%PREDICTED_DUP_PARTIAL\t%PREDICTED_INTERGENIC\t%PREDICTED_INTRAGENIC_EXON_DUP\t%PREDICTED_INTRONIC\t%PREDICTED_INV_SPAN\t%PREDICTED_LOF\t%PREDICTED_MSV_EXON_OVERLAP\t%PREDICTED_NEAREST_TSS\t%PREDICTED_NONCODING_BREAKPOINT\t%PREDICTED_NONCODING_SPAN\t%PREDICTED_PARTIAL_DISPERSED_DUP\t%PREDICTED_PARTIAL_EXON_DUP\t%PREDICTED_PROMOTER\t%PREDICTED_TSS_DUP\t%PREDICTED_UTR\t[%SAMPLE,]\n' tmp_001.vcf.gz >> ufc.structural_variants.txt

  bcftools view -i 'ID=@rare_svs.0001.tsv' aou.vcf.gz -Oz -o tmp_0001.vcf.gz
  bcftools query -i 'GT="alt"' -f '%ID\t0.001\t%SVLEN\t%PREDICTED_BREAKEND_EXONIC\t%PREDICTED_COPY_GAIN\t%PREDICTED_DUP_PARTIAL\t%PREDICTED_INTERGENIC\t%PREDICTED_INTRAGENIC_EXON_DUP\t%PREDICTED_INTRONIC\t%PREDICTED_INV_SPAN\t%PREDICTED_LOF\t%PREDICTED_MSV_EXON_OVERLAP\t%PREDICTED_NEAREST_TSS\t%PREDICTED_NONCODING_BREAKPOINT\t%PREDICTED_NONCODING_SPAN\t%PREDICTED_PARTIAL_DISPERSED_DUP\t%PREDICTED_PARTIAL_EXON_DUP\t%PREDICTED_PROMOTER\t%PREDICTED_TSS_DUP\t%PREDICTED_UTR\t[%SAMPLE,]\n' tmp_0001.vcf.gz >> ufc.structural_variants.txt
  >>>

  output {
    File out1 = "rare_svs.001.tsv"
    File out2 = "rare_svs.0001.tsv"
    File out3 = "ufc.structural_variants.txt"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

