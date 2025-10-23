version 1.0

workflow SV_ANALYSIS_1_PATHOGENIC_SV {

  input {
    File sv_vcf              # Input VCF with SVs
    File cpg_list            # Text file with list of CPG gene symbols
    File aou_subjects
  }

  # Step 1: Filter rare SVs
  call T1_filter_rare_svs {
    input:
      sv_vcf = sv_vcf,
      aou_subjects = aou_subjects
  }

  # Step 2: Annotate SVs with overlapping genes
  call annotate_svs {
    input:
      rare_vcf = filter_rare_svs.rare_vcf
      gene_list = tsg_list #Tumor Suppressors
  }

  # Step 3: Compute per-patient CPG overlap percentages
  call compute_cpg_overlap {
    input:
      annotated_vcf = annotate_svs.annotated_vcf
      cpg_list = cpg_list
  }

  output {
    File rare_vcf = filter_rare_svs.rare_vcf
    File annotated_vcf = annotate_svs.annotated_vcf
    File patient_cpg_summary = compute_cpg_overlap.patient_cpg_summary
  }
}

#--------------------------------------
task T1_filter_rare_svs {
  input {
    File sv_vcf
    File aou_subjects
  }

  command <<<
  set -euxo pipefail
  bcftools view -S ~{aou_subjects} ~{sv_vcf} -Oz -o tmp.vcf.gz

  # Filter VCF for SVs with gnomAD AF < 0.01 and cohort AF < 0.01
  bcftools view -i 'INFO/gnomad_v3_1_sv_POPMAX_AF<0.01 && AF<0.01' tmp.vcf.gz -o rare_svs.vcf.gz -Oz
  >>>

  output {
    File out1 = "rare_svs.vcf.gz"
  }

  runtime {
    docker: "vanallenlab/bcftools"
    memory: "4G"
    preemptible: 3
  }
}

#--------------------------------------
task T2_assign_svs {
  input {
    File rare_vcf
    File gene_list
  }

  command <<<
    # Annotate SVs with overlapping genes using bedtools intersect
    bcftools query -f '%CHROM:%POS\_%INFO/END:%ALT\t%INFO/PREDICTED_LOF\t[%SAMPLE,]\n' ~{rare_vcf} > sv_temp.tsv

    awk 'NR==FNR {cpg[$1]; next} $2 in cpg' ~{gene_list} sv_tmp.tsv > ufc.pathogenic_svs.tsv
  >>>

  output {
    File annotated_vcf = "ufc.pathogenic_svs.tsv"
  }

  runtime {
    docker: "biocontainers/bedtools:v2.30.0_cv2"
    memory: "4G"
    cpu: 1
  }
}