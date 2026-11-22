# The Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow BRAT1_FORMAL_ANALYSIS {
  input {
    File vcf = "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/BRAT1_region.vcf.bgz"
    File vcf_idx = "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/BRAT1_region.vcf.bgz.tbi"
    String callset_name = "brat1_haplotype_ld_log_reg"
    String ID = "chr7_2555487_C_A"
  }

  call T0_extract_BRAT1_gt {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx
  }
  
  call T1_create_bed_files

  scatter (bed_file in T1_create_bed_files.out1) {
    call T2_shard_vcf {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        bed_file = bed_file
    }
    call T3_plink_log_reg {
      input:
        vcf = T2_shard_vcf.out1,
        brat1_alt = T0_extract_BRAT1_gt.alt,
        brat1_hom_ref = T0_extract_BRAT1_gt.hom_ref
    }
  }

  call Tasks.concatenateFiles_noheader as concat_logreg{
    input:
      files = T3_plink_log_reg.out1,
      callset_name = callset_name
  }
  call T5_filter_vcf {
    input:
      vcf = vcf,
      association_results = concat_logreg.out2
  }
  call T6_beagle_phasing {
    input:
      vcf = T5_filter_vcf.out1
  }
  #call T7_find_haplotype{
  #  input:
  #    vcf = T6_beagle_phasing.phased_vcf,
  #    brat1_het = T0_extract_BRAT1_gt.het,
  #}
  call T7_chunk_genotypes {
    input:
      vcf = T6_beagle_phasing.phased_vcf
  }
  scatter (tsv in T7_chunk_genotypes.out1){
    call T8_count_snps {
      input:
        tsv = tsv,
        association_results = concat_logreg.out2
    }
  }
  call T9_find_haplotype{
    input:
      vcf = T6_beagle_phasing.phased_vcf,
      brat1_het = T0_extract_BRAT1_gt.alt
  }
  call T10_paste_output {
    input:
      haplotype_scores = T8_count_snps.out1,
      phased_samples = T9_find_haplotype.out1
  }
}

task T0_extract_BRAT1_gt {
  input {
    File vcf
    File vcf_idx
  }

  command <<<
    set -euo pipefail

    # Extract the BRAT1 region and isolate the variant of interest
    bcftools view \
      --threads 4 \
      -r chr7:2555486-2555488 \
      -i 'ID="chr7_2555487_C_A"' \
      -Oz \
      -o brat1.vcf.gz \
      ~{vcf}

    # Index the resulting single-variant VCF
    bcftools index --threads 4 -t brat1.vcf.gz

    # Extract sample IDs and genotypes
    bcftools query \
      -f '[%SAMPLE\t%GT\n]' \
      brat1.vcf.gz > brat1_genotypes.tsv

    # Classify genotypes
    awk -F'\t' '
      BEGIN {
        OFS = "\t"
      }
      {
        sample = $1
        gt = $2

        if (gt == "0/0" || gt == "0|0") {
          print sample > "brat1_hom_ref.list"
        }
        else if (gt == "0/1" || gt == "1/0" ||
                 gt == "0|1" || gt == "1|0") {
          print sample > "brat1_het.list"
          print sample > "brat1_alt.list"
        }
        else if (gt == "1/1" || gt == "1|1") {
          print sample > "brat1_hom_alt.list"
          print sample > "brat1_alt.list"
        }
        else {
          print sample > "brat1_other_gt.list"
        }
      }
    ' brat1_genotypes.tsv

  cat brat1_hom_ref.list brat1_alt.list > brat1_cohort.list
  >>>

  output {

    File brat1_cohort = "brat1_cohort.list"
    File out_vcf = "brat1.vcf.gz"
    File out_index = "brat1.vcf.gz.tbi"

    File hom_ref = "brat1_hom_ref.list"
    File het = "brat1_het.list"
    File hom_alt = "brat1_hom_alt.list"
    File alt = "brat1_alt.list"
    File other_gt = "brat1_other_gt.list"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bcftools"
    cpu: 4
    memory: "12 GB"
  }
}

task T1_create_bed_files {
  input {
    Int start = 1555487
    Int end = 3555487
    Int shards = 100
  }

  command <<<
  set -euxo pipefail

    fixed_start=2555486
    fixed_end=2555488

    for ((i=0; i<~{shards}; i++)); do
      shard_start=$((~{start} + i * (~{end} - ~{start}) / ~{shards}))
      shard_end=$((~{start} + (i + 1) * (~{end} - ~{start}) / ~{shards}))

      printf -v shard_num "%03d" "$i"

      printf "chr7\t%s\t%s\nchr7\t%s\t%s\n" \
        "$shard_start" "$shard_end" \
        "$fixed_start" "$fixed_end" \
        > "chr7_shard_${shard_num}.bed"
    done
  >>>
  output{
    Array[File] out1 = glob("*.bed")
  }
  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/ubuntu"
  }
}

task T2_shard_vcf {
  input {
    File vcf
    File vcf_idx
    File bed_file
  }

  command <<<
    set -euxo pipefail

    bcftools view \
      --threads 4 \
      -R ~{bed_file} \
      -Oz \
      -o shard.vcf.gz \
      ~{vcf}

    bcftools index \
      --threads 4 \
      -t shard.vcf.gz
  >>>

  output {
    File out1 = "shard.vcf.gz"
    File out2 = "shard.vcf.gz.tbi"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bcftools"
    cpu: 4
    preemptible:3
    memory: "16 GB"
  }
}

task T3_plink_association {
  input {
    File vcf
    File brat1_alt
    File brat1_hom_ref
  }

  command <<<
    set -euo pipefail

    # Convert VCF to PLINK binary format
    plink \
      --vcf ~{vcf} \
      --make-bed \
      --out brat1

    # Create phenotype file
    awk '{print $1, $1, 2}' ~{brat1_alt} > brat1.pheno
    awk '{print $1, $1, 1}' ~{brat1_hom_ref} >> brat1.pheno

    # Run chi-square association test
    plink \
      --bfile brat1 \
      --pheno brat1.pheno \
      --allow-no-sex \
      --assoc \
      --out tmp

  awk 'BEGIN {OFS="\t"} {print $2,$5,$6,$8,$9,$10}' tmp.assoc | grep -v NA > brat1_assoc.assoc

  >>>

  output {
    File out1 = "brat1_assoc.assoc"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bioinformatics_tools"
    cpu: 4
    memory: "12 GB"
    preemptible: 3
  }
}

task T3_plink_log_reg {
  input {
    File vcf
    File pcs = "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/AoU_CDRv9_pcs.tsv.gz"
    File brat1_alt
    File brat1_hom_ref
  }

  command <<<
    set -euxo pipefail

    # Prepare PCs
    zcat ~{pcs} | awk 'BEGIN {OFS="\t"} {
        printf "%s", $1
        for (i=3; i<=NF; i++) printf "\t%s", $i
        print ""
    }' > pcs.covar

    sed -i '1s/.*/#IID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10/' pcs.covar

    # Convert VCF to PLINK 2 format (preserves exact VCF REF and ALT)
    plink2 \
      --vcf ~{vcf} \
      --make-pgen \
      --out brat1

    # Create phenotype file (header added for PLINK 2 compatibility)
    echo -e "#IID\tPHENO" > brat1.pheno
    awk '{print $1 "\t2"}' ~{brat1_alt} >> brat1.pheno
    awk '{print $1 "\t1"}' ~{brat1_hom_ref} >> brat1.pheno

    # Run ancestry-adjusted log reg association test
    plink2 \
      --pfile brat1 \
      --pheno brat1.pheno \
      --pheno-name PHENO \
      --allow-no-sex \
      --covar pcs.covar \
      --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
      --glm dominant firth-fallback hide-covar omit-ref \
      --out brat1_ancestry_adjusted

    plink2 \
      --pfile brat1 \
      --pheno brat1.pheno \
      --pheno-name PHENO \
      --keep-if PHENO == 2 \
      --freq cols=chrom,pos,ref,alt1,alt1freq,nobs \
      --out brat1_cases

    plink2 \
      --pfile brat1 \
      --pheno brat1.pheno \
      --pheno-name PHENO \
      --keep-if PHENO == 1 \
      --freq cols=chrom,pos,ref,alt1,alt1freq,nobs \
      --out brat1_controls

  awk 'BEGIN {OFS="\t"} {print $3,$6,$7}' brat1_cases.afreq > cases.afreq
  awk 'BEGIN {OFS="\t"} {print $3,$6,$7}' brat1_controls.afreq > controls.afreq

  cut -f2,3,9,10,12,13,14,15,16,17 brat1_ancestry_adjusted.PHENO.glm.logistic.hybrid > brat1.logistic

  python3 <<CODE
  import pandas as pd
  df_cases_freq = pd.read_csv("cases.afreq", sep=r"\s+", index_col=False)
  df_controls_freq = pd.read_csv("controls.afreq", sep=r"\s+", index_col=False)
  df_log = pd.read_csv("brat1.logistic", sep=r"\s+", index_col=False)
  df_cases_freq = df_cases_freq[["ID", "ALT1_FREQ"]].rename(columns={"ALT1_FREQ": "CARRIER_ALT1_FREQ"})
  df_controls_freq = df_controls_freq[["ID", "ALT1_FREQ"]].rename(columns={"ALT1_FREQ": "NONCARRIER_ALT1_FREQ"})
  df_merged = df_log.merge(df_cases_freq, on="ID", how="left")
  df_merged = df_merged.merge(df_controls_freq, on="ID", how="left")
  df_merged.to_csv("brat1.logistic", sep="\t", index=False, na_rep="NA")
  CODE
  >>>

  output {
    File out1 = "brat1.logistic"
    File out2 = "cases.afreq"
    File out3 = "controls.afreq"
    File out4 = "brat1_ancestry_adjusted.PHENO.glm.logistic.hybrid"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bioinformatics_tools"
    cpu: 4
    memory: "12 GB"
    preemptible: 3
  }
}

task T5_filter_vcf {
  input {
    File vcf
    File association_results
  }

  command <<<
  set -euxo pipefail

  python3 <<CODE
  import pandas as pd
  df = pd.read_csv("~{association_results}", sep=r"\s+")

  filtered_variants = df[
    (df["P"] < 1e-50) &
    (df["CARRIER_ALT1_FREQ"] >= 0.35) &
    (df["CARRIER_ALT1_FREQ"] <= 0.65)
  ]["ID"].tolist()

  # Always include the BRAT1 variant of interest
  filtered_variants.append("chr7_2555487_C_A")

  # Write output
  with open("filtered_variants.list", "w") as f:
    f.write("\n".join(filtered_variants) + "\n")

  CODE
  bcftools view --threads 4 -i 'ID=@filtered_variants.list' -Oz -o filtered.vcf.bgz ~{vcf}
  >>>

  output {
    File out1 = "filtered.vcf.bgz"
    File out2 = "filtered_variants.list"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bioinformatics_tools"
    cpu: 4
    memory: "12 GB"
    preemptible: 3
  }
}

task T6_beagle_phasing {
  input {
    File vcf
    Int memory_gb = 12
  }

  command <<<
    set -euxo pipefail

    java -Xmx~{memory_gb}g \
      -jar /opt/beagle/beagle.jar \
      gt=~{vcf} \
      out=brat1_phased
  >>>

  output {
    File phased_vcf = "brat1_phased.vcf.gz"
    File log = "brat1_phased.log"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/beagle"
    cpu: 4
    memory: "16 GB"
    preemptible: 3
  }
}

task T7_chunk_genotypes {
  input {
    File vcf
  }

  command <<<
  set -euxo pipefail

  bcftools index -t ~{vcf}
  echo "indexed_file"

  # Prepare File
  bcftools query -f 'ID[\t%SAMPLE]\n' ~{vcf} > samples.tmp
  head -n 1 samples.tmp > tmp.tsv
  rm samples.tmp

  bcftools query -f '%ID[\t%GT]\n' ~{vcf} >> tmp.tsv
  grep -v "chr7_2555487_C_A" tmp.tsv > brat1_genotypes.tsv
  echo "wrote genotypes file"

  n=$(head -n 1 brat1_genotypes.tsv | awk -F'\t' '{print NF-1}')

  for ((start=2; start<=n+1; start+=2000)); do
    end=$((start+1999))
    if [ "$end" -gt "$((n+1))" ]; then
        end=$((n+1))
    fi
    cut -f1,"${start}-${end}" brat1_genotypes.tsv > "chunk_$(( (start-2)/2000 + 1 )).tsv"
  done

  >>>

  output {
    Array[File] out1 = glob("chunk*.tsv")
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bioinformatics_tools"
    memory: "8 GB"
    preemptible: 3
  }
}

task T8_count_snps {
  input {
    File tsv
    File association_results
  }

  command <<<
  python3 <<CODE
  import pandas as pd
  import numpy as np

  df = pd.read_csv("~{tsv}", sep=r"\s+",index_col=False)
  df_assoc = pd.read_csv("~{association_results}", sep=r"\s+",index_col=False)
  print("read file into pandas")

  df = df.merge(df_assoc[["ID", "OR"]], how="left")
  df["OR"] = pd.to_numeric(df["OR"], errors="coerce")

  # Determine which allele to count
  target = np.where(df["OR"] > 1, 1, 0)
  
  # Calculate variant weights
  #weights = np.abs(np.log(df["OR"].values))

  # Extract the two phased alleles
  hap1 = df.iloc[:, 1:-1].map(lambda x: x.split("|")[0])
  hap2 = df.iloc[:, 1:-1].map(lambda x: x.split("|")[1])
  print("step 1")

  hap1 = hap1.apply(pd.to_numeric, errors="coerce")
  hap2 = hap2.apply(pd.to_numeric, errors="coerce")
  print("step 2")

  # Count enriched alleles
  hap1_sum = (hap1.values == target[:, None]).sum(axis=0)
  hap2_sum = (hap2.values == target[:, None]).sum(axis=0)

  #hap1_sum = ((hap1.values == target[:, None]) * weights[:, None]).sum(axis=0)
  #hap2_sum = ((hap2.values == target[:, None]) * weights[:, None]).sum(axis=0)

  #scores = pd.DataFrame([hap1_sum, hap2_sum],columns=df.columns[1:-2])
  scores = pd.DataFrame([hap1_sum, hap2_sum], columns=df.columns[1:-1])
  print ("step 4")
  scores.to_csv("brat1_haplotype_scores.tsv", sep="\t", index=False, header=True)
  CODE
  >>>
  output {
    File out1 = "brat1_haplotype_scores.tsv"
  }
  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/pydata_stack"
    memory: "4 GB"
    preemptible: 3
  }
}

task T10_paste_output {
  input {
    Array[File] haplotype_scores
    File phased_samples
  }
  command <<<
  paste -d '\t' ~{sep=' ' haplotype_scores} > tmp.tsv

  awk '
  {
      for (i = 1; i <= NF; i++) {
          a[i] = (NR == 1 ? $i : a[i] "\t" $i)
      }
  }
  END {
      for (i = 1; i <= NF; i++)
          print a[i]
  }
  ' tmp.tsv > haplotype_scores.tsv

  python3 <<CODE
  import pandas as pd
  
  df = pd.read_csv("haplotype_scores.tsv",sep='\t',header=None,index_col=False,names=["person_id","left_hap","right_hap"],dtype={"person_id": str})
  df_phased = pd.read_csv("~{phased_samples}",sep='\t',index_col=False, dtype={"person_id": str})

  df_merged = pd.merge(df, df_phased, on="person_id",how="left")
  df_merged.to_csv("AoU_BRAT1_haplotype_phased_scores.tsv",sep='\t',index=False) 
  CODE
  >>>
  output {
    File out1 = "AoU_BRAT1_haplotype_phased_scores.tsv"
  }
  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bioinformatics_tools"
    preemptible: 3
  }
}

task T9_find_haplotype{
  input {
    File vcf
    File brat1_het
  }

  command <<<
  set -euxo pipefail

  bcftools index -t ~{vcf}
  bcftools view -S ~{brat1_het} ~{vcf} -Oz -o brat1_het.vcf.gz

  echo -e "person_id\tphase" > brat1_carrier_phase.tsv
  bcftools query \
      -i 'ID="chr7_2555487_C_A"' \
      -f '[%SAMPLE\t%GT\n]' \
      brat1_het.vcf.gz |
  awk '
      $2 == "0|1" {print $1 "\tright"}
      $2 == "1|0" {print $1 "\tleft"}
  ' >> brat1_carrier_phase.tsv

  >>>

  output {
    File out1 = "brat1_carrier_phase.tsv"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bioinformatics_tools"
    memory: "8 GB"
    preemptible: 3
  }
}
