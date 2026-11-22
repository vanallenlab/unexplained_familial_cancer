# The Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow GENE_REGION_FORMAL_ANALYSIS {
  input {
    File vcf
    File vcf_idx
    String gene_name = "BRAT1"
    String chromosome = "chr7"
    Int region_start = 1555487
    Int region_end = 3555487
    Int variant_position = 2555487
    String variant_id = "chr7_2555487_C_A"
  }

  call T0_extract_variant_gt {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      chromosome = chromosome,
      variant_position = variant_position,
      variant_id = variant_id,
      gene_name = gene_name
  }
  
  call T1_create_bed_files {
    input:
      chromosome = chromosome,
      region_start = region_start,
      region_end = region_end,
      variant_position = variant_position
  }

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
        carrier_alt = T0_extract_variant_gt.alt,
        hom_ref = T0_extract_variant_gt.hom_ref,
        gene_name = gene_name
    }
  }

  call Tasks.concatenateFiles_noheader as concat_logreg{
    input:
      files = T3_plink_log_reg.out1,
      callset_name = gene_name
  }
  call T5_filter_vcf {
    input:
      vcf = vcf,
      association_results = concat_logreg.out2,
      variant_id = variant_id,
      gene_name = gene_name
  }
  call T6_beagle_phasing {
    input:
      vcf = T5_filter_vcf.out1,
      gene_name = gene_name
  }
  #call T7_find_haplotype{
  #  input:
  #    vcf = T6_beagle_phasing.phased_vcf,
  #    brat1_het = T0_extract_variant_gt.het,
  #}
  call T7_chunk_genotypes {
    input:
      vcf = T6_beagle_phasing.phased_vcf,
      variant_id = variant_id,
      gene_name = gene_name
  }
  scatter (tsv in T7_chunk_genotypes.out1){
    call T8_count_snps {
      input:
        tsv = tsv,
        association_results = concat_logreg.out2,
        gene_name = gene_name
    }
  }
  call T9_find_haplotype{
    input:
      vcf = T6_beagle_phasing.phased_vcf,
      carrier_het = T0_extract_variant_gt.het,
      variant_id = variant_id,
      gene_name = gene_name
  }
  call T10_paste_output {
    input:
      haplotype_scores = T8_count_snps.out1,
      phased_samples = T9_find_haplotype.out1,
      gene_name = gene_name
  }
}

task T0_extract_variant_gt {
  input {
    File vcf
    File vcf_idx
    String chromosome
    Int variant_position
    String variant_id
    String gene_name
  }

  command <<<
    set -euo pipefail

    # Extract the target variant from the supplied chromosome and position
    bcftools view \
      --threads 4 \
      -r ~{chromosome}:$((~{variant_position} - 1))-$((~{variant_position} + 1)) \
      -i 'ID="~{variant_id}"' \
      -Oz \
      -o ~{gene_name}.vcf.gz \
      ~{vcf}

    # Index the resulting single-variant VCF
    bcftools index --threads 4 -t ~{gene_name}.vcf.gz

    # Extract sample IDs and genotypes
    bcftools query \
      -f '[%SAMPLE\t%GT\n]' \
      ~{gene_name}.vcf.gz > ~{gene_name}_genotypes.tsv

    # Classify genotypes
    awk -F'\t' '
      BEGIN {
        OFS = "\t"
      }
      {
        sample = $1
        gt = $2

        if (gt == "0/0" || gt == "0|0") {
          print sample > "~{gene_name}_hom_ref.list"
        }
        else if (gt == "0/1" || gt == "1/0" ||
                 gt == "0|1" || gt == "1|0") {
          print sample > "~{gene_name}_het.list"
          print sample > "~{gene_name}_alt.list"
        }
        else if (gt == "1/1" || gt == "1|1") {
          print sample > "~{gene_name}_hom_alt.list"
          print sample > "~{gene_name}_alt.list"
        }
        else {
          print sample > "~{gene_name}_other_gt.list"
        }
      }
    ' ~{gene_name}_genotypes.tsv

  cat ~{gene_name}_hom_ref.list ~{gene_name}_alt.list > ~{gene_name}_cohort.list
  >>>

  output {

    File cohort = "~{gene_name}_cohort.list"
    File out_vcf = "~{gene_name}.vcf.gz"
    File out_index = "~{gene_name}.vcf.gz.tbi"

    File hom_ref = "~{gene_name}_hom_ref.list"
    File het = "~{gene_name}_het.list"
    File hom_alt = "~{gene_name}_hom_alt.list"
    File alt = "~{gene_name}_alt.list"
    File other_gt = "~{gene_name}_other_gt.list"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bcftools"
    cpu: 4
    memory: "12 GB"
  }
}

task T1_create_bed_files {
  input {
    String chromosome
    Int region_start
    Int region_end
    Int variant_position
    Int shards = 100
  }

  command <<<
  set -euxo pipefail

    if (( ~{region_start} >= ~{region_end} )); then
      echo "ERROR: region_start must be less than region_end" >&2
      exit 1
    fi
    if (( ~{variant_position} < ~{region_start} || ~{variant_position} > ~{region_end} )); then
      echo "ERROR: variant_position must fall within the supplied region" >&2
      exit 1
    fi

    fixed_start=$((~{variant_position} - 1))
    fixed_end=$((~{variant_position} + 1))

    for ((i=0; i<~{shards}; i++)); do
      shard_start=$((~{region_start} + i * (~{region_end} - ~{region_start}) / ~{shards}))
      shard_end=$((~{region_start} + (i + 1) * (~{region_end} - ~{region_start}) / ~{shards}))

      printf -v shard_num "%03d" "$i"

      printf "~{chromosome}\t%s\t%s\n~{chromosome}\t%s\t%s\n" \
        "$shard_start" "$shard_end" \
        "$fixed_start" "$fixed_end" \
        > "~{chromosome}_shard_${shard_num}.bed"
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

task T3_plink_log_reg {
  input {
    File vcf
    File pcs = "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/AoU_CDRv9_pcs.tsv.gz"
    File carrier_alt
    File hom_ref
    String gene_name
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
      --out ~{gene_name}

    # Create phenotype file (header added for PLINK 2 compatibility)
    echo -e "#IID\tPHENO" > ~{gene_name}.pheno
    awk '{print $1 "\t2"}' ~{carrier_alt} >> ~{gene_name}.pheno
    awk '{print $1 "\t1"}' ~{hom_ref} >> ~{gene_name}.pheno

    # Run ancestry-adjusted log reg association test
    plink2 \
      --pfile ~{gene_name} \
      --pheno ~{gene_name}.pheno \
      --pheno-name PHENO \
      --allow-no-sex \
      --covar pcs.covar \
      --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
      --glm dominant firth-fallback hide-covar omit-ref \
      --out ~{gene_name}_ancestry_adjusted

    plink2 \
      --pfile ~{gene_name} \
      --pheno ~{gene_name}.pheno \
      --pheno-name PHENO \
      --keep-if PHENO == 2 \
      --freq cols=chrom,pos,ref,alt1,alt1freq,nobs \
      --out ~{gene_name}_cases

    plink2 \
      --pfile ~{gene_name} \
      --pheno ~{gene_name}.pheno \
      --pheno-name PHENO \
      --keep-if PHENO == 1 \
      --freq cols=chrom,pos,ref,alt1,alt1freq,nobs \
      --out ~{gene_name}_controls

  awk 'BEGIN {OFS="\t"} {print $3,$6,$7}' ~{gene_name}_cases.afreq > ~{gene_name}_cases_freq.tsv
  awk 'BEGIN {OFS="\t"} {print $3,$6,$7}' ~{gene_name}_controls.afreq > ~{gene_name}_controls_freq.tsv

  cut -f2,3,9,10,12,13,14,15,16,17 ~{gene_name}_ancestry_adjusted.PHENO.glm.logistic.hybrid > ~{gene_name}.logistic

  python3 <<CODE
  import pandas as pd
  df_cases_freq = pd.read_csv("~{gene_name}_cases_freq.tsv", sep=r"\s+", index_col=False)
  df_controls_freq = pd.read_csv("~{gene_name}_controls_freq.tsv", sep=r"\s+", index_col=False)
  df_log = pd.read_csv("~{gene_name}.logistic", sep=r"\s+", index_col=False)
  df_cases_freq = df_cases_freq[["ID", "ALT1_FREQ"]].rename(columns={"ALT1_FREQ": "CARRIER_ALT1_FREQ"})
  df_controls_freq = df_controls_freq[["ID", "ALT1_FREQ"]].rename(columns={"ALT1_FREQ": "NONCARRIER_ALT1_FREQ"})
  df_merged = df_log.merge(df_cases_freq, on="ID", how="left")
  df_merged = df_merged.merge(df_controls_freq, on="ID", how="left")
  df_merged.to_csv("~{gene_name}.logistic", sep="\t", index=False, na_rep="NA")
  CODE
  >>>

  output {
    File out1 = "~{gene_name}.logistic"
    File out2 = "~{gene_name}_cases_freq.tsv"
    File out3 = "~{gene_name}_controls_freq.tsv"
    File out4 = "~{gene_name}_ancestry_adjusted.PHENO.glm.logistic.hybrid"
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
    String variant_id
    String gene_name
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

  # Always include the specified variant of interest
  filtered_variants.append("~{variant_id}")

  # Write output
  with open("~{gene_name}_filtered_variants.list", "w") as f:
    f.write("\n".join(filtered_variants) + "\n")

  CODE
  bcftools view --threads 4 -i 'ID=@~{gene_name}_filtered_variants.list' -Oz -o ~{gene_name}_filtered.vcf.bgz ~{vcf}
  >>>

  output {
    File out1 = "~{gene_name}_filtered.vcf.bgz"
    File out2 = "~{gene_name}_filtered_variants.list"
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
    String gene_name
    Int memory_gb = 12
  }

  command <<<
    set -euxo pipefail

    java -Xmx~{memory_gb}g \
      -jar /opt/beagle/beagle.jar \
      gt=~{vcf} \
      out=~{gene_name}_phased
  >>>

  output {
    File phased_vcf = "~{gene_name}_phased.vcf.gz"
    File log = "~{gene_name}_phased.log"
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
    String variant_id
    String gene_name
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
  grep -v "~{variant_id}" tmp.tsv > ~{gene_name}_genotypes.tsv
  echo "wrote genotypes file"

  n=$(head -n 1 ~{gene_name}_genotypes.tsv | awk -F'\t' '{print NF-1}')

  for ((start=2; start<=n+1; start+=2000)); do
    end=$((start+1999))
    if [ "$end" -gt "$((n+1))" ]; then
        end=$((n+1))
    fi
    cut -f1,"${start}-${end}" ~{gene_name}_genotypes.tsv > "~{gene_name}_chunk_$(( (start-2)/2000 + 1 )).tsv"
  done

  >>>

  output {
    Array[File] out1 = glob("~{gene_name}_chunk*.tsv")
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
    String gene_name
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
  scores.to_csv("~{gene_name}_haplotype_scores.tsv", sep="\t", index=False, header=True)
  CODE
  >>>
  output {
    File out1 = "~{gene_name}_haplotype_scores.tsv"
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
    String gene_name
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
  df_merged.to_csv("AoU_~{gene_name}_haplotype_phased_scores.tsv",sep='\t',index=False) 
  CODE
  >>>
  output {
    File out1 = "AoU_~{gene_name}_haplotype_phased_scores.tsv"
  }
  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bioinformatics_tools"
    preemptible: 3
  }
}

task T9_find_haplotype{
  input {
    File vcf
    File carrier_het
    String variant_id
    String gene_name
  }

  command <<<
  set -euxo pipefail

  bcftools index -t ~{vcf}
  bcftools view -S ~{carrier_het} ~{vcf} -Oz -o ~{carrier_het}.vcf.gz

  echo -e "person_id\tphase" > ~{gene_name}_carrier_phase.tsv
  bcftools query \
      -i 'ID="~{variant_id}"' \
      -f '[%SAMPLE\t%GT\n]' \
      ~{carrier_het}.vcf.gz |
  awk '
      $2 == "0|1" {print $1 "\tright"}
      $2 == "1|0" {print $1 "\tleft"}
  ' >> ~{gene_name}_carrier_phase.tsv

  >>>

  output {
    File out1 = "~{gene_name}_carrier_phase.tsv"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/bioinformatics_tools"
    memory: "8 GB"
    preemptible: 3
  }
}
