# The Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow STEP_R4_CDRv9_PRS {
  input {
    String vcf_dir = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/vcf"
    String GOOGLE_PROJECT = "wb-cordial-diamond-9893"
    Int lines_per_file
    Int n_batches_to_test
    String output_dir
    #Array[File] test_vcf_list = ["gs://dataproc-temp-wb-cordial-diamond-9893/prs_test_vcfs.list"]
    Array[String] PGS_IDS
  }

  call T1_gather_cdrv9_vcfs {
    input:
      vcf_dir = vcf_dir,
      lines_per_file = lines_per_file,
      n_batches_to_test = n_batches_to_test,
      GOOGLE_PROJECT = GOOGLE_PROJECT
  }
      

  # Process PGS file
  scatter (PGS_ID in PGS_IDS){  
    call T1_prepare_scores_file {
      input:
        prs_file = "gs://rw-migration-aou-rw-cf3e7774/UFC_REFERENCE_FILES/pgs_files/~{PGS_ID}_hmPOS_GRCh38.txt.gz",
        PGS_ID = PGS_ID
    }
  }

  call T1_merge_bed_files {
    input:
      mini_beds = T1_prepare_scores_file.out2
  }
  #scatter (vcf_list in test_vcf_list){
  scatter(vcf_list in T1_gather_cdrv9_vcfs.out1){
    call T2_calculate_prs {
      input:
        vcf_list = vcf_list,
        prs_files = T1_prepare_scores_file.out1,
        bed_file = T1_merge_bed_files.out1,
        GOOGLE_PROJECT = GOOGLE_PROJECT
    }
  }

  call T4_combine_prs{
    input:
      prs_files = select_all(T2_calculate_prs.out1)
  }

  #call Tasks.concatenateFiles {
  #  input:
  #    files = T2_calculate_prs.out2,
  #    output_name = "prs_failure_log" 
  #}

  call Tasks.copy_file_to_storage as copy1{
    input:
      text_file = T4_combine_prs.out1,
      output_dir = output_dir 
  }

  #call Tasks.copy_file_to_storage as copy2{
  #  input:
  #    text_file = concatenateFiles.out1,
  #    output_dir = output_dir
  #}

}

task T1_gather_cdrv9_vcfs {
  input {
    String vcf_dir
    String GOOGLE_PROJECT
    Int lines_per_file
    Int n_batches_to_test = 0
  }

  command <<<
    set -euxo pipefail

    gsutil -u ~{GOOGLE_PROJECT} ls "~{vcf_dir}/*.bgz" > all_vcfs.list

    split -l ~{lines_per_file} -d -a 8 all_vcfs.list vcf_batch_

    for f in vcf_batch_*; do
      n=${f#vcf_batch_}
      mv "$f" "vcf_batch_${n}.list"
    done

    if [[ ~{n_batches_to_test} -gt 0 ]]; then
      ls vcf_batch*.list | tail -n +$((~{n_batches_to_test} + 1)) | xargs -r rm
    fi

    #for f in vcf_batch_*.list; do
    #    [[ "$f" == "vcf_batch_00000105.list" ||
    #       "$f" == "vcf_batch_00000106.list" ||
    #       "$f" == "vcf_batch_00000455.list" ]] || rm "$f"
    #done

    echo "Total VCFs: $(wc -l < all_vcfs.list)"
    echo "Total batches: $(ls vcf_batch*.list | wc -l)"
  >>>

  output {
    Array[File] out1 = glob("vcf_batch*.list")
  }

  runtime {
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
  }
}

task T1_prepare_scores_file {
  input{
    File prs_file
    String PGS_ID
  }
  command <<<
  python3 <<CODE
  import pandas as pd

  usecols = ['hm_source','hm_chr', 'hm_pos', 'effect_allele', 'other_allele', 'effect_weight']

  df = pd.read_csv("~{prs_file}", comment = "#", index_col=False, sep='\t',dtype={'hm_chr': str, 'hm_pos': str},usecols=usecols)
  df = df.dropna(subset=['hm_chr', 'hm_pos','effect_allele', 'other_allele', 'effect_weight'])
  #df = df[df['hm_source'] == "liftover"]

  # Create BED file for GRCh38 positions
  bed_df = df[['hm_chr', 'hm_pos']].drop_duplicates().copy()
  bed_df['hm_chr'] = 'chr' + bed_df['hm_chr'].astype(str)
  bed_df['start'] = bed_df['hm_pos'].astype(int) - 1
  bed_df['end'] = bed_df['hm_pos'].astype(int) + 1
  bed_df[['hm_chr', 'start', 'end']].to_csv("~{PGS_ID}.bed", sep="\t", index=False, header=False)

  # Create primary and flipped IDs
  df1 = df[['hm_chr', 'hm_pos', 'effect_allele', 'other_allele', 'effect_weight']].copy()
  df1['ID'] = "chr" + df1['hm_chr'] + "_" + df1['hm_pos'] + "_" + df1['effect_allele'] + "_" + df1['other_allele']

  df2 = df[['hm_chr', 'hm_pos', 'effect_allele', 'other_allele', 'effect_weight']].copy()
  df2['ID'] = "chr" + df2['hm_chr'] + "_" + df2['hm_pos'] + "_" + df2['other_allele'] + "_" + df2['effect_allele']

  # Combine both versions
  score_df = pd.concat([
      df1[['ID', 'effect_allele', 'effect_weight']],
      df2[['ID', 'effect_allele', 'effect_weight']]
  ])

  # Drop duplicates just in case
  score_df = score_df.drop_duplicates()

  # Save to file
  score_df.to_csv("~{PGS_ID}.score", sep="\t", index=False, header=False)
  CODE
  >>>
  output {
    File out1 = "~{PGS_ID}.score"
    File out2 = "~{PGS_ID}.bed"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
    memory: "8GB"
  }
}


task T1_merge_bed_files {
  input {
    Array[File] mini_beds
  }

  command <<<
  set -euxo pipefail

  cat ~{sep=' ' mini_beds} | awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print}' > combined.unsorted.bed

  printf 'chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchr20\nchr21\nchr22\nchrX\nchrY\n' > genome.txt

  bedtools sort -i combined.unsorted.bed -g genome.txt > combined.sorted.bed

  bedtools merge -i combined.sorted.bed > combined.bed
  >>>

  output {
    File out1 = "combined.bed"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    memory: "16 GB"
  }
}

task T2_calculate_prs {
  input {
    File vcf_list
    Array[File] prs_files
    File bed_file
    String GOOGLE_PROJECT
  }

  command <<<
    set -euo pipefail

    while read -r vcf_file; do
      vcf_basename=$(basename "$vcf_file" .vcf.bgz)


      # 1. Copy VCF Locally
      gsutil -u ~{GOOGLE_PROJECT} cp "$vcf_file" input.vcf.gz

      bcftools view \
          -T ~{bed_file} \
          -Ou input.vcf.gz |
      bcftools annotate \
          -x FORMAT,^FORMAT/GT \
          -Ou |
      bcftools norm \
          -m -any \
          -Ou |
      bcftools annotate \
          --set-id '%CHROM\_%POS\_%REF\_%ALT' \
          -Oz -o tmp.vcf.gz

      # Check whether the VCF contains any variant records.
      n_variants=$(bcftools view -H tmp.vcf.gz | wc -l)

      if [[ "$n_variants" -eq 0 ]]; then
          echo "No variants in $vcf_file; skipping scoring."
          continue
      fi

      mkdir -p SCORES

      for prs_file in ~{sep=' ' prs_files}; do

        prs_basename=$(basename "$prs_file" | sed 's/\.[^.]*$//')

        if plink --vcf tmp.vcf.gz \
          --score "$prs_file" sum \
          --double-id \
          --out tmp \
          > "plink.log" 2>&1; then
          :
        else
          exit_code=$?
          {
            echo "FAILED (exit code ${exit_code}): VCF=${vcf_basename}, PGS=${prs_basename}"
            cat plink.log
            echo "----------------------------------------"
           } >> plink_failures.log
        fi

        cat plink.log
        # If plink returned anything (were there overlapping snps between vcf and score file)
        if [ -f tmp.profile ]; then

          output_file="SCORES/${prs_basename}.intermediate.pgs"

          if [ ! -f "$output_file" ]; then

              # Initialize file with person IDs and the first score column
              awk '
                NR == 1 { printf "person_id\tPGS_col\n"; next }
                { print $2 "\t" $6 }
              ' tmp.profile > "$output_file"

          else

              # Append the score column using the same header
              awk '
                NR == 1 { print "PGS_col"; next }
                { print $6 }
              ' tmp.profile > tmp.score

              paste "$output_file" tmp.score > tmp.combined
              mv tmp.combined "$output_file"

              rm -f tmp.score

          fi

        fi

        rm -f \
          tmp.profile \
          tmp.nosex \
          tmp.nopred \
          tmp.log

      done

      rm -f  tmp.vcf.gz

    done < ~{vcf_list}

    # Put all Sharded scores into one column (per PRS)
    #mkdir -p SCORES_SUM
    #for input_file in SCORES/*.intermediate.pgs; do

        # Handle the case where no matching files exist
    #    [ -f "$input_file" ] || continue
         
    #    filename=$(basename "$input_file")
    #    pgs_id="${filename%%.intermediate.pgs}"
    #    output_file="SCORES_SUM/$filename"

    #    awk -v pgs_id="$pgs_id" '
    #        BEGIN { OFS = "\t" }

    #        NR == 1 {
    #            print "person_id", pgs_id
    #            next
    #        }

    #        {
    #            sum = 0
    #            for (i = 2; i <= NF; i++) {
    #                sum += $i
    #            }
    #            print $1, sum
    #        }
    #    ' "$input_file" > "$output_file"

    #   rm "$input_file"
    #done

    python3 <<CODE
    import glob
    import os

    os.makedirs("SCORES_SUM", exist_ok=True)

    input_files = sorted(glob.glob("SCORES/*.intermediate.pgs"))

    for filename in input_files:
        pgs_id = os.path.basename(filename).removesuffix(".intermediate.pgs")
        output_file = os.path.join("SCORES_SUM", os.path.basename(filename))

        total_scores = {}

        with open(filename) as f:
            next(f)  # Skip header

            for line in f:
                fields = line.rstrip("\n").split("\t")
                sample = fields[0]

                # Sum all shard-level score columns for this person
                for score_str in fields[1:]:
                    score = float(score_str)
                    total_scores[sample] = total_scores.get(sample, 0.0) + score

        # Write one summed score per person
        with open(output_file, "w") as out:
            out.write(f"person_id\t{pgs_id}\n")

            for sample, score in sorted(total_scores.items()):
                out.write(f"{sample}\t{round(score, 16)}\n")

        os.remove(filename)

    CODE

    # Put all PGS scores into one file, if any exist
    shopt -s nullglob
    score_files=(SCORES_SUM/*.intermediate.pgs)

    if (( ${#score_files[@]} > 0 )); then
        output_file="SCORES_ALL.tsv"
        first=true

        for file in "${score_files[@]}"; do
            if [ "$first" = true ]; then
                cp "$file" "$output_file"
                first=false
            else
                cut -f2 "$file" > tmp.score
                paste "$output_file" tmp.score > tmp.combined
                mv tmp.combined "$output_file"
                rm tmp.score
            fi
        done

        gzip "$output_file"
        echo "Combined and compressed PGS scores."
    else
        echo "WARNING: No summed PGS files found; skipping combination and gzip."
    fi
    # Put all PGS scores into one file
    #output_file="SCORES_ALL.tsv"
    #first=true

    #for file in $(ls SCORES_SUM/*.intermediate.pgs | sort); do
    #    if [ "$first" = true ]; then
    #        cp "$file" "$output_file"
    #        first=false
    #    else
    #        cut -f2 "$file" > tmp.score
    #        paste "$output_file" tmp.score > tmp.combined
    #        mv tmp.combined "$output_file"
    #        rm tmp.score
    #    fi
    #done

    #gzip SCORES_ALL.tsv
  >>>

  output {
    File? out1 = "SCORES_ALL.tsv.gz"
    #File? out2 = "plink_failures.log"
  }

  runtime {
    docker: "vanallenlab/bioinformatics_tools"
    preemptible: 3
    memory: "8 GB"
  }
}

task T2_calculate_prs_og {
  input {
    File vcf_list
    Array[File] prs_files
    File sample_list
    String GOOGLE_PROJECT
  }

  command <<<
    set -euo pipefail

    while read -r vcf_file; do
      vcf_basename=$(basename "$vcf_file" .vcf.bgz)


      # 1. Copy VCF locally
      gsutil -u ~{GOOGLE_PROJECT} cp "$vcf_file" input.vcf.gz

      bcftools view \
          -S ~{sample_list} \
          -Ou input.vcf.gz |
      bcftools annotate \
          -x FORMAT,^FORMAT/GT \
          -Ou |
      bcftools norm \
          -m -any \
          -Ou |
      bcftools annotate \
          --set-id '%CHROM\_%POS\_%REF\_%ALT' \
          -Oz -o tmp.vcf.gz

      for prs_file in ~{sep=' ' prs_files}; do

        prs_basename=$(basename "$prs_file" | sed 's/\.[^.]*$//')

        # Make directory to organize score files
        mkdir -p "${prs_basename}"

        plink --vcf tmp.vcf.gz \
          --score "$prs_file" sum \
          --double-id \
          --out tmp \
          > /dev/null 2>&1 || true

        if [ -f tmp.profile ]; then
          awk '{print $2 "\t" $6}' tmp.profile \
            > "${prs_basename}/${vcf_basename}.${prs_basename}.profile"
        else
          echo -e "IID\tSCORE" \
            > "${prs_basename}/${vcf_basename}.${prs_basename}.profile"
        fi

        rm -f \
          tmp.profile \
          tmp.nosex \
          tmp.nopred \
          tmp.log

      done

      rm -f  tmp.vcf.gz

    done < ~{vcf_list}

    # Archive results for each PRS
    for prs_file in ~{sep=' ' prs_files}; do

      prs_basename=$(basename "$prs_file" | sed 's/\.[^.]*$//')

      tar -czf "${prs_basename}.tar.gz" "${prs_basename}/"

    done
  >>>

  output {
    Array[File] out1 = glob("*.tar.gz")
  }

  runtime {
    docker: "vanallenlab/bioinformatics_tools"
    preemptible: 3
    memory: "8 GB"
  }
}

task T3_sum_prs_batch_scores {
  input {
    Array[File] zipped_prs_files
  }
  command <<<
    
  set -euo pipefail

  # Extract all PGS archives
  for archive in ~{sep=' ' zipped_prs_files}; do
    tar -xzf "$archive"
  done

  python3 <<CODE
  import glob
  import os

  # Loop through each PGS directory
  pgs_dirs = sorted(glob.glob("PGS*/"))

  for pgs_dir in pgs_dirs:
    pgs_id = os.path.basename(os.path.normpath(pgs_dir))

    # Find all shard-level profile files for this PGS
    pgs_files = sorted(glob.glob(os.path.join(pgs_dir, "*.profile")))

    # This will hold cumulative scores
    total_scores = {}

    for filename in pgs_files:
        with open(filename) as f:
            next(f)  # skip header

            for line in f:
                sample, score_str = line.strip().split('\t')
                score = float(score_str)

                total_scores[sample] = (total_scores.get(sample, 0.0) + score)

    # Write final output
    with open(f"{pgs_id}.intermediate.pgs", "w") as out:
        out.write(f"person_id\t{pgs_id}\n")

        for sample, score in sorted(total_scores.items()):
            out.write(f"{sample}\t{round(score, 16)}\n")
  CODE

  # Compress final PGS files
  for pgs_file in *.intermediate.pgs; do
      tar -czf "${pgs_file}.tar.gz" "$pgs_file"
      rm -f "$pgs_file"
  done
  >>>
  output {
    Array[File] out1 = glob("*.intermediate.pgs.tar.gz")
  }
  
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task T3_sum_intermediate_scores {
  input {
    File scores_archive
  }

  command <<<
    set -euo pipefail

    # Extract the archive
    tar -xzf ~{scores_archive}

    python3 <<CODE
    import glob
    import os
    import math

    pgs_files = sorted(glob.glob("SCORES/*.intermediate.pgs"))

    all_scores = {}

    for filename in pgs_files:
        pgs_id = os.path.basename(filename).replace(".intermediate.pgs", "")
        print(f"Processing {pgs_id}")

        total_scores = {}

        with open(filename) as f:
            header = next(f).strip().split("\t")
            score_columns = len(header) - 1

            for line in f:
                if not line.strip():
                    continue

                fields = line.strip().split("\t")
                person_id = fields[0]
                scores = [float(x) for x in fields[1:] if x != "NA"]

                total_scores[person_id] = math.fsum(scores)

        for person_id, score in total_scores.items():
            if person_id not in all_scores:
                all_scores[person_id] = {}
            all_scores[person_id][pgs_id] = score

        print(f"  Processed {len(total_scores)} people")
        print(f"  Summed {score_columns} score columns")

    # Write one combined file
    pgs_ids = [
        os.path.basename(f).replace(".intermediate.pgs", "")
        for f in pgs_files
    ]

    with open("batch_prs.tsv", "w") as out:
        out.write("person_id\t" + "\t".join(pgs_ids) + "\n")

        for person_id in sorted(all_scores):
            scores = [
                f"{all_scores[person_id][pgs_id]:.16f}"
                for pgs_id in pgs_ids
            ]
            out.write(person_id + "\t" + "\t".join(scores) + "\n")

    CODE

    gzip batch_prs.tsv
  >>>

  output {
    File out1 = "batch_prs.tsv.gz"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}


task T4_combine_prs {
  input {
    Array[File] prs_files
  }

  command <<<
  set -euo pipefail

  python3 <<CODE
  import gzip

  prs_files = "~{sep=' ' prs_files}".split()

  all_scores = {}
  pgs_ids = []

  for filename in prs_files:
      print(f"Processing {filename}", flush=True)

      with gzip.open(filename, "rt") as f:
          header = f.readline().strip().split("\t")
          file_pgs_ids = header[1:]

          # Track every unique PGS column across files
          for pgs_id in file_pgs_ids:
              if pgs_id not in pgs_ids:
                  pgs_ids.append(pgs_id)

          # Read and accumulate scores
          for line in f:
              if not line.strip():
                  continue

              fields = line.rstrip("\n").split("\t")
              person_id = fields[0]
              scores = fields[1:]

              if person_id not in all_scores:
                  all_scores[person_id] = {}

              for i in range(len(file_pgs_ids)):
                  pgs_id = file_pgs_ids[i]
                  score = float(scores[i])

                  if pgs_id not in all_scores[person_id]:
                      all_scores[person_id][pgs_id] = 0.0

                  all_scores[person_id][pgs_id] += score


  # Write the combined output
  with gzip.open("CDRv9_prs.tsv.gz", "wt") as out:
      out.write("person_id\t" + "\t".join(pgs_ids) + "\n")

      for person_id in sorted(all_scores):
          scores = [ f"{all_scores[person_id].get(pgs_id, 0.0):.16f}" for pgs_id in pgs_ids]

          out.write(person_id + "\t" + "\t".join(scores) + "\n")

  print(f"Processed {len(prs_files)} files")
  print(f"Found {len(pgs_ids)} unique PGSs")
  print(f"Processed {len(all_scores)} people")
  CODE
  >>>

  output {
    File out1 = "CDRv9_prs.tsv.gz"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
    disks: "local-disk 50 HDD" 
  }
}
