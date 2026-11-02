#!/usr/bin/env bash

# The Van Allen Lab Unexplained Familial Cancers (UFC) Study
# Copyright (c) 2025-Present, Ryan L. Collins, Noah Fields, and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to refine GATK-SV outputs for all UFC samples

# Note that the below code makes *heavy* use of code from the DFCI G2C Database
# This code is available/documented here:
# https://github.com/collins-genomics/pancan_germline_wgs

# Note also that this code was ported from the All of Us Workbench v1.0 to v2.0 (Verily Pre)


##############################
# SETUP INTERACTIVE TERMINAL #
##############################

# This setup makes heavy assumptions about standards and helper functions defined 
# for the G2C project, which we repurpose here for convenience

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://rw-migration-aou-rw-84a0039b
export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-db439520"
gcloud storage cp $MAIN_WORKSPACE_BUCKET/code/scripts/configure_verily_vm.sh ./ && \
. configure_verily_vm.sh && \
rm configure_verily_vm.sh
unset WN

# Set up local directory structure
for dir in scratch staging cromshell cromshell/inputs cromshell/job_ids cromshell/progress; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done

# Copy legacy mingq code
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/misc/legacy_mingq_wdl code/wdl/gatk-sv/

# Create dependencies .zip for GATK-SV module submissions
cd code/wdl/gatk-sv && \
zip -r gatksv.dependencies.zip . && \
mv gatksv.dependencies.zip ~/ && \
cd ~

# Create dependencies .zip for minGQ module submissions
cd code/wdl/gatk-sv/legacy_mingq_wdl && \
zip -r mingq.dependencies.zip . && \
mv mingq.dependencies.zip ~/ && \
cd ~

# Make tarball with SV annotation dependencies
cd code/wdl/gatk-sv/svannotation_v1.1 && \
zip svannotation.dependencies.zip *.wdl && \
mv svannotation.dependencies.zip ~/ && \
cd ~

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json


##########################
# LAUNCH CROMWELL SERVER #
##########################

# Launch Cromwell in server mode on the same Verily VM
# Note that this must be executed from a separate terminal after 
# configuring the environment per the above

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://rw-migration-aou-rw-84a0039b
export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-db439520"

cd ~

# Source .bashrc and bash utility functions
. code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh

# Ensure we are using Cromwell v92 to expose default machine types
if ! [ -e ~/bin/cromwell-92.jar ]; then
  mkdir ~/bin
  wget -P ~/bin/ https://github.com/broadinstitute/cromwell/releases/download/92/cromwell-92.jar
fi
export CROMWELL_JAR=~/bin/cromwell-92.jar

# Set terminal timezone to Eastern US time
export TZDIR=$(python -c "import tzdata, os; print(os.path.join(os.path.dirname(tzdata.__file__), 'zoneinfo'))")
export TZ=America/New_York
date +"%Y-%m-%d %H:%M:%S %Z %z"

# Reassign $WORKSPACE_BUCKET to match AoU v1-migrated URIs
export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-db439520"

# Create cromwell config file
if ! [ -e /home/jupyter/.cromwell ]; then
  mkdir /home/jupyter/.cromwell
fi
if [ -e /home/jupyter/.cromwell/cromwell.conf ]; then
  rm /home/jupyter/.cromwell/cromwell.conf
  rm /home/jupyter/.cromwell/cromwell.override.conf
fi
wb cromwell generate-config \
  --google-bucket-name=$WORKSPACE_BUCKET \
  --dir=/home/jupyter/.cromwell
cat << EOF > /home/jupyter/.cromwell/cromwell.override.conf
include "cromwell.conf"

google {
  cloud-sdk-image-url = "gcr.io/google.com/cloudsdktool/google-cloud-cli:slim"
}

backend {
  default = "GCPBATCH"
  providers {

    # Disables the Local backend
    Local.config.root = "/dev/null"

    GCPBATCH {
     
        actor-factory = "cromwell.backend.google.batch.GcpBatchBackendLifecycleActorFactory"

      config {
        project = "${GOOGLE_PROJECT}"
        concurrent-job-limit = 3000
        root = "${WORKSPACE_BUCKET}/workflows/cromwell-executions"

                        
          virtual-private-cloud {
          network-name = "projects/${GOOGLE_PROJECT}/global/networks/network"
          subnetwork-name = "projects/${GOOGLE_PROJECT}/regions/*/subnetworks/subnetwork"
                   
        }

        batch {
          auth = "application_default"
          compute-service-account = "${PET_SA_EMAIL}"
          location = "us-central1"
        }

        default-runtime-attributes {
          noAddress: true
        }
        
        filesystems {
          gcs {
              auth = "application_default"
          }
        }
      }
    }
  }
}

call-caching {
  enabled = true
}

system {
  job-rate-control {
    jobs = 10
    per = 10 seconds
  }
  workflow-heartbeats {
    ttl = 20 minutes
    write-failure-shutdown-duration = 15 minutes
    write-batch-size = 250
  }
}

engine {
  filesystems {
    gcs {
      auth = "application_default"
    }
  }
}

database {
  profile = "slick.jdbc.HsqldbProfile$"

  db {
    driver = "org.hsqldb.jdbcDriver"

    url = """
      jdbc:hsqldb:file:/home/jupyter/.cromwell/db/cromwell;
      shutdown=false;
      hsqldb.default_table_type=cached;
      hsqldb.tx=mvcc;
      hsqldb.result_max_memory_rows=10000;
      hsqldb.large_data=true;
      hsqldb.script_format=3
    """

    connectionTimeout = 1200000
    numThreads = 8
  }

  insert-batch-size = 2000
}
EOF

# Ensure cromwell local database exists
DBDIR=/home/jupyter/.cromwell/db
mkdir -p $DBDIR
if [[ -d "$DBDIR" ]]; then
  size_gb=$( du -sBG "$DBDIR" | cut -f1 | tr -d 'G' )
  echo "Cromwell DB size: $size_gb GB"
fi

# Create cromshell config file
if [ ! -e /home/jupyter/.cromshell ]; then
  mkdir /home/jupyter/.cromshell
fi
cat << EOF > /home/jupyter/.cromshell/cromshell_config.json
{
  "cromwell_server": "http://localhost:8000",
  "requests_timeout": 5
}
EOF

# Launch cromwell in server mode
echo -e "\nNow launching cromwell server using $CROMWELL_JAR"
echo -e "Cromwell version detected as: $( java -jar $CROMWELL_JAR --version )"
java \
  -Xms12G \
  -Xmx48G \
  -XX:+UseG1GC \
  -Dconfig.file=/home/jupyter/.cromwell/cromwell.override.conf \
  -Xlog:gc*:file=/home/jupyter/.cromwell/gc.log:time \
  -XX:+HeapDumpOnOutOfMemoryError \
  -XX:HeapDumpPath=/home/jupyter/.cromwell/ \
  -jar $CROMWELL_JAR \
  server


#################
# UFC FUNCTIONS #
#################
# Resets a QC plotting directory
reset_qc_plot_dir() {
  # Check inputs
  if [ $# -lt 3 ]; then
    echo "Must provide staging directory, QC workflow suffix, and QC metric staging directory as positional arguments"
    return 2
  fi
  staging_dir=$1
  wid_suffix=$2
  gs_metric_dir=$3
  unset prev_stats
  if [ $# -eq 4 ]; then
    prev_stats=$4
  fi

  cat << EOF > $staging_dir/main_keys.list
size_distrib
af_distrib
size_vs_af_distrib
all_svs_bed
common_svs_bed
genotype_distrib
ld_stats
EOF

  cat << EOF > $staging_dir/bench_keys.list
site_benchmark_ppv_by_freqs
site_benchmark_sensitivity_by_freqs
site_benchmark_common_sv_ppv_beds
site_benchmark_common_sv_sens_beds
twin_genotype_benchmark_distribs
trio_mendelian_violation_distribs
EOF

  # Clear old input arrays
  while read key; do  
    fname=$staging_dir/$key.uris.list
    if [ -e $fname ]; then rm $fname; fi
  done < $staging_dir/main_keys.list
  while read key; do
    for subset in giab_easy giab_hard; do
      fname=$staging_dir/$key.$subset.uris.list
      if [ -e $fname ]; then rm $fname; fi
    done
  done < $staging_dir/bench_keys.list
  for suffix in af_distribution size_distribution; do
    fname=$staging_dir/gnomAD_$suffix.uris.list
    if [ -e $fname ]; then rm $fname; fi
  done
  for key in sample_benchmark_ppv_distribs sample_benchmark_sensitivity_distribs; do
    for dset in external_srwgs external_lrwgs; do
      for subset in giab_easy giab_hard; do
        fname=$staging_dir/$key.$subset.$dset.uris.list
        if [ -e $fname ]; then rm $fname; fi
      done
    done
  done

  # Build input arrays
  for k in $( seq 1 22 ) X Y; do
    
    # Localize output tracker json and get URIs for QC metrics
    json_fname=Collect$wid_suffix.chr$k.outputs.json
    gsutil cp $gs_metric_dir/chr$k/$json_fname $staging_dir/
    while read key; do
      jq .\"CollectVcfQcMetrics.$key\" $staging_dir/$json_fname \
      | fgrep -xv "null" | tr -d '"' \
      >> $staging_dir/$key.uris.list
    done < $staging_dir/main_keys.list
    while read key; do
      for subset in giab_easy giab_hard; do
        jq .\"CollectVcfQcMetrics.$key\" $staging_dir/$json_fname \
        | fgrep -xv "null" | tr -d '"[]' | sed 's/,$/\n/g' \
        | sed '/^$/d' | awk '{ print $1 }' | fgrep $subset \
        >> $staging_dir/$key.$subset.uris.list
      done
    done < $staging_dir/bench_keys.list

    # Due to delisting behavior of manage_chromshards.py, external sample benchmark
    # results need to be parsed in a custom manner as below
    for key in sample_benchmark_ppv_distribs sample_benchmark_sensitivity_distribs; do
      for dset in external_srwgs external_lrwgs; do
        for subset in giab_easy giab_hard; do
          jq .\"CollectVcfQcMetrics.$key\" $staging_dir/$json_fname \
          | fgrep -xv "null" | tr -d '"[]' | sed 's/,$/\n/g' \
          | sed '/^$/d' | awk '{ print $1 }' | fgrep $dset | fgrep $subset \
          >> $staging_dir/$key.$subset.$dset.uris.list
        done
      done
    done

    # Clear local copy of output tracker json
    rm $staging_dir/$json_fname

    # Add precomputed gnomAD v4.1 reference distributions to file lists
    for suffix in af_distribution size_distribution; do
      echo "gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/chr$k/gnomad.v4.1.gatksv.chr$k.$suffix.merged.tsv.gz" \
      >> $staging_dir/gnomAD_$suffix.uris.list
    done
  done

# Write input .json
cat << EOF | python -m json.tool > $staging_dir/Plot$wid_suffix.inputs.json
{
  "PlotVcfQcMetrics.af_distribution_tsvs": $( collapse_txt $staging_dir/af_distrib.uris.list ),
  "PlotVcfQcMetrics.all_sv_beds": $( collapse_txt $staging_dir/all_svs_bed.uris.list ),
  "PlotVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PlotVcfQcMetrics.benchmark_interval_names": ["Easy", "Hard"],
  "PlotVcfQcMetrics.common_af_cutoff": 0.01,
  "PlotVcfQcMetrics.common_sv_beds": $( collapse_txt $staging_dir/common_svs_bed.uris.list ),
  "PlotVcfQcMetrics.custom_qc_target_metrics": "$WORKSPACE_BUCKET/data/misc/dfci-ufc.sv.qc_targets.tsv",
  "PlotVcfQcMetrics.deduplicate": true,
  "PlotVcfQcMetrics.g2c_analysis_docker": "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/g2c_analysis:ec08fc1",
  "PlotVcfQcMetrics.linux_docker": "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/ubuntu:plucky-20251001",
  "PlotVcfQcMetrics.output_prefix": "dfci-ufc.sv.v2.$wid_suffix",
  "PlotVcfQcMetrics.peak_ld_stat_tsvs": $( collapse_txt $staging_dir/ld_stats.uris.list ),
  "PlotVcfQcMetrics.PlotSiteBenchmarking.gcp_machine_type": "n2d-standard-8",
  "PlotVcfQcMetrics.PlotSiteMetrics.gcp_machine_type": "n2d-standard-8",
  "PlotVcfQcMetrics.ref_af_distribution_tsvs": $( collapse_txt $staging_dir/gnomAD_af_distribution.uris.list ),
  "PlotVcfQcMetrics.ref_size_distribution_tsvs": $( collapse_txt $staging_dir/gnomAD_size_distribution.uris.list ),
  "PlotVcfQcMetrics.ref_cohort_prefix": "gnomAD_v4.1_gatksv",
  "PlotVcfQcMetrics.ref_cohort_plot_title": "gnomAD-SV v4.1",
  "PlotVcfQcMetrics.sample_ancestry_labels": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.qc_ancestry.tsv",
  "PlotVcfQcMetrics.sample_phenotype_labels": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.qc_phenotype.tsv",
  "PlotVcfQcMetrics.sample_benchmark_dataset_prefixes": ["external_srwgs", "external_lrwgs"],
  "PlotVcfQcMetrics.sample_benchmark_dataset_titles": ["external srWGS", "external lrWGS"],
  "PlotVcfQcMetrics.sample_benchmark_ppv_distribs": [[ $( collapse_txt $staging_dir/sample_benchmark_ppv_distribs.giab_easy.external_srwgs.uris.list ),
                                                       $( collapse_txt $staging_dir/sample_benchmark_ppv_distribs.giab_hard.external_srwgs.uris.list ) ],
                                                     [ $( collapse_txt $staging_dir/sample_benchmark_ppv_distribs.giab_easy.external_lrwgs.uris.list ),
                                                       $( collapse_txt $staging_dir/sample_benchmark_ppv_distribs.giab_hard.external_lrwgs.uris.list ) ]],
  "PlotVcfQcMetrics.sample_benchmark_sensitivity_distribs": [[ $( collapse_txt $staging_dir/sample_benchmark_sensitivity_distribs.giab_easy.external_srwgs.uris.list ),
                                                               $( collapse_txt $staging_dir/sample_benchmark_sensitivity_distribs.giab_hard.external_srwgs.uris.list ) ],
                                                             [ $( collapse_txt $staging_dir/sample_benchmark_sensitivity_distribs.giab_easy.external_lrwgs.uris.list ),
                                                               $( collapse_txt $staging_dir/sample_benchmark_sensitivity_distribs.giab_hard.external_lrwgs.uris.list ) ]],
  "PlotVcfQcMetrics.sample_genotype_distribution_tsvs": $( collapse_txt $staging_dir/genotype_distrib.uris.list ),
  "PlotVcfQcMetrics.site_benchmark_common_sv_ppv_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_sv_ppv_beds.giab_easy.uris.list ),
                                                           $( collapse_txt $staging_dir/site_benchmark_common_sv_ppv_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_common_sv_sens_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_sv_sens_beds.giab_easy.uris.list ),
                                                            $( collapse_txt $staging_dir/site_benchmark_common_sv_sens_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_ppv_by_freqs": [[ $( collapse_txt $staging_dir/site_benchmark_ppv_by_freqs.giab_easy.uris.list ),
                                                     $( collapse_txt $staging_dir/site_benchmark_ppv_by_freqs.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_sensitivity_by_freqs": [[ $( collapse_txt $staging_dir/site_benchmark_sensitivity_by_freqs.giab_easy.uris.list ),
                                                             $( collapse_txt $staging_dir/site_benchmark_sensitivity_by_freqs.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_dataset_prefixes": ["gnomad-sv_v4.1"],
  "PlotVcfQcMetrics.site_benchmark_dataset_titles": ["gnomAD-SV v4.1"],
  "PlotVcfQcMetrics.size_distribution_tsvs": $( collapse_txt $staging_dir/size_distrib.uris.list ),
  "PlotVcfQcMetrics.size_vs_af_distribution_tsvs": $( collapse_txt $staging_dir/size_vs_af_distrib.uris.list ),
  "PlotVcfQcMetrics.trio_mendelian_violation_distribs": [ $( collapse_txt $staging_dir/trio_mendelian_violation_distribs.giab_easy.uris.list ),
                                                          $( collapse_txt $staging_dir/trio_mendelian_violation_distribs.giab_hard.uris.list ) ],
  "PlotVcfQcMetrics.twin_genotype_benchmark_distribs": [ $( collapse_txt $staging_dir/twin_genotype_benchmark_distribs.giab_easy.uris.list ),
                                                         $( collapse_txt $staging_dir/twin_genotype_benchmark_distribs.giab_hard.uris.list ) ]
}
EOF

  # Add previous stats to input .json if optioned
  if ! [ -z $prev_stats ]; then
    echo -e "{\"PlotVcfQcMetrics.previous_stats\": \"$prev_stats\"}" \
    > $staging_dir/prev_stats.json
    code/scripts/update_json.py \
      -i $staging_dir/Plot$wid_suffix.inputs.json \
      -u $staging_dir/prev_stats.json \
      -o cromshell/inputs/Plot$wid_suffix.inputs.json
  else
    cp \
      $staging_dir/Plot$wid_suffix.inputs.json \
      cromshell/inputs/Plot$wid_suffix.inputs.json
  fi
}


############################################
# SUBSET G2C GATK-SV OUTPUT TO UFC SAMPLES #
############################################

# Make staging directory
staging_dir=staging/SubsetSamples
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Get list of all sample IDs present in G2C SV callset after SV imputation
gsutil -m cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup/chr22/ConcatVcfs/dfci-g2c.v1.chr22.imputed.vcf.gz \
| gunzip -c | head -n5000 | bcftools query -l \
> $staging_dir/all_g2c_samples.post_gatksv.samples.list

# Map Noah UFC IDs onto G2C sample IDs
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
  gs://rw-migration-aou-rw-cf3e7774/UFC_REFERENCE_FILES/cohorts/ufc_subjects.list \
  $staging_dir/
zcat $staging_dir/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
| awk -v FS="\t" '{ if ($NF=="True" && $3 ~ /aou|ceph|mesa|ufc/) print }' \
| sort -k2,2 \
| join -1 2 -2 1 -t $'\t' \
  - <( sort $staging_dir/ufc_subjects.list ) \
| sort -Vk1,1 -k2,2V -k3,3V \
> $staging_dir/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.ufc_samples.preheader.tsv
paste \
  <( cut -f2 $staging_dir/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.ufc_samples.preheader.tsv ) \
  <( cut --complement -f2 $staging_dir/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.ufc_samples.preheader.tsv ) \
| cat <( zcat $staging_dir/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz | head -n1 ) - \
| gzip -c \
> $staging_dir/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.ufc_samples.tsv.gz
zcat $staging_dir/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.ufc_samples.tsv.gz \
| sed '1d' | cut -f1 | sort -Vk1,1 \
> $staging_dir/all_ufc.g2c_ids.samples.list

# Define list of G2C samples to exclude for UFC callset
zcat $staging_dir/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
| sed '1d' | cut -f1 | sort -Vk1,1 \
| fgrep -xvf $staging_dir/all_ufc.g2c_ids.samples.list \
> $staging_dir/dfci-g2c.exclude_for_ufc_callset.samples.list
gsutil -m cp \
  $staging_dir/dfci-g2c.exclude_for_ufc_callset.samples.list \
  $WORKSPACE_BUCKET/data/sample_info/

# Write template .json for sample exclusion WDL
cat << EOF > $staging_dir/ExcludeSamplesFromVcf.inputs.template.json
{
  "ExcludeSamplesFromVcf.docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "ExcludeSamplesFromVcf.exclude_samples_list": "$WORKSPACE_BUCKET/data/sample_info/dfci-g2c.exclude_for_ufc_callset.samples.list",
  "ExcludeSamplesFromVcf.outfile_prefix": "dfci-ufc.gatksv.v2",
  "ExcludeSamplesFromVcf.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/sv_gt_cleanup/\$CONTIG/ConcatVcfs/dfci-g2c.v1.\$CONTIG.imputed.vcf.gz"
}
EOF

# Exclude unneeded samples from each chromosome's SV VCF
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/ExcludeSamplesFromVcf.wdl \
  --input-json-template $staging_dir/ExcludeSamplesFromVcf.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/raw_gatksv_vcfs \
  --status-tsv cromshell/progress/dfci-ufc.v2.ExcludeSamplesFromVcf.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.v2" \
  --outer-gate 30 \
  --max-attempts 2


##############################
# COLLECT INITIAL QC METRICS #
##############################

# Reaffirm staging directory
staging_dir=staging/initial_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectVcfQcMetrics.inputs.template.json
{
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 100,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.01,
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:e4eaf92",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "ubuntu:plucky-20251001",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 10000,
  "CollectVcfQcMetrics.output_prefix": "dfci-ufc.sv.v2.initial_qc.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.gcp_machine_type": "n2d-standard-4",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 15.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.ref_build": "hg38",
  "CollectVcfQcMetrics.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "CollectVcfQcMetrics.ref_fasta_idx" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "CollectVcfQcMetrics.sample_benchmark_dataset_names": ["external_srwgs", "external_lrwgs"],
  "CollectVcfQcMetrics.sample_benchmark_id_maps": [["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"],
                                                   ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"]],
  "CollectVcfQcMetrics.sample_benchmark_vcfs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz"],
                                                ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz"]],
  "CollectVcfQcMetrics.sample_benchmark_vcf_idxs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"],
                                                    ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"]],
  "CollectVcfQcMetrics.shard_vcf": false,
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4_gatksv"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.gatksv.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$WORKSPACE_BUCKET/data/sample_info/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs_array": ["$WORKSPACE_BUCKET/data/raw_gatksv_vcfs/\$CONTIG/ExcludeSamples/dfci-ufc.gatksv.v2.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs_array": ["$WORKSPACE_BUCKET/data/raw_gatksv_vcfs/\$CONTIG/ExcludeSamples/dfci-ufc.gatksv.v2.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectVcfQcMetrics.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs \
  --name CollectInitialVcfQcMetrics \
  --status-tsv cromshell/progress/dfci-ufc.sv.v2.CollectVcfQcMetrics.initial_qc.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0.1 \
  --max-attempts 3


###########################
# PLOT INITIAL QC METRICS #
###########################

# Reaffirm staging directory
staging_dir=staging/initial_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Make file for target number of SVs per genome
echo -e "variants_per_genome.sv:median\t9124" > $staging_dir/dfci-ufc.sv.qc_targets.tsv
echo -e "variants_per_genome.all:median\t9124" >> $staging_dir/dfci-ufc.sv.qc_targets.tsv
gsutil -m cp \
  $staging_dir/dfci-ufc.sv.qc_targets.tsv \
  $WORKSPACE_BUCKET/data/misc/

# Reset inputs
reset_qc_plot_dir \
  $staging_dir \
  InitialVcfQcMetrics \
  $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs

# Submit QC visualization workflow
cromshell --no_turtle -t 120 -mc submit --no-validation --do-not-flatten-wdls \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotInitialVcfQcMetrics.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-ufc.sv.v2.PlotInitialVcfQcMetrics.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotInitialVcfQcMetrics.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs/PlotQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotInitialVcfQcMetrics.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs/PlotQc/

# Clear Cromwell execution & output buckets
gsutil -m ls $( cat cromshell/job_ids/dfci-ufc.sv.v2.PlotInitialVcfQcMetrics.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/workflows/cromwell-*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


#########################################
# APPLY HIGH SENSITIVITY AOU GT FILTERS #
#########################################

# Clear cromwell cache
~/code/scripts/rotate_cromwell_cache.sh delete

# Make staging directory
staging_dir=staging/OnlyFilterGenotypes
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template .json for module 19
cat << EOF > $staging_dir/OnlyFilterGenotypes.inputs.template.json
{
  "FilterGenotypes.FilterVcf.script": "$MAIN_WORKSPACE_BUCKET/code/scripts/apply_sl_filter.g2c_imputed_integrated.py",
  "FilterGenotypes.no_call_rate_cutoff": 1,
  "FilterGenotypes.output_prefix": "dfci-ufc.v2.\$CONTIG",
  "FilterGenotypes.ploidy_table": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/17/dfci-g2c.v1.ploidy.tsv",
  "FilterGenotypes.sl_filter_args": "--small-del-threshold -53 --medium-del-threshold -12 --small-dup-threshold -105 --medium-dup-threshold -81 --ins-threshold -97",
  "FilterGenotypes.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "FilterGenotypes.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-06-27-v1.0.4-63e6c81e",
  "FilterGenotypes.vcf": "$WORKSPACE_BUCKET/data/raw_gatksv_vcfs/\$CONTIG/ExcludeSamples/dfci-ufc.gatksv.v2.vcf.gz"
}
EOF

# Run filtering component of GATK-SV module 19 on all UFC SV VCFs
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/OnlyFilterGenotypes.wdl \
  --input-json-template $staging_dir/OnlyFilterGenotypes.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/OnlyFilterGenotypes \
  --dependencies-zip gatksv.dependencies.zip \
  --status-tsv cromshell/progress/OnlyFilterGenotypes.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.v2" \
  --outer-gate 30 \
  --submission-gate 0.1 \
  --max-attempts 2


############################################
# COLLECT QC METRICS AFTER FILTERGENOTYPES #
############################################

# Clear cromwell cache
~/code/scripts/rotate_cromwell_cache.sh delete

# Reaffirm staging directory
staging_dir=staging/post_filtergenotypes
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectVcfQcMetrics.postFilterGenotypes.inputs.template.json
{
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 100,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.01,
  "CollectVcfQcMetrics.extra_vcf_preprocessing_commands": " | bcftools view -f .,PASS,MULTIALLELIC --no-update --no-version ",
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:e4eaf92",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "ubuntu:plucky-20251001",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 10000,
  "CollectVcfQcMetrics.output_prefix": "dfci-ufc.sv.v2.post_filtergenotypes.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.gcp_machine_type": "n2d-standard-4",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 15.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.ref_build": "hg38",
  "CollectVcfQcMetrics.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "CollectVcfQcMetrics.ref_fasta_idx" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "CollectVcfQcMetrics.sample_benchmark_dataset_names": ["external_srwgs", "external_lrwgs"],
  "CollectVcfQcMetrics.sample_benchmark_id_maps": [["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"],
                                                   ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"]],
  "CollectVcfQcMetrics.sample_benchmark_vcfs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz"],
                                                ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz"]],
  "CollectVcfQcMetrics.sample_benchmark_vcf_idxs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"],
                                                    ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"]],
  "CollectVcfQcMetrics.shard_vcf": false,
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4_gatksv"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.gatksv.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$WORKSPACE_BUCKET/data/sample_info/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs_array": ["$WORKSPACE_BUCKET/data/OnlyFilterGenotypes/\$CONTIG/SanitizeHeader/dfci-ufc.v2.\$CONTIG.filter_genotypes.sanitized.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs_array": ["$WORKSPACE_BUCKET/data/OnlyFilterGenotypes/\$CONTIG/SanitizeHeader/dfci-ufc.v2.\$CONTIG.filter_genotypes.sanitized.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectVcfQcMetrics.postFilterGenotypes.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $WORKSPACE_BUCKET/qc/post_filtergenotypes \
  --name CollectQcPostFilterGenotypes \
  --status-tsv cromshell/progress/dfci-ufc.sv.v2.CollectVcfQcMetrics.post_filtergenotypes.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0.1 \
  --max-attempts 2


#########################################
# PLOT QC METRICS AFTER FILTERGENOTYPES #
#########################################

# Reaffirm staging directory
staging_dir=staging/post_filtergenotypes
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Reset inputs
reset_qc_plot_dir \
  $staging_dir \
  QcPostFilterGenotypes \
  $WORKSPACE_BUCKET/qc/post_filtergenotypes \
  $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs/PlotQc/dfci-ufc.sv.v2.InitialVcfQcMetrics.all_qc_summary_metrics.tsv

# Submit QC visualization workflow
cromshell --no_turtle -t 120 -mc submit --no-validation --do-not-flatten-wdls \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotQcPostFilterGenotypes.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostFilterGenotypes.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostFilterGenotypes.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $WORKSPACE_BUCKET/qc/post_filtergenotypes/PlotQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostFilterGenotypes.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $WORKSPACE_BUCKET/qc/post_filtergenotypes/PlotQc/

# Clear Cromwell execution & output buckets
gsutil -m ls $( cat cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostFilterGenotypes.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/workflows/cromwell-*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


###############################
# COLLECT MINGQ TRAINING DATA #
###############################

# Clear cromwell cache
~/code/scripts/rotate_cromwell_cache.sh delete

# Make staging directory
staging_dir=staging/MinGQ
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template .json for minGQ part 1
cat << EOF > $staging_dir/MinGQPart1.inputs.template.json
{
  "Module07FilterGTsPart1.CollectTrioSVdat_PCRMINUS.script": "$MAIN_WORKSPACE_BUCKET/code/scripts/gather_trio_genos.g2c_imputed.py",
  "Module07FilterGTsPart1.CollectTrioSVdat_PCRPLUS.script": "$MAIN_WORKSPACE_BUCKET/code/scripts/gather_trio_genos.g2c_imputed.py",
  "Module07FilterGTsPart1.contiglist": "gs://dfci-g2c-refs/hg38/contig_fais/\$CONTIG.fai",
  "Module07FilterGTsPart1.filter_metric": "GQ",
  "Module07FilterGTsPart1.gather_trio_geno_options": ["--fill-incomplete", "--default-value-homref 1", "--default-value-other 0"],
  "Module07FilterGTsPart1.gcloud_sdk_docker": "google/cloud-sdk:581.0.0",
  "Module07FilterGTsPart1.getAFs.famfile": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "Module07FilterGTsPart1.max_shards_per_chrom_step1": 100,
  "Module07FilterGTsPart1.min_records_per_shard_step1": 50000,
  "Module07FilterGTsPart1.prefix": "dfci-ufc.sv.v2.\$CONTIG",
  "Module07FilterGTsPart1.revise_mei_svtypes": true,
  "Module07FilterGTsPart1.runtime_attr_GatherTrioData": {"disk_gb" : 250, "boot_disk_gb" : 20, "cpu_cores" : 4, "mem_gb" : 15},
  "Module07FilterGTsPart1.sv_base_mini_docker": "vanallenlab/sv-pipeline-base:gnomad_rf_99fb8c",
  "Module07FilterGTsPart1.sv_pipeline_base_docker": "vanallenlab/sv-pipeline-base:gnomad_rf_99fb8c",
  "Module07FilterGTsPart1.sv_pipeline_docker": "vanallenlab/sv-pipeline:2022-07-26-v0.23-beta-662e9fdd",
  "Module07FilterGTsPart1.sv_pipeline_docker_for_mei_revise": "vanallenlab/sv-pipeline:filter_gt_029fa3b",
  "Module07FilterGTsPart1.sv_pipeline_updates_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2022-09-27-v0.26.2-beta-1f1a3b3a",
  "Module07FilterGTsPart1.trios_famfile": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "Module07FilterGTsPart1.vcf": "$WORKSPACE_BUCKET/data/OnlyFilterGenotypes/\$CONTIG/SanitizeHeader/dfci-ufc.v2.\$CONTIG.filter_genotypes.sanitized.vcf.gz",
  "Module07FilterGTsPart1.vcf_idx": "$WORKSPACE_BUCKET/data/OnlyFilterGenotypes/\$CONTIG/SanitizeHeader/dfci-ufc.v2.\$CONTIG.filter_genotypes.sanitized.vcf.gz.tbi"
}
EOF

# Collect minGQ optimization data for all chromosomes
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/legacy_mingq_wdl/Module07FilterGTsPart1CollectData.wdl \
  --input-json-template $staging_dir/MinGQPart1.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/MinGQPart1 \
  --dependencies-zip mingq.dependencies.zip \
  --name MinGQPart1 \
  --status-tsv cromshell/progress/MinGQPart1.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --max-attempts 3


#####################
# TRAIN MINGQ MODEL #
#####################

# Reaffirm staging directory
staging_dir=staging/MinGQ
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write list of training data tarballs from all contigs
gsutil -m ls \
  $WORKSPACE_BUCKET/data/MinGQPart1/chr*/GatherTrioData_PCRMINUS/dfci-ufc.sv.v2.chr*.PCRMINUS.tar.gz \
| sort -V \
> $staging_dir/MinGQ.training_tarball.uris.list

# Write .json of inputs (run this once across all chromosomes)
cat << EOF > cromshell/inputs/MinGQPart2.inputs.json
{
  "Module07FilterGTsPart2.PCRMINUS_cleaned_trios_famfile": "$WORKSPACE_BUCKET/data/MinGQPart1/chr1/SplitFamfile_PCRMINUS/dfci-ufc.sv.v2.chr1.PCRMINUS.cleaned_trios.fam",
  "Module07FilterGTsPart2.PCRMINUS_trio_tarballs": $( collapse_txt $staging_dir/MinGQ.training_tarball.uris.list ),
  "Module07FilterGTsPart2.gcloud_sdk_docker": "google/cloud-sdk",
  "Module07FilterGTsPart2.optimize_excludeEV": ";;RD,SR;;",
  "Module07FilterGTsPart2.optimize_excludeFILTERs": ";PESR_GT_OVERDISPERSION;PESR_GT_OVERDISPERSION,HIGH_SR_BACKGROUND;;",
  "Module07FilterGTsPart2.optimize_includeEV": "RD;SR;;;",
  "Module07FilterGTsPart2.optimize_includeFILTERs": "PESR_GT_OVERDISPERSION;HIGH_SR_BACKGROUND;;;",
  "Module07FilterGTsPart2.optimize_includeSVTYPEs": "DEL;DUP;INS;MEI;INV,CPX,CTX;BND;DEL,DUP,INS,MEI,INV,CPX,CTX,BND",
  "Module07FilterGTsPart2.optimize_maxFreqs": "0.001;0.01;0.1;1;1",
  "Module07FilterGTsPart2.optimize_maxSizes": "250;1000;5000;25000;400000000;400000000",
  "Module07FilterGTsPart2.optimize_metric": "GQ",
  "Module07FilterGTsPart2.optimize_minFreqs": "0;0.001;0.01;0.1;0",
  "Module07FilterGTsPart2.optimize_minSizes": "0;250;1000;5000;25000;0",
  "Module07FilterGTsPart2.prefix": "dfci-ufc.sv.v2",
  "Module07FilterGTsPart2.roc_max_fdr_PCRMINUS": 0.05,
  "Module07FilterGTsPart2.roc_max_fdr_PCRPLUS": 0.05,
  "Module07FilterGTsPart2.roc_max_metric": 100,
  "Module07FilterGTsPart2.roc_min_metric": 0,
  "Module07FilterGTsPart2.roc_shards": 500,
  "Module07FilterGTsPart2.roc_step_metric": 1,
  "Module07FilterGTsPart2.runtime_attr_roc_single": {"disk_gb" : 50},
  "Module07FilterGTsPart2.sv_base_mini_docker": "vanallenlab/sv-base-mini:gnomad_rf_4760ac",
  "Module07FilterGTsPart2.sv_pipeline_base_docker": "vanallenlab/sv-pipeline-base:gnomad_rf_043721",
  "Module07FilterGTsPart2.sv_pipeline_base_docker_buildTree": "vanallenlab/sv-pipeline-base:gnomad_rf_abd6c1",
  "Module07FilterGTsPart2.sv_pipeline_docker": "vanallenlab/sv-pipeline:2022-07-26-v0.23-beta-662e9fdd"
}
EOF

# Launch minGQ training workflow
cromshell --no_turtle -t 120 -mc submit --no-validation --do-not-flatten-wdls \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip mingq.dependencies.zip \
  code/wdl/gatk-sv/legacy_mingq_wdl/Module07FilterGTsPart2TrainModel.wdl \
  cromshell/inputs/MinGQPart2.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-ufc.sv.v2.MinGQPart2.job_ids.list

# Monitor minGQ training workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.MinGQPart2.job_ids.list ) 5

# Once workflow is complete, stage output
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.MinGQPart2.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I $WORKSPACE_BUCKET/data/MinGQPart2/

# Clear Cromwell execution & output buckets
gsutil -m ls $( cat cromshell/job_ids/dfci-ufc.sv.v2.MinGQPart2.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/workflows/cromwell-*/Module07FilterGTsPart2/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


#####################
# APPLY MINGQ MODEL #
#####################

# Reaffirm staging directory
staging_dir=staging/MinGQ
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template .json for minGQ part 3
cat << EOF > $staging_dir/MinGQPart3.inputs.template.json
{
  "Module07FilterGTsPart3.apply_filter_PCRMINUS.script": "$MAIN_WORKSPACE_BUCKET/code/scripts/gatksv_filter_GTs_by_metric.g2c_imputed.py",
  "Module07FilterGTsPart3.apply_filter_PCRPLUS.script": "$MAIN_WORKSPACE_BUCKET/code/scripts/gatksv_filter_GTs_by_metric.g2c_imputed.py",
  "Module07FilterGTsPart3.CombineVcfs.generate_index": true,
  "Module07FilterGTsPart3.PCRMINUS_lookup_table": "$WORKSPACE_BUCKET/data/MinGQPart2/dfci-ufc.sv.v2.PCRMINUS.minGQ.filter_lookup_table.txt",
  "Module07FilterGTsPart3.PCRMINUS_vcf_idx_lists": \$CONTIG_VCF_IDXS,
  "Module07FilterGTsPart3.PCRMINUS_vcf_lists": \$CONTIG_VCFS,
  "Module07FilterGTsPart3.allow_overlaps_merge": true,
  "Module07FilterGTsPart3.filter_GT_options": ["--fail-missing-scores", "--filter-homalt", "--annotate-ncr", "--dropEmpties"],
  "Module07FilterGTsPart3.max_noCallRate": \$CONTIG_NCR,
  "Module07FilterGTsPart3.naive_merge": false,
  "Module07FilterGTsPart3.prefix": "dfci-ufc.sv.v2",
  "Module07FilterGTsPart3.runtime_attr_CombineVcfs": {"mem_gb" : 7.5, "boot_disk_gb" : 20},
  "Module07FilterGTsPart3.sort_after_merge": false,
  "Module07FilterGTsPart3.sv_base_mini_docker": "vanallenlab/sv-base-mini:gnomad_rf_4760ac",
  "Module07FilterGTsPart3.sv_pipeline_base_docker": "vanallenlab/sv-pipeline-base:gnomad_rf_99fb8c",
  "Module07FilterGTsPart3.sv_pipeline_docker": "vanallenlab/sv-pipeline:2022-07-26-v0.23-beta-662e9fdd"
}
EOF

# Build chromosome-specific override json of VCFs and VCF indexes
for k in $( seq 1 22 ) X Y; do
  echo "chr$k"
done > contig_lists/hg38.primary.contigs.list
echo "{ " > $staging_dir/MinGQPart3.contig_variable_overrides.json
while read contig; do
  if [ $contig == "chrX" ]; then
    ncr=0.4376095 # This is equivalent to all males plus 10% of females
  elif [ $contig == "chrY" ]; then
    ncr=0.6614273 # This is equivalent to all females plus 10% of males
  else
    ncr=0.1
  fi
  echo "\"$contig\" : { \"CONTIG_NCR\" : $ncr },"
done < contig_lists/hg38.primary.contigs.list \
| paste -s -d\  | sed 's/,$//g' \
>> $staging_dir/MinGQPart3.contig_variable_overrides.json
echo " }" >> $staging_dir/MinGQPart3.contig_variable_overrides.json
add_contig_vcfs_to_chromshard_overrides_json \
  $staging_dir/MinGQPart3.contig_variable_overrides.json \
  $WORKSPACE_BUCKET/data/MinGQPart1 \
  Module07FilterGTsPart1.PCRMINUS_vcfs \
  Module07FilterGTsPart1.PCRMINUS_vcf_idxs \
  contig_lists/hg38.primary.contigs.list

# Apply trained minGQ model to each chromosome
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/legacy_mingq_wdl/Module07FilterGTsPart3ApplyModel.wdl \
  --input-json-template $staging_dir/MinGQPart3.inputs.template.json \
  --contig-variable-overrides $staging_dir/MinGQPart3.contig_variable_overrides.json \
  --staging-bucket $WORKSPACE_BUCKET/data/MinGQPart3 \
  --dependencies-zip mingq.dependencies.zip \
  --name MinGQPart3 \
  --status-tsv cromshell/progress/MinGQPart3.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --max-attempts 3 \
  --contig-gate 25


##################################
# COLLECT QC METRICS AFTER MINGQ #
##################################

# Reaffirm staging directory
staging_dir=staging/post_mingq
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectVcfQcMetrics.postMinGQ.inputs.template.json
{
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 100,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.01,
  "CollectVcfQcMetrics.extra_vcf_preprocessing_commands": " | bcftools view -f .,PASS,MULTIALLELIC --no-update --no-version ",
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:1aac84d",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "ubuntu:plucky-20251001",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 10000,
  "CollectVcfQcMetrics.output_prefix": "dfci-ufc.sv.v2.post_mingq.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.gcp_machine_type": "n2d-standard-4",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 15.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.ref_build": "hg38",
  "CollectVcfQcMetrics.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "CollectVcfQcMetrics.ref_fasta_idx" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "CollectVcfQcMetrics.sample_benchmark_dataset_names": ["external_srwgs", "external_lrwgs"],
  "CollectVcfQcMetrics.sample_benchmark_id_maps": [["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"],
                                                   ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"]],
  "CollectVcfQcMetrics.sample_benchmark_vcfs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz"],
                                                ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz"]],
  "CollectVcfQcMetrics.sample_benchmark_vcf_idxs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"],
                                                    ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"]],
  "CollectVcfQcMetrics.shard_vcf": false,
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4_gatksv"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.gatksv.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$WORKSPACE_BUCKET/data/sample_info/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs_array": ["$WORKSPACE_BUCKET/data/MinGQPart3/\$CONTIG/CombineVcfs/dfci-ufc.sv.v2.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs_array": ["$WORKSPACE_BUCKET/data/MinGQPart3/\$CONTIG/CombineVcfs/dfci-ufc.sv.v2.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectVcfQcMetrics.postMinGQ.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $WORKSPACE_BUCKET/qc/post_mingq \
  --name CollectQcPostMinGQ \
  --status-tsv cromshell/progress/dfci-ufc.sv.v2.CollectVcfQcMetrics.post_mingq.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0.1 \
  --max-attempts 3


###############################
# PLOT QC METRICS AFTER MINGQ #
###############################

# Reaffirm staging directory
staging_dir=staging/post_mingq
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Reset inputs
reset_qc_plot_dir \
  $staging_dir \
  QcPostMinGQ \
  $WORKSPACE_BUCKET/qc/post_mingq \
  $WORKSPACE_BUCKET/qc/post_filtergenotypes/PlotQc/dfci-ufc.sv.v2.QcPostFilterGenotypes.all_qc_summary_metrics.tsv

# Submit QC visualization workflow
cromshell --no_turtle -t 120 -mc submit --no-validation --do-not-flatten-wdls \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotQcPostMinGQ.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostMinGQ.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostMinGQ.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $WORKSPACE_BUCKET/qc/post_mingq/PlotQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostMinGQ.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $WORKSPACE_BUCKET/qc/post_mingq/PlotQc/

# Clear Cromwell execution & output buckets
gsutil -m ls $( cat cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostMinGQ.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/workflows/cromwell-*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


##########################
# POSTHOC CLEANUP PART 1 #
##########################

# Reaffirm staging directory
staging_dir=staging/posthoc_part1
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write flat list of files needed for posthoc cleanup part 1
cat << EOF > $staging_dir/script_files.txt
$WORKSPACE_BUCKET/refs/MANE.GRCh38.v1.2.ensembl_genomic.gtf.gz
EOF

# Write template .json for posthoc cleanup part 1
cat << EOF > $staging_dir/PosthocCleanupPart1.inputs.template.json
{
  "ApplyScriptSingleVcf.ApplyScript.contig": "\$CONTIG",
  "ApplyScriptSingleVcf.bcftools_docker": "vanallenlab/g2c_analysis:e4eaf92",
  "ApplyScriptSingleVcf.exec_prefix": "python ",
  "ApplyScriptSingleVcf.script": "$WORKSPACE_BUCKET/code/scripts/vcf_cleanup.part_1.before_outlier_exclusion.py",
  "ApplyScriptSingleVcf.script_files": $( collapse_txt $staging_dir/script_files.txt ),
  "ApplyScriptSingleVcf.script_options": "--gtf MANE.GRCh38.v1.2.ensembl_genomic.gtf.gz",
  "ApplyScriptSingleVcf.vcf": "$WORKSPACE_BUCKET/data/MinGQPart3/\$CONTIG/CombineVcfs/dfci-ufc.sv.v2.vcf.gz",
  "ApplyScriptSingleVcf.vcf_idx": "$WORKSPACE_BUCKET/data/MinGQPart3/\$CONTIG/CombineVcfs/dfci-ufc.sv.v2.vcf.gz.tbi"
}
EOF

# Apply posthoc cleanup step 1 to each chromosome
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/ApplyScriptSingleVcf.wdl \
  --input-json-template $staging_dir/PosthocCleanupPart1.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/PosthocCleanupPart1 \
  --dependencies-zip g2c.dependencies.zip \
  --name PosthocCleanupPart1 \
  --status-tsv cromshell/progress/PosthocCleanupPart1.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0 \
  --max-attempts 3


###########################################
# COLLECT QC METRICS AFTER POSTHOC PART 1 #
###########################################

# Reaffirm staging directory
staging_dir=staging/post_cleanup_part1
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectVcfQcMetrics.postCleanupPart1.inputs.template.json
{
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 100,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.01,
  "CollectVcfQcMetrics.extra_vcf_preprocessing_commands": " | bcftools view -f .,PASS,MULTIALLELIC --no-update --no-version ",
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:e4eaf92",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "ubuntu:plucky-20251001",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 10000,
  "CollectVcfQcMetrics.output_prefix": "dfci-ufc.sv.v2.posthoc_part1.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.gcp_machine_type": "n2d-standard-4",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 15.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.ref_build": "hg38",
  "CollectVcfQcMetrics.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "CollectVcfQcMetrics.ref_fasta_idx" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "CollectVcfQcMetrics.sample_benchmark_dataset_names": ["external_srwgs", "external_lrwgs"],
  "CollectVcfQcMetrics.sample_benchmark_id_maps": [["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"],
                                                   ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"]],
  "CollectVcfQcMetrics.sample_benchmark_vcfs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz"],
                                                ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz"]],
  "CollectVcfQcMetrics.sample_benchmark_vcf_idxs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"],
                                                    ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"]],
  "CollectVcfQcMetrics.shard_vcf": false,
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4_gatksv"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.gatksv.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$WORKSPACE_BUCKET/data/sample_info/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs_array": ["$WORKSPACE_BUCKET/data/PosthocCleanupPart1/\$CONTIG/ApplyScript/dfci-ufc.sv.v2.\$CONTIG.out.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs_array": ["$WORKSPACE_BUCKET/data/PosthocCleanupPart1/\$CONTIG/ApplyScript/dfci-ufc.sv.v2.\$CONTIG.out.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectVcfQcMetrics.postCleanupPart1.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $WORKSPACE_BUCKET/qc/post_cleanupPart1 \
  --name CollectQcPostCleanupPart1 \
  --status-tsv cromshell/progress/dfci-ufc.sv.v2.CollectVcfQcMetrics.post_cleanupPart1.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0.1 \
  --max-attempts 3


########################################
# PLOT QC METRICS AFTER CLEANUP PART 1 #
########################################

# Reaffirm staging directory
staging_dir=staging/post_cleanup_part1
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Reset inputs
reset_qc_plot_dir \
  $staging_dir \
  QcPostCleanupPart1 \
  $WORKSPACE_BUCKET/qc/post_cleanupPart1 \
  $WORKSPACE_BUCKET/qc/post_mingq/PlotQc/dfci-ufc.sv.v2.QcPostMinGQ.all_qc_summary_metrics.tsv

# Submit QC visualization workflow
cromshell --no_turtle -t 120 -mc submit --no-validation --do-not-flatten-wdls \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotQcPostCleanupPart1.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostCleanupPart1.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostCleanupPart1.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $WORKSPACE_BUCKET/qc/post_cleanupPart1/PlotQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostCleanupPart1.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $WORKSPACE_BUCKET/qc/post_cleanupPart1/PlotQc/

# Clear Cromwell execution & output buckets
gsutil -m ls $( cat cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostCleanupPart1.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/workflows/cromwell-*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


######################
# MARK BATCH EFFECTS #
######################

# Reaffirm staging directory
staging_dir=staging/batch_effects
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write GATK-SV batch membership file for European ancestry samples
gsutil -m cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
| gunzip -c | awk '{ if ($6=="EUR") print $1 }' \
> $staging_dir/g2c.EUR.samples.list
gsutil -m cat \
  $WORKSPACE_BUCKET/data/PosthocCleanupPart1/chr22/ApplyScript/dfci-ufc.sv.v2.chr22.out.vcf.gz \
| bcftools query -l \
> $staging_dir/ufc.samples.list
while read bid; do
  gsutil cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/batch_info/sample_lists/$bid.samples.list \
  | fgrep -xf $staging_dir/g2c.EUR.samples.list \
  | fgrep -xf $staging_dir/ufc.samples.list \
  | awk -v OFS="\t" -v bid="$bid" '{ print bid, $1 }'
done < <( gsutil cat $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/batch_info/dfci-g2c.gatk-sv.batches.list ) \
| sort -Vk1,1 -k2,2V \
> $staging_dir/dfci-ufc.v2.sv.EUR.batch_assignments.tsv
gsutil -m cp \
  $staging_dir/dfci-ufc.v2.sv.EUR.batch_assignments.tsv \
  $WORKSPACE_BUCKET/data/misc/

# Write template .json for batch effect detection
cat << EOF > $staging_dir/MarkBatchEffects.inputs.template.json
{
  "MarkBatchEffects.custom_info_flag_id": "BATCH_EFFECTS_MASKED",
  "MarkBatchEffects.flag_failing_records": false,
  "MarkBatchEffects.g2c_pipeline_docker": "vanallenlab/g2c_analysis:e4eaf92",
  "MarkBatchEffects.group_membership_tsv": "$WORKSPACE_BUCKET/data/misc/dfci-ufc.v2.sv.EUR.batch_assignments.tsv",
  "MarkBatchEffects.lower_nonref_gt_freq": 0.02,
  "MarkBatchEffects.min_samples_per_group": 20,
  "MarkBatchEffects.out_vcf_prefix": "dfci-g2c.sv.v2.\$CONTIG",
  "MarkBatchEffects.ref_fai": "gs://dfci-g2c-refs/hg38/contig_fais/\$CONTIG.fai",
  "MarkBatchEffects.upper_nonref_gt_freq": 0.20,
  "MarkBatchEffects.vcf": "$WORKSPACE_BUCKET/data/PosthocCleanupPart1/\$CONTIG/ApplyScript/dfci-ufc.sv.v2.\$CONTIG.out.vcf.gz",
  "MarkBatchEffects.vcf_idx": "$WORKSPACE_BUCKET/data/PosthocCleanupPart1/\$CONTIG/ApplyScript/dfci-ufc.sv.v2.\$CONTIG.out.vcf.gz.tbi"
}
EOF

# Run & manage batch effect detection
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/MarkBatchEffects.wdl \
  --input-json-template $staging_dir/MarkBatchEffects.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/MarkBatchEffects \
  --dependencies-zip g2c.dependencies.zip \
  --status-tsv cromshell/progress/MarkBatchEffects.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0 \
  --max-attempts 2


#######################
# MARK COHORT EFFECTS #
#######################

# Reaffirm staging directory
staging_dir=staging/cohort_effects
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write GATK-SV cohort membership file for European ancestry samples
gsutil -m cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
| gunzip -c | awk -v FS="\t" -v OFS="\t" '{ if ($6=="EUR") print $3, $1 }' \
| fgrep -wf staging/batch_effects/ufc.samples.list \
> $staging_dir/dfci-ufc.v2.sv.EUR.cohort_assignments.tsv
gsutil -m cp \
  $staging_dir/dfci-ufc.v2.sv.EUR.cohort_assignments.tsv \
  $WORKSPACE_BUCKET/data/misc/

# Write template .json for cohort effect detection
cat << EOF > $staging_dir/MarkCohortEffects.inputs.template.json
{
  "MarkBatchEffects.custom_info_flag_id": "COHORT_EFFECTS_MASKED",
  "MarkBatchEffects.flag_failing_records": false,
  "MarkBatchEffects.g2c_pipeline_docker": "vanallenlab/g2c_analysis:e4eaf92",
  "MarkBatchEffects.group_membership_tsv": "$WORKSPACE_BUCKET/data/misc/dfci-ufc.v2.sv.EUR.cohort_assignments.tsv",
  "MarkBatchEffects.lower_nonref_gt_freq": 0.02,
  "MarkBatchEffects.min_samples_per_group": 20,
  "MarkBatchEffects.out_vcf_prefix": "dfci-ufc.sv.v2.\$CONTIG",
  "MarkBatchEffects.ref_fai": "gs://dfci-g2c-refs/hg38/contig_fais/\$CONTIG.fai",
  "MarkBatchEffects.upper_nonref_gt_freq": 0.20,
  "MarkBatchEffects.vcf": "$WORKSPACE_BUCKET/data/MarkBatchEffects/\$CONTIG/ConcatVcfs/dfci-g2c.sv.v2.\$CONTIG.vcf.gz",
  "MarkBatchEffects.vcf_idx": "$WORKSPACE_BUCKET/data/MarkBatchEffects/\$CONTIG/ConcatVcfs/dfci-g2c.sv.v2.\$CONTIG.vcf.gz.tbi"
}
EOF

# Run & manage cohort effect detection
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/MarkBatchEffects.wdl \
  --input-json-template $staging_dir/MarkCohortEffects.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/MarkCohortEffects \
  --dependencies-zip g2c.dependencies.zip \
  --status-tsv cromshell/progress/MarkCohortEffects.progress.tsv \
  --name MarkCohortEffects \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0 \
  --max-attempts 2


############################
# OUTLIER SAMPLE EXCLUSION #
############################

# Reaffirm staging directory
staging_dir=staging/outlier_exclusion
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write input .json for SV counting task
for file in $staging_dir/vcfs.list $staging_dir/vcf_idxs.list; do
  if [ -e $file ]; then rm $file; fi
done
for k in $( seq 1 22 ) X Y; do
  echo "$WORKSPACE_BUCKET/data/MarkCohortEffects/chr$k/ConcatVcfs/dfci-ufc.sv.v2.chr$k.vcf.gz" \
  >> $staging_dir/vcfs.list
  echo "$WORKSPACE_BUCKET/data/MarkCohortEffects/chr$k/ConcatVcfs/dfci-ufc.sv.v2.chr$k.vcf.gz.tbi" \
  >> $staging_dir/vcf_idxs.list
done
cat << EOF > cromshell/inputs/count_svs_posthoc.inputs.json
{
  "CountSvsPerSample.bcftools_view_options": "-f PASS,MULTIALLELIC",
  "CountSvsPerSample.g2c_pipeline_docker": "vanallenlab/g2c_analysis:e4eaf92",
  "CountSvsPerSample.output_prefix": "dfci-ufc.v2.postCleanupPart1",
  "CountSvsPerSample.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-01-14-v1.0.1-88dbd052",
  "CountSvsPerSample.vcfs": $( collapse_txt $staging_dir/vcfs.list ),
  "CountSvsPerSample.vcf_idxs": $( collapse_txt $staging_dir/vcf_idxs.list )
}
EOF

# Submit SV counting task
cromshell --no_turtle -t 120 -mc submit --no-validation --do-not-flatten-wdls \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip g2c.dependencies.zip \
  code/wdl/gatk-sv/CountSvsPerSample.wdl \
  cromshell/inputs/count_svs_posthoc.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/count_svs_posthoc.job_ids.list

# Monitor SV counting task
monitor_workflow $( tail -n1 cromshell/job_ids/count_svs_posthoc.job_ids.list )

# Once complete, download SV counts per sample, collapse CPX+INV, and exclude CTX and mCNVs
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/count_svs_posthoc.job_ids.list ) \
| awk '{ print $NF }' | gsutil cp -I $staging_dir/
sed 's/\tINV\t\|\tCPX\t/\tINV_CPX\t/g' \
  $staging_dir/dfci-ufc.v2.postCleanupPart1.counts.tsv \
| awk -v FS="\t" -v OFS="\t" '{ if ($2 !~ /CTX|BND|CNV/) print }' \
> $staging_dir/dfci-ufc.v2.postCleanupPart1.counts.subsetted.tsv
code/scripts/sum_svcounts.py \
  --outfile $staging_dir/dfci-ufc.v2.postCleanupPart1.counts.subsetted.collapsed.tsv \
  $staging_dir/dfci-ufc.v2.postCleanupPart1.counts.subsetted.tsv

# Regenerate ancestry labels split by cohort
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $staging_dir/
pop_idx=$( zcat $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
           | sed -n '1p' | sed 's/\t/\n/g' \
           | awk '{ if ($1=="intake_qc_pop") print NR }' )
zcat $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
| sed '1d' \
| awk -v idx=$pop_idx -v FS="\t" -v OFS="\t" '{ if ($3!="aou") $3="oth"; print $1, $idx"_"$3 }' \
| cat <( echo -e "sample_id\tlabel" ) - \
> $staging_dir/dfci-g2c.intake_pop_labels.aou_split.tsv

# Define outliers
code/scripts/define_variant_count_outlier_samples.R \
  --counts-tsv $staging_dir/dfci-ufc.v2.postCleanupPart1.counts.subsetted.collapsed.tsv \
  --sample-labels-tsv $staging_dir/dfci-g2c.intake_pop_labels.aou_split.tsv \
  --n-iqr 3 \
  --plot \
  --plot-title-prefix "GATKSV" \
  --out-prefix $staging_dir/dfci-ufc.sv.v2.posthoc_outliers

# Compress and archive outlier data for future reference
cd $staging_dir && \
tar -czvf dfci-ufc.sv.v2.posthoc_outliers.tar.gz dfci-ufc.sv.v2.posthoc_outliers* && \
gsutil -m cp \
  dfci-ufc.sv.v2.posthoc_outliers.tar.gz \
  dfci-ufc.sv.v2.posthoc_outliers.outliers.samples.list \
  $WORKSPACE_BUCKET/data/posthoc_outliers/ && \
cd ~

# Write template input .json for outlier exclusion & hard filter task
cat << EOF > $staging_dir/OutlierExclusion.inputs.template.json
{
  "PosthocHardFilterPart2.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PosthocHardFilterPart2.exclude_samples_list": "$WORKSPACE_BUCKET/data/posthoc_outliers/dfci-ufc.sv.v2.posthoc_outliers.outliers.samples.list",
  "PosthocHardFilterPart2.vcf": "$WORKSPACE_BUCKET/data/MarkCohortEffects/\$CONTIG/ConcatVcfs/dfci-ufc.sv.v2.\$CONTIG.vcf.gz",
  "PosthocHardFilterPart2.vcf_idx": "$WORKSPACE_BUCKET/data/MarkCohortEffects/\$CONTIG/ConcatVcfs/dfci-ufc.sv.v2.\$CONTIG.vcf.gz.tbi"
}
EOF

# Submit outlier exclusion
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/PosthocHardFilterPart2.wdl \
  --input-json-template $staging_dir/OutlierExclusion.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/OutlierExclusionWorkflow \
  --name OutlierExclusion \
  --status-tsv cromshell/progress/dfci-ufc.v2.OutlierExclusion.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0 \
  --max-attempts 3


######################################################
# COLLECT QC AFTER BATCH/COHORT EFFECTS AND OUTLIERS #
######################################################

# Reaffirm staging directory
staging_dir=staging/post_bfx_outliers
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectVcfQcMetrics.postBatchFxOutliers.inputs.template.json
{
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 100,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.01,
  "CollectVcfQcMetrics.extra_vcf_preprocessing_commands": " | bcftools view -f .,PASS,MULTIALLELIC --no-update --no-version ",
  "CollectVcfQcMetrics.g2c_analysis_docker": "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/g2c_analysis:ec08fc1",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/ubuntu:plucky-20251001",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 10000,
  "CollectVcfQcMetrics.output_prefix": "dfci-ufc.sv.v2.post_bfx_outliers.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 15.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.ref_build": "hg38",
  "CollectVcfQcMetrics.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "CollectVcfQcMetrics.ref_fasta_idx" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "CollectVcfQcMetrics.sample_benchmark_dataset_names": ["external_srwgs", "external_lrwgs"],
  "CollectVcfQcMetrics.sample_benchmark_id_maps": [["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"],
                                                   ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"]],
  "CollectVcfQcMetrics.sample_benchmark_vcfs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz"],
                                                ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz"]],
  "CollectVcfQcMetrics.sample_benchmark_vcf_idxs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"],
                                                    ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"]],
  "CollectVcfQcMetrics.shard_vcf": false,
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4_gatksv"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.gatksv.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$WORKSPACE_BUCKET/data/sample_info/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs_array": ["$WORKSPACE_BUCKET/data/OutlierExclusionWorkflow/\$CONTIG/HardFilterPart2/dfci-ufc.sv.v2.\$CONTIG.posthoc_filtered.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs_array": ["$WORKSPACE_BUCKET/data/OutlierExclusionWorkflow/\$CONTIG/HardFilterPart2/dfci-ufc.sv.v2.\$CONTIG.posthoc_filtered.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectVcfQcMetrics.postBatchFxOutliers.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $WORKSPACE_BUCKET/qc/postBatchFxOutliers \
  --name CollectQcPostBatchFxOutliers \
  --status-tsv cromshell/progress/dfci-ufc.sv.v2.CollectVcfQcMetrics.post_bfx_outliers.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0.1 \
  --max-attempts 4


###################################################
# PLOT QC AFTER BATCH/COHORT EFFECTS AND OUTLIERS #
###################################################

# Reaffirm staging directory
staging_dir=staging/post_bfx_outliers
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Reset inputs
reset_qc_plot_dir \
  $staging_dir \
  QcPostBatchFxOutliers \
  $WORKSPACE_BUCKET/qc/postBatchFxOutliers \
  $WORKSPACE_BUCKET/qc/post_cleanupPart1/PlotQc/dfci-ufc.sv.v2.QcPostCleanupPart1.all_qc_summary_metrics.tsv

# Submit QC visualization workflow
cromshell --no_turtle -t 120 -mc submit --no-validation --do-not-flatten-wdls \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotQcPostBatchFxOutliers.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostBatchFxOutliers.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostBatchFxOutliers.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $WORKSPACE_BUCKET/qc/post_batchFxOutliers/PlotQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostBatchFxOutliers.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $WORKSPACE_BUCKET/qc/post_batchFxOutliers/PlotQc/

# Clear Cromwell execution & output buckets
gsutil -m ls $( cat cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostBatchFxOutliers.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/workflows/cromwell-*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


# ##################
# # VISUALIZE CNVS #
# ##################

# # TODO: implement this (maybe not necessary?)


##########################
# POSTHOC CLEANUP PART 2 #
##########################

# Reaffirm staging directory
staging_dir=staging/posthoc_part2
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write flat list of files needed for posthoc cleanup part 2
cat << EOF > $staging_dir/script_files.txt
$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped
EOF

# Write template .json for posthoc cleanup part 2
cat << EOF > $staging_dir/PosthocCleanupPart2.inputs.template.json
{
  "ApplyScriptSingleVcf.ApplyScript.contig": "\$CONTIG",
  "ApplyScriptSingleVcf.bcftools_docker": "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/g2c_analysis:1aac84d",
  "ApplyScriptSingleVcf.exec_prefix": "python ",
  "ApplyScriptSingleVcf.script": "$WORKSPACE_BUCKET/code/scripts/vcf_cleanup.part_2.after_outlier_exclusion.py",
  "ApplyScriptSingleVcf.script_files": $( collapse_txt $staging_dir/script_files.txt ),
  "ApplyScriptSingleVcf.script_options": "--ped-file dfci-g2c.all_samples.ped --version-number 2.0",
  "ApplyScriptSingleVcf.vcf": "$WORKSPACE_BUCKET/data/OutlierExclusionWorkflow/\$CONTIG/HardFilterPart2/dfci-ufc.sv.v2.\$CONTIG.posthoc_filtered.vcf.gz",
  "ApplyScriptSingleVcf.vcf_idx": "$WORKSPACE_BUCKET/data/OutlierExclusionWorkflow/\$CONTIG/HardFilterPart2/dfci-ufc.sv.v2.\$CONTIG.posthoc_filtered.vcf.gz.tbi"
}
EOF

# Apply posthoc cleanup step 2 to each chromosome
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/ApplyScriptSingleVcf.wdl \
  --input-json-template $staging_dir/PosthocCleanupPart2.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/PosthocCleanupPart2 \
  --dependencies-zip g2c.dependencies.zip \
  --name PosthocCleanupPart2 \
  --status-tsv cromshell/progress/PosthocCleanupPart2.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0 \
  --max-attempts 2


###########################################
# COLLECT QC METRICS AFTER POSTHOC PART 2 #
###########################################

# Reaffirm staging directory
staging_dir=staging/post_cleanup_part2
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectVcfQcMetrics.postCleanupPart2.inputs.template.json
{
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 100,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.01,
  "CollectVcfQcMetrics.extra_vcf_preprocessing_commands": " | bcftools view -f .,PASS,MULTIALLELIC --no-update --no-version ",
  "CollectVcfQcMetrics.g2c_analysis_docker": "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/vanallenlab/g2c_analysis:ec08fc1",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/ubuntu:plucky-20251001",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 10000,
  "CollectVcfQcMetrics.output_prefix": "dfci-ufc.sv.v2.post_mingq.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.gcp_machine_type": "n2d-standard-4",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 15.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.ref_build": "hg38",
  "CollectVcfQcMetrics.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "CollectVcfQcMetrics.ref_fasta_idx" : "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "CollectVcfQcMetrics.sample_benchmark_dataset_names": ["external_srwgs", "external_lrwgs"],
  "CollectVcfQcMetrics.sample_benchmark_id_maps": [["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"],
                                                   ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.1KGP_id_map.tsv",
                                                    "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/dfci-g2c.v1.AoU_id_map.tsv"]],
  "CollectVcfQcMetrics.sample_benchmark_vcfs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz"],
                                                ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz",
                                                 "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz"]],
  "CollectVcfQcMetrics.sample_benchmark_vcf_idxs": [["gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/sv/1KGP.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/srwgs/sv/AoU.srWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"],
                                                    ["gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi",
                                                     "$MAIN_WORKSPACE_BUCKET/refs/aou/dense_vcfs/lrwgs/sv/AoU.lrWGS.sv.cleaned.\$CONTIG.vcf.gz.tbi"]],
  "CollectVcfQcMetrics.shard_vcf": false,
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4_gatksv"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.gatksv.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$WORKSPACE_BUCKET/data/sample_info/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs_array": ["$WORKSPACE_BUCKET/data/PosthocCleanupPart2/\$CONTIG/ApplyScript/dfci-ufc.sv.v2.\$CONTIG.posthoc_filtered.\$CONTIG.out.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs_array": ["$WORKSPACE_BUCKET/data/PosthocCleanupPart2/\$CONTIG/ApplyScript/dfci-ufc.sv.v2.\$CONTIG.posthoc_filtered.\$CONTIG.out.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectVcfQcMetrics.postCleanupPart2.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $WORKSPACE_BUCKET/qc/post_cleanupPart2 \
  --name CollectQcPostCleanupPart2 \
  --status-tsv cromshell/progress/dfci-ufc.sv.v2.CollectVcfQcMetrics.post_cleanupPart2.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0.1 \
  --max-attempts 3


########################################
# PLOT QC METRICS AFTER POSTHOC PART 2 #
########################################

# Reaffirm staging directory
staging_dir=staging/post_cleanup_part2
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Reset inputs
reset_qc_plot_dir \
  $staging_dir \
  QcPostCleanupPart2 \
  $WORKSPACE_BUCKET/qc/post_cleanupPart2 \
  $WORKSPACE_BUCKET/qc/post_batchFxOutliers/PlotQc/dfci-ufc.sv.v2.QcPostBatchFxOutliers.all_qc_summary_metrics.tsv

# Submit QC visualization workflow
cromshell --no_turtle -t 120 -mc submit --no-validation --do-not-flatten-wdls \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotQcPostCleanupPart2.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostCleanupPart2.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostCleanupPart2.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $WORKSPACE_BUCKET/qc/post_cleanupPart2/PlotQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostCleanupPart2.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $WORKSPACE_BUCKET/qc/post_cleanupPart2/PlotQc/

# Clear Cromwell execution & output buckets
gsutil -m ls $( cat cromshell/job_ids/dfci-ufc.sv.v2.PlotQcPostCleanupPart2.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/workflows/cromwell-*/PlotVcfQcMetrics/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage


################
# ANNOTATE VCF #
################

# Reaffirm staging directory
staging_dir=staging/annotation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template .json for callset annotation
cat << EOF > $staging_dir/AnnotateVcf.inputs.template.json
{
  "AnnotateVcf.contig_list": "gs://dfci-g2c-refs/hg38/contig_lists/\$CONTIG.list",
  "AnnotateVcf.external_af_population": ["ALL","AFR","AMR","EAS","EUR","MID","FIN","ASJ","RMI","SAS","AMI"],
  "AnnotateVcf.external_af_ref_bed": "gs://gatk-sv-resources-public/gnomad_AF/gnomad_v4_SV.Freq.tsv.gz",
  "AnnotateVcf.external_af_ref_prefix": "gnomad_v4.1_sv",
  "AnnotateVcf.gatk_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT",
  "AnnotateVcf.par_bed": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.par.bed",
  "AnnotateVcf.prefix": "dfci-ufc.sv.v2.1.\$CONTIG",
  "AnnotateVcf.protein_coding_gtf": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/gencode.v47.basic.protein_coding.canonical.gtf",
  "AnnotateVcf.runtime_attr_svannotate": {"cpu_cores" : 4, "mem_gb": 7.5},
  "AnnotateVcf.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "AnnotateVcf.sv_per_shard": 5000,
  "AnnotateVcf.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-09-02-v1.0.5-631368eb",
  "AnnotateVcf.vcf": "$WORKSPACE_BUCKET/data/PosthocCleanupPart2/\$CONTIG/ApplyScript/dfci-ufc.sv.v2.\$CONTIG.posthoc_filtered.\$CONTIG.out.vcf.gz"
}
EOF

# Annotate each chromosome
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/svannotation_v1.1/AnnotateVcf.wdl \
  --input-json-template $staging_dir/AnnotateVcf.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/AnnotateVcf_v2 \
  --dependencies-zip svannotation.dependencies.zip \
  --status-tsv cromshell/progress/AnnotateVcf.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v2" \
  --outer-gate 30 \
  --submission-gate 0 \
  --max-attempts 2

