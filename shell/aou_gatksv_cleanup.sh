#!/usr/bin/env bash

# The Van Allen Lab Unexplained Familial Cancers (UFC) Study
# Copyright (c) 2025-Present, Ryan L. Collins, Noah Fields, and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to refine GATK-SV outputs for all UFC samples

# Note that the below code makes *heavy* use of code from the DFCI G2C Database
# This code is available/documented here:
# https://github.com/vanallenlab/pancan_germline_wgs


#########
# SETUP #
#########

# Note that this setup makes heavy assumptions about standards and helper 
# functions defined for the G2C project, which we repurpose here for convenience

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in scratch staging cromshell cromshell/inputs cromshell/job_ids cromshell/progress; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}
find code/ -name "*.R" | xargs -I {} chmod a+x {}

# Source .bashrc and bash utility functions
. code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh

# Ensure Cromwell/Cromshell are configured
code/scripts/setup_cromshell.py

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Create dependencies .zip for workflow submissions
cd code/wdl/pancan_germline_wgs && \
zip g2c.dependencies.zip *.wdl && \
mv g2c.dependencies.zip ~/ && \
cd ~

# Install necessary packages
. code/refs/install_packages.sh python R

# Download workspace-specific contig lists
gsutil cp -r \
  gs://dfci-g2c-refs/hg38/contig_lists \
  ./


############################################
# SUBSET G2C GATK-SV OUTPUT TO UFC SAMPLES #
############################################

# Make staging directory
staging_dir=staging/SubsetSamples
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Get list of all sample IDs present in G2C SV callset after variant recalibration
gsutil -m cat \
  gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45/dfci-g2c-callsets/gatk-sv/module-outputs/18/chr22/ConcatVcfs/dfci-g2c.v1.chr22.concordance.vcf.gz \
| gunzip -c | head -n5000 | bcftools query -l \
> $staging_dir/all_g2c_samples.post_gatksv.samples.list

# Map Noah UFC IDs onto G2C sample IDs
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/dfci-g2c.sample_meta.posthoc_outliers.ceph_update.tsv.gz \
  gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_4_RANDOM_FOREST/inputs/ufc_subjects.list \
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
  "ExcludeSamplesFromVcf.outfile_prefix": "dfci-ufc.gatksv.v1",
  "ExcludeSamplesFromVcf.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/18/\$CONTIG/ConcatVcfs/dfci-g2c.v1.\$CONTIG.concordance.vcf.gz"
}
EOF

# Exclude unneeded samples from each chromosome's SV VCF
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/ExcludeSamplesFromVcf.wdl \
  --input-json-template $staging_dir/ExcludeSamplesFromVcf.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/raw_gatksv_vcfs \
  --status-tsv cromshell/progress/ExcludeSamplesFromVcf.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.v1" \
  --outer-gate 30 \
  --max-attempts 2


#########################################
# APPLY HIGH SENSITIVITY AOU GT FILTERS #
#########################################

# Make staging directory
staging_dir=staging/FilterGenotypes
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template .json for module 19
  cat << EOF > $staging_dir/FilterGenotypes.inputs.template.json
{
    "FilterGenotypes.gatk_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT",
    "FilterGenotypes.genome_tracks": ["gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-RepeatMasker.bed.gz",
                                      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-Segmental-Dups.bed.gz",
                                      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-Simple-Repeats.bed.gz",
                                      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38_umap_s100.bed.gz",
                                      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38_umap_s24.bed.gz"],
    "FilterGenotypes.gq_recalibrator_model_file": "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/gatk-sv-recalibrator.aou_phase_1.v1.model",
    "FilterGenotypes.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
    "FilterGenotypes.no_call_rate_cutoff": 1,
    "FilterGenotypes.output_prefix": "dfci-ufc.v1.\$CONTIG",
    "FilterGenotypes.ped_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
    "FilterGenotypes.ploidy_table": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/17/dfci-g2c.v1.ploidy.tsv",
    "FilterGenotypes.primary_contigs_fai": "gs://dfci-g2c-refs/hg38/contig_fais/\$CONTIG.fai",
    "FilterGenotypes.recalibrate_gq_args": ["--keep-homvar false","--keep-homref true","--keep-multiallelic true","--skip-genotype-filtering true","--min-samples-to-estimate-allele-frequency -1"],
    "FilterGenotypes.run_qc": false,
    "FilterGenotypes.runtime_override_plot_qc_per_family": {"mem_gb" : 15, "disk_gb" : 100},
    "FilterGenotypes.sl_filter_args": "--small-del-threshold -53 --medium-del-threshold -12 --small-dup-threshold -105 --medium-dup-threshold -81 --ins-threshold -97",
    "FilterGenotypes.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
    "FilterGenotypes.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-06-27-v1.0.4-63e6c81e",
    "FilterGenotypes.vcf": "$WORKSPACE_BUCKET/data/raw_gatksv_vcfs/\$CONTIG/ConcatVcfs/dfci-g2c.v1.\$CONTIG.concordance.vcf.gz"
}
EOF

# Run module 19 on all UFC SV VCFs
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/FilterGenotypes.wdl \
  --input-json-template $staging_dir/FilterGenotypes.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/module19 \
  --status-tsv cromshell/progress/FilterGenotypes.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.v1" \
  --outer-gate 30 \
  --max-attempts 2


#########################
# TRAIN AND APPLY MINGQ #
#########################

# TODO: implement this


#############
# OTHER TBD #
#############

# TODO: implement this

