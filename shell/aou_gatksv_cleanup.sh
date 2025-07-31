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
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/misc/legacy_mingq_wdl code/wdl/gatk-sv/
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
cp vcf-qc/*.wdl ./ && \
zip g2c.dependencies.zip *.wdl && \
mv g2c.dependencies.zip ~/ && \
cd ~

# Create dependencies .zip for GATK-SV module submissions
cd code/wdl/gatk-sv && \
zip gatksv.dependencies.zip *.wdl && \
mv gatksv.dependencies.zip ~/ && \
cd ~

# Create dependencies .zip for minGQ module submissions
cd code/wdl/gatk-sv/legacy_mingq_wdl && \
zip mingq.dependencies.zip *.wdl && \
mv mingq.dependencies.zip ~/ && \
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


##############################
# COLLECT INITIAL QC METRICS #
##############################

# Reaffirm staging directory
staging_dir=staging/initial_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectVcfQcMetrics.inputs.template.json
{
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 100,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.01,
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:10c3529",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 4568,
  "CollectVcfQcMetrics.output_prefix": "dfci-ufc.v1.initial_qc.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 15.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
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
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.snv.sites.bed.gz"],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.indel.sites.bed.gz"],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs": ["$WORKSPACE_BUCKET/data/raw_gatksv_vcfs/\$CONTIG/ExcludeSamples/DFCI_UFC_WGS.gatksv.v1.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs": ["$WORKSPACE_BUCKET/data/raw_gatksv_vcfs/\$CONTIG/ExcludeSamples/DFCI_UFC_WGS.gatksv.v1.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectVcfQcMetrics.inputs.template.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs \
  --name CollectInitialVcfQcMetrics \
  --status-tsv cromshell/progress/dfci-ufc.v1.CollectVcfQcMetrics.initial_qc.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.v1" \
  --outer-gate 30 \
  --max-attempts 3


###########################
# PLOT INITIAL QC METRICS #
###########################

# TODO: IMPLEMENT THIS


#########################################
# APPLY HIGH SENSITIVITY AOU GT FILTERS #
#########################################

# Make staging directory
staging_dir=staging/FilterGenotypes
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template .json for module 19
cat << EOF > $staging_dir/FilterGenotypes.inputs.template.json
{
  "FilterGenotypes.gatk_docker": "us.gcr.io/broad-dsde-methods/markw/gatk:mw-tb-form-sv-filter-training-data-899360a",
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
  "FilterGenotypes.vcf": "$WORKSPACE_BUCKET/data/raw_gatksv_vcfs/\$CONTIG/ExcludeSamples/DFCI_UFC_WGS.gatksv.v1.vcf.gz"
}
EOF

# Run module 19 on all UFC SV VCFs
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/FilterGenotypes.wdl \
  --input-json-template $staging_dir/FilterGenotypes.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/module19 \
  --dependencies-zip gatksv.dependencies.zip \
  --status-tsv cromshell/progress/FilterGenotypes.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.v1" \
  --outer-gate 30 \
  --max-attempts 2


############################################
# COLLECT QC METRICS AFTER FILTERGENOTYPES #
############################################

# TODO: IMPLEMENT THIS


#########################################
# PLOT QC METRICS AFTER FILTERGENOTYPES #
#########################################

# TODO: IMPLEMENT THIS


###############################
# COLLECT MINGQ TRAINING DATA #
###############################

# Make staging directory
staging_dir=staging/MinGQ
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template .json for minGQ part 1
cat << EOF > $staging_dir/MinGQPart1.inputs.template.json
{
  "Module07FilterGTsPart1.contiglist": "gs://dfci-g2c-refs/hg38/contig_fais/\$CONTIG.fai",
  "Module07FilterGTsPart1.filter_metric": "GQ",
  "Module07FilterGTsPart1.gather_trio_geno_options": ["--fill-incomplete", "--default-value-homref 1", "--default-value-other 0"],
  "Module07FilterGTsPart1.gcloud_sdk_docker": "gcr.io/google.com/cloudsdktool/google-cloud-cli@sha256:c31cf1fa2832affd0c875c89909245b9e794035472e2ec088c7c136a89df107e",
  "Module07FilterGTsPart1.getAFs.famfile": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "Module07FilterGTsPart1.max_shards_per_chrom_step1": 100,
  "Module07FilterGTsPart1.min_records_per_shard_step1": 50000,
  "Module07FilterGTsPart1.prefix": "dfci-ufc.v1.\$CONTIG",
  "Module07FilterGTsPart1.revise_mei_svtypes": true,
  "Module07FilterGTsPart1.runtime_attr_GatherTrioData": {"disk_gb" : 250, "boot_disk_gb" : 20, "cpu_cores" : 4, "mem_gb" : 15},
  "Module07FilterGTsPart1.sv_base_mini_docker": "vanallenlab/sv-pipeline-base:gnomad_rf_99fb8c",
  "Module07FilterGTsPart1.sv_pipeline_base_docker": "vanallenlab/sv-pipeline-base:gnomad_rf_99fb8c",
  "Module07FilterGTsPart1.sv_pipeline_docker": "vanallenlab/sv-pipeline:2022-07-26-v0.23-beta-662e9fdd",
  "Module07FilterGTsPart1.sv_pipeline_docker_for_mei_revise": "vanallenlab/sv-pipeline:filter_gt_029fa3b",
  "Module07FilterGTsPart1.sv_pipeline_updates_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2022-09-27-v0.26.2-beta-1f1a3b3a",
  "Module07FilterGTsPart1.trios_famfile": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "Module07FilterGTsPart1.vcf": "$WORKSPACE_BUCKET/data/module19/\$CONTIG/SanitizeHeader/dfci-ufc.v1.\$CONTIG.filter_genotypes.sanitized.vcf.gz",
  "Module07FilterGTsPart1.vcf_idx": "$WORKSPACE_BUCKET/data/module19/\$CONTIG/SanitizeHeader/dfci-ufc.v1.\$CONTIG.filter_genotypes.sanitized.vcf.gz.tbi"
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
  --workflow-id-log-prefix "dfci-ufc.v1" \
  --outer-gate 30 \
  --max-attempts 2


#####################
# TRAIN MINGQ MODEL #
#####################

# TODO: implement this


#####################
# APPLY MINGQ MODEL #
#####################

# TODO: implement this


##################################
# COLLECT QC METRICS AFTER MINGQ #
##################################

# TODO: IMPLEMENT THIS


###############################
# PLOT QC METRICS AFTER MINGQ #
###############################

# TODO: IMPLEMENT THIS


#############
# OTHER TBD #
#############

# TODO: implement this

