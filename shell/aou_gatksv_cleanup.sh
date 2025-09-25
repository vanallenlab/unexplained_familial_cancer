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
find code/ -name "*.sh" | xargs -I {} chmod a+x {}

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

# Create dependencies .zip for generic G2C workflow submissions
cd code/wdl/pancan_germline_wgs && \
zip -r g2c.dependencies.zip . && \
mv g2c.dependencies.zip ~/ && \
cd ~

# Create dependencies .zip for QC workflow submissions
cd code/wdl/pancan_germline_wgs/vcf-qc && \
zip qc.dependencies.zip *.wdl && \
mv qc.dependencies.zip ~/ && \
cd ~

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

# Install necessary packages
. code/refs/install_packages.sh python R

# Download workspace-specific contig lists
gsutil cp -r \
  gs://dfci-g2c-refs/hg38/contig_lists \
  ./


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
cat << EOF | python -m json.tool > cromshell/inputs/Plot$wid_suffix.inputs.json
{
  "PlotVcfQcMetrics.af_distribution_tsvs": $( collapse_txt $staging_dir/af_distrib.uris.list ),
  "PlotVcfQcMetrics.all_sv_beds": $( collapse_txt $staging_dir/all_svs_bed.uris.list ),
  "PlotVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PlotVcfQcMetrics.benchmark_interval_names": ["Easy", "Hard"],
  "PlotVcfQcMetrics.common_af_cutoff": 0.01,
  "PlotVcfQcMetrics.common_sv_beds": $( collapse_txt $staging_dir/common_svs_bed.uris.list ),
  "PlotVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:7d85acf",
  "PlotVcfQcMetrics.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "PlotVcfQcMetrics.output_prefix": "dfci-ufc.sv.v1.$wid_suffix",
  "PlotVcfQcMetrics.peak_ld_stat_tsvs": $( collapse_txt $staging_dir/ld_stats.uris.list ),
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
  "PlotVcfQcMetrics.site_benchmark_dataset_prefixes": ["gnomad_v4.1_gatksv"],
  "PlotVcfQcMetrics.site_benchmark_dataset_titles": ["gnomAD-SV v4.1"],
  "PlotVcfQcMetrics.size_distribution_tsvs": $( collapse_txt $staging_dir/size_distrib.uris.list ),
  "PlotVcfQcMetrics.size_vs_af_distribution_tsvs": $( collapse_txt $staging_dir/size_vs_af_distrib.uris.list ),
  "PlotVcfQcMetrics.trio_mendelian_violation_distribs": [ $( collapse_txt $staging_dir/trio_mendelian_violation_distribs.giab_easy.uris.list ),
                                                          $( collapse_txt $staging_dir/trio_mendelian_violation_distribs.giab_hard.uris.list ) ],
  "PlotVcfQcMetrics.twin_genotype_benchmark_distribs": [ $( collapse_txt $staging_dir/twin_genotype_benchmark_distribs.giab_easy.uris.list ),
                                                         $( collapse_txt $staging_dir/twin_genotype_benchmark_distribs.giab_hard.uris.list ) ]
}
EOF
}


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
  "ExcludeSamplesFromVcf.vcf": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/CollapseRedundantSvs/\$CONTIG/RC3/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.identical.reclustered.vcf.gz"
}
EOF

# Exclude unneeded samples from each chromosome's SV VCF
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/ExcludeSamplesFromVcf.wdl \
  --input-json-template $staging_dir/ExcludeSamplesFromVcf.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/raw_gatksv_vcfs \
  --status-tsv cromshell/progress/dfci-ufc.v1.ExcludeSamplesFromVcf.progress.tsv \
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
  "CollectVcfQcMetrics.all_samples_fam_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped",
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": 100,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.01,
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:75e54bf",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 10000,
  "CollectVcfQcMetrics.output_prefix": "dfci-ufc.sv.v1.initial_qc.\$CONTIG",
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
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4_gatksv"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.gatksv.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$WORKSPACE_BUCKET/data/sample_info/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs": ["$WORKSPACE_BUCKET/data/raw_gatksv_vcfs/\$CONTIG/ExcludeSamples/dfci-ufc.gatksv.v1.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs": ["$WORKSPACE_BUCKET/data/raw_gatksv_vcfs/\$CONTIG/ExcludeSamples/dfci-ufc.gatksv.v1.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectVcfQcMetrics.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs \
  --name CollectInitialVcfQcMetrics \
  --status-tsv cromshell/progress/dfci-ufc.sv.v1.CollectVcfQcMetrics.initial_qc.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v1" \
  --outer-gate 30 \
  --submission-gate 0.1 \
  --max-attempts 3


###########################
# PLOT INITIAL QC METRICS #
###########################

# Reaffirm staging directory
staging_dir=staging/initial_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Reset inputs
reset_qc_plot_dir \
  $staging_dir \
  InitialVcfQcMetrics \
  $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs

# Submit QC visualization workflow
cromshell --no_turtle -t 120 -mc submit --no-validation \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotInitialVcfQcMetrics.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-ufc.sv.v1.PlotInitialVcfQcMetrics.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v1.PlotInitialVcfQcMetrics.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs/PlotQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v1.PlotInitialVcfQcMetrics.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs/PlotQc/

# TODO: UPDATE THIS
# # Clear Cromwell execution & output buckets for patch jobs
# gsutil -m ls $( cat cromshell/job_ids/GnarlyJointGenotypingPart1.inputs.$contig.patch.job_ids.list \
#                 | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwell/*/GnarlyJointGenotypingPart1/" \
#                   '{ print bucket_prefix$1"/**" }' ) \
# > uris_to_delete.list
# cleanup_garbage


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
  "FilterGenotypes.vcf": "$WORKSPACE_BUCKET/data/raw_gatksv_vcfs/\$CONTIG/ExcludeSamples/dfci-ufc.gatksv.v1.vcf.gz"
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
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:75e54bf",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 4568,
  "CollectVcfQcMetrics.output_prefix": "dfci-ufc.sv.v1.post_filtergenotypes.\$CONTIG",
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
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4_gatksv"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": [],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": [],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.gatksv.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$WORKSPACE_BUCKET/data/sample_info/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.twins_tsv": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.tsv",
  "CollectVcfQcMetrics.vcfs": ["$WORKSPACE_BUCKET/data/module19/\$CONTIG/SanitizeHeader/dfci-ufc.v1.\$CONTIG.filter_genotypes.sanitized.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs": ["$WORKSPACE_BUCKET/data/module19/\$CONTIG/SanitizeHeader/dfci-ufc.v1.\$CONTIG.filter_genotypes.sanitized.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectVcfQcMetrics.postFilterGenotypes.inputs.template.json \
  --dependencies-zip qc.dependencies.zip \
  --staging-bucket $WORKSPACE_BUCKET/qc/post_filtergenotypes \
  --name CollectQcPostFilterGenotypes \
  --status-tsv cromshell/progress/dfci-ufc.sv.v1.CollectVcfQcMetrics.post_filtergenotypes.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.sv.v1" \
  --outer-gate 30 \
  --submission-gate 0.1 \
  --max-attempts 3


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
  $WORKSPACE_BUCKET/qc/post_filtergenotypes

# Submit QC visualization workflow
cromshell --no_turtle -t 120 -mc submit --no-validation \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip qc.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotQcPostFilterGenotypes.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-ufc.sv.v1.PlotQcPostFilterGenotypes.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v1.PlotQcPostFilterGenotypes.job_ids.list ) 5

# Once workflow is complete, stage output
gsutil -m rm -rf $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs/PlotQc
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-ufc.sv.v1.PlotQcPostFilterGenotypes.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $WORKSPACE_BUCKET/qc/raw_gatksv_vcfs/PlotQc/

# TODO: UPDATE THIS
# # Clear Cromwell execution & output buckets for patch jobs
# gsutil -m ls $( cat cromshell/job_ids/GnarlyJointGenotypingPart1.inputs.$contig.patch.job_ids.list \
#                 | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwell/*/GnarlyJointGenotypingPart1/" \
#                   '{ print bucket_prefix$1"/**" }' ) \
# > uris_to_delete.list
# cleanup_garbage


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
  "Module07FilterGTsPart1.prefix": "dfci-ufc.sv.v1.\$CONTIG",
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
  --workflow-id-log-prefix "dfci-ufc.sv.v1" \
  --outer-gate 30 \
  --max-attempts 2


#####################
# TRAIN MINGQ MODEL #
#####################

# Reaffirm staging directory
staging_dir=staging/MinGQ
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write list of training data tarballs from all contigs
gsutil -m ls \
  $WORKSPACE_BUCKET/data/MinGQPart1/chr*/GatherTrioData_PCRMINUS/dfci-ufc.v1.chr*.PCRMINUS.tar.gz \
| sort -V \
> $staging_dir/MinGQ.training_tarball.uris.list

# Write .json of inputs (run this once across all chromosomes)
cat << EOF > cromshell/inputs/MinGQPart2.inputs.json
{
  "Module07FilterGTsPart2.PCRMINUS_cleaned_trios_famfile": "$WORKSPACE_BUCKET/data/MinGQPart1/chr1/SplitFamfile_PCRMINUS/dfci-ufc.v1.chr1.PCRMINUS.cleaned_trios.fam",
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
  "Module07FilterGTsPart2.prefix": "dfci-ufc.v1",
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
cromshell --no_turtle -t 120 -mc submit --no-validation \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip mingq.dependencies.zip \
  code/wdl/gatk-sv/legacy_mingq_wdl/Module07FilterGTsPart2TrainModel.wdl \
  cromshell/inputs/MinGQPart2.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-ufc.v1.MinGQPart2.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-ufc.v1.MinGQPart2.job_ids.list )

# Once workflow is complete, stage output
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-ufc.v1.MinGQPart2.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I $WORKSPACE_BUCKET/data/MinGQPart2/

# # Clear Cromwell execution & output buckets
# gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v1.InferTwins.job_ids.list \
#                 | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwell-*/InferTwins/" \
#                   '{ print bucket_prefix$1"/**" }' ) \
# > uris_to_delete.list
# cleanup_garbage


#####################
# APPLY MINGQ MODEL #
#####################

# Reaffirm staging directory
staging_dir=staging/MinGQ
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template .json for minGQ part 3
cat << EOF > $staging_dir/MinGQPart3.inputs.template.json
{
  "Module07FilterGTsPart3.CombineVcfs.generate_index": true,
  "Module07FilterGTsPart3.PCRMINUS_lookup_table": "${this.minGQ_PCRMINUS_lookup_table_postGQR_lenient}",
  "Module07FilterGTsPart3.PCRMINUS_vcf_idx_lists": "${this.minGQ_PCRMINUS_vcf_idx_shards_postGQR_lenient}",
  "Module07FilterGTsPart3.PCRMINUS_vcf_lists": "${this.minGQ_PCRMINUS_vcf_shards_postGQR_lenient}",
  "Module07FilterGTsPart3.allow_overlaps_merge": true,
  "Module07FilterGTsPart3.filter_GT_options": ["--fail-missing-scores", "--filter-homalt", "--annotate-ncr", "--dropEmpties"],
  "Module07FilterGTsPart3.max_noCallRate": 0.05,
  "Module07FilterGTsPart3.naive_merge": false,
  "Module07FilterGTsPart3.prefix": "dfci-ufc.v1",
  "Module07FilterGTsPart3.runtime_attr_CombineVcfs": {"mem_gb" : 7.5, "boot_disk_gb" : 20},
  "Module07FilterGTsPart3.sort_after_merge": false,
  "Module07FilterGTsPart3.sv_base_mini_docker": "vanallenlab/sv-base-mini:gnomad_rf_4760ac",
  "Module07FilterGTsPart3.sv_pipeline_base_docker": "vanallenlab/sv-pipeline-base:gnomad_rf_99fb8c",
  "Module07FilterGTsPart3.sv_pipeline_docker": "vanallenlab/sv-pipeline:2022-07-26-v0.23-beta-662e9fdd"
}
EOF

# Apply trained minGQ model to each chromosome
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-sv/legacy_mingq_wdl/Module07FilterGTsPart3ApplyModel.wdl \
  --input-json-template $staging_dir/MinGQPart3.inputs.template.json \
  --staging-bucket $WORKSPACE_BUCKET/data/MinGQPart3 \
  --dependencies-zip mingq.dependencies.zip \
  --name MinGQPart3 \
  --status-tsv cromshell/progress/MinGQPart3.progress.tsv \
  --workflow-id-log-prefix "dfci-ufc.v1" \
  --outer-gate 30 \
  --max-attempts 2


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

