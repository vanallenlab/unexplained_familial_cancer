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


############################################
# SUBSET G2C GATK-SV OUTPUT TO UFC SAMPLES #
############################################

# TODO: implement this


#########################################
# APPLY HIGH SENSITIVITY AOU GT FILTERS #
#########################################

# TODO: implement this


#########################
# TRAIN AND APPLY MINGQ #
#########################

# TODO: implement this


#############
# OTHER TBD #
#############

# TODO: implement this

