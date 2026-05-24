# Diverse Mediators of Cancer Predisposition Uncovered by Germline Whole Genome Sequencing of Unexplained Familial Cancers

## Loss-of-Function (LoF) SV analysis


## Synopsis    

This subrepository contains the working code and scripts used to detect, genotype, filter, and annotate germline variants from germline WGS across cancer types. This code also contains code to analyze and graph results.

A preprint of this study is available on medRxiv: https://www.medrxiv.org/content/10.64898/2026.05.08.26352653v1
--- 

## LoF SV script Functions
 1. dfci-ufc.aou.prepare_fig2A - We match rare (AF <1%) LoF SVs to cancer types 
 2. dfci-ufc.aou.sv_fishers_step1 - We create a list of rare, ultra-rare (AF <0.1%), and singleton LoF SVs and match to all samples
 3. dfci-ufc.aou.sv_fishers_step2 - We use the results of step1 to analyze all 18 cancer cohorts for enrichments of rare, ultra-rare, and singleton SVs.