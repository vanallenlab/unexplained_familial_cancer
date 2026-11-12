# The Unexplained Familial Cancer (UFC) Project

## WDLs used for UFC project

Copyright (c) 2023-Present, Noah Fields , Ryan L. Collins, and the Van Allen laboratory at Dana-Farber Cancer Institute.
Distributed under terms of the GNU GPL v2.0 License (see LICENSE).

**Note: this repository is under active development. More documentation will be added as the project evolves.**

---  

 ### DOCKER PREFIX
 1. us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/

This repository contains the working code and scripts used to process genotype information.

--- 
 ### Steps to recreate a cohort while including MBD4, CTNNA1, MSH3, RPS20, HOXB13, and GREM1 (154 total) 
 1. Create Cohort 1 - Using the AoU CDRv9 Exome Hail MT, extract CPGs w/ a 2kb buffer
 2. Do work to 1) VEP annotate these files, and 2) curate a genotype only list of VCFs
 3. Create Cohort 3 - Obtain a list of individuals who have a PGV in any of 154 PGVs Output - "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/"

 ### Steps Involved in Replication (R)
 1. Step R0 - An auxillary WDL to take vcfs with unannotated IDs and annotate them with 'CHROM_POS_REF_ALT'.
 2. Step R1 - Given a list of Genes, output files listing who carries PGVs in which genes. (Example: TP53.tier1_0001.list holds all individuals with a Tier 1 variant (AF < 0.1%) in TP53); Output - "gs://dataproc-staging-wb-cordial-diamond-9893/step_r1_output/"
 3. Step R3 - XX
 4. Step R4 - Calculate polygenic risk X scores on Y sample IDs. Output - "gs://dataproc-staging-wb-cordial-diamond-9893/step_r4_output/"
 5. Step R5 - Residualize PCs on AoU generated PCs. Useful for PRS distribution (130d41a7-b164-499a-b116-1508a8c252e1)
 6. Step R6 - Get the carriers of nominated PGVs (using same AF AND Tier)
 
 ### BRAT1 Analysis
 1. Location of BRAT1 region VCF and index - 1MB buffer (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/BRAT1_region.vcf.bgz)
 2. File containing variants of interest for AoU PCs (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/loadings.ht.tar.gz)
 3. BRAT1 LD Results (Chi-Squared) - (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/brat1_haplotype_ld.tsv.gz)
 
 ### WDL Dependency Packages
 1. Ufc_utilities.zip (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/Ufc_utlities.zip)
 2. Vep_utilities.zip (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/Vep_utlities.zip)
 
 ### PRS files
 1. RAW PRS for everyone in AoU CDRv9 (gs://dataproc-staging-wb-cordial-diamond-9893/AoU_Cancer_PRS/CDRv9_prs.tsv.gz)
 2. Ancestry Adjusted PRS in AoU CDRv9 (gs://dataproc-staging-wb-cordial-diamond-9893/AoU_Cancer_PRS/CDRv9_prs.5pc_adjusted.tsv.gz)
 
 ### Key Files For Replication
 1. 154 genes denoted as CPGs for this analysis (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/ufc_cpg.list)
 2. 154 genes denoted as CPGs for this analysis (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/ufc_cpg.bed)
 3. Metadata for CDRv7 (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/CDRv7.metadata)
 4. Metadata for UFC CDRv7 (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/CDRv7.ufc_original.metadata)
 5. Metadata for CDRv9 (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/CDRv9.metadata)
 6. Metadata for Replication UFC CDRv9 Replication Cohort (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/CDRv9.ufc_replication.metadata)
 7. Dataset containing inferred_sex, intake_qc_pop (Ryan Curated) (gs://rw-migration-aou-rw-84a0039b/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz)
 8. AoU generated PCs 1 - 10 for everybody (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/AoU_CDRv9_pcs.tsv.gz)
 9. List of Genes to replicate (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/12_genes_to_replicate.list; some genes have already been processed)
 
 ### Relatedness Info
 1. Discovery UFC Principal Components (gs://rw-migration-aou-rw-cf3e7774/STEP_11_GENETIC_RELATEDNESS/ufc.eigenvec)
 2. Discovery UFC LD pruned variants (gs://rw-migration-aou-rw-cf3e7774/STEP_11_GENETIC_RELATEDNESS/ufc_1000G_snps.ld_pruned.vcf.gz)
 3. All of Us Flagged Related Samples (gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv)
 
 
 ### VCF Folders (sharded by CPGs)
  1. CDRv9 No Genotypes, VEP Annotated (gs://dataproc-staging-wb-cordial-diamond-9893/vep_annotated_cpgs_sept2026/CDRv9/)
  2. CDRv9 No Genotypes, No annotations (gs://dataproc-staging-wb-cordial-diamond-9893/unannotated_cpgs_sept2026/)
  3. CDRv9 Genotypes, No annotations (gs://dataproc-staging-wb-cordial-diamond-9893/unannotated_genotyped_cpgs_sept2026/CDRv9/)
  4. CDRv7 UFC Genotypes, No Annotations (gs://rw-migration-aou-rw-cf3e7774/STEP_8_FILTER_TO_TP_VARIANTS/sharded_vcfs/)
  5. CDRv7 UFC Genotypes, VEP Annotated (gs://rw-migration-aou-rw-cf3e7774/STEP_9_RUN_VEP/)
  
  ### Important List of Samples
  1. AoU UFC Samples (gs://rw-migration-aou-rw-cf3e7774/UFC_REFERENCE_FILES/cohorts/ufc_subjects.aou.list)
  2. TSV of Individiduals (CDV7 - Discovery Cohort) with PGVs and what genes and variants they have (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/ufc_subjects.aou.pgv_carriers_in_154_cpgs.tsv)
  2. List of Individuals (CDRv7 - Discovery Cohort) to exclude for PGVs in 154 CPGs (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/ufc_subjects.aou.pgv_carriers_in_154_cpgs.list)
  3. List of Individuals (CDRv9) with PGVs in 154 CPGs (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/oct11_2026_pgv_carriers.CDRv9.list.tsv)
  4. List of Individuals passing SV QC (gs://rw-migration-aou-rw-cf3e7774/all_ufc_samples.march4.sv_qc_pass.2695.list)
  5. List of All Samples in v7 (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/v7.subjects)
  6. List of All Samples in v9 (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/v9.subjects)
  7. List of Samples in v9 but not v7 (n=294,373) (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/v9_only.subjects)
  8. List of Individuals w/ PGVs in BRCA1 in CDRv9 (gs://dataproc-staging-wb-cordial-diamond-9893/step_r1_output/BRCA1.tier1_0001.list; all of them)
  
  ### SAIGE FILES To SPEED THINGS UP
  1. Concatenated VCF for Rare Variants used by SAIGE-GENE+ (AF<1%) (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/ufc_rare_variants.001.vcf.gz)
  2. Concatenated VCF for Rare Variants used by SAIGE-GENE+ (AF<0.1%) (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/ufc_rare_variants.0001.vcf.gz)
  3. UFC LD Plink BED (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/ufc_ld.bed)
  4. UFC LD Plink BIM (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/ufc_ld.bim)
  5. UFC LD Plink FAM (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/ufc_ld.fam)
  6. Sparse Matrix GRM (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx)
  7. Sparse Matrix IDs (gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt)
  
  
  # Researcher Account Activity Confirmation Link
  1. https://redcap.pmi-ops.org/surveys/?s=K4PA9H8E979AHJEP

  
