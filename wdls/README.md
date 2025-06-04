# The Unexplained Familial Cancer (UFC) Project

## WDLs used for UFC project

Copyright (c) 2023-Present, Noah Fields , Ryan L. Collins, and the Van Allen laboratory at Dana-Farber Cancer Institute.
Distributed under terms of the GNU GPL v2.0 License (see LICENSE).

**Note: this repository is under active development. More documentation will be added as the project evolves.**

---  

## Synopsis    

This repository contains the working code and scripts used to process genotype information.

--- 

## Steps Involved in QC 
 1. Step1 - We filter out samples that deemed low quality in step1.
 2. Step2 - Use the GnarlyJointGenotyper to Joint Genotype
 3. Step3 - Hard Filter out variants that don't pass sufficient QC metrics.
 4. Step4 - Run a random forest model to decide assign true positive probabilities to each alternate allele.

 ### Assess the RF model that we trained and find a true positive threshold for snvs and indels
 5. Step5 - Count the total number of variants per genome and the number of de novo variants in trios where possible.
 6. Step6 - Assess the sensitivity of how many variants
 7. Step7 - Assess Hardy-Weinberg Equilibrium for variants in the cohort and unrelated european variants from AllofUs sub-cohort of the study.

 ### Create Final VCF
 8. Step8 - Based on the chosen true positive threshold for SNVs and Indels, we will filter to alternate alleles that have a true positive threshold above the chosen threshold (SNVs and Indels have their own threshold).
 9. Step9 - Run Variant Ensemble Prediction (VEP) on the VCFs.

 ### Create Summary Stats for VCF
 10. Step10 - Count the number of autosomal coding variants for each subject in the cohort.
 11. Step11 - Create PCs for the AllofUs arm of the cohort as well as the kinship coefficients.
 14. Step14 - Compare the allele frequencies for gnomAD v3.1 against in cohort allele frequencies stratified on ancestry.


 ### Filter to Subjects for Analysis
 12. Step12 - Create specific cohorts for each cancer type analysis.
 13. Step13 - Find which samples in our cohort already have a pathogenic variant.

--- 