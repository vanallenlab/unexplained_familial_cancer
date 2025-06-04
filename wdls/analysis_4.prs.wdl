# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_4_PRS {
  input {
    String step_8_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_8_FILTER_TO_TP_VARIANTS/sharded_vcfs"
    File sample_data
    File prs_file
    File step_11_pcs = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_11_GENETIC_RELATEDNESS/ufc.eigenvec"
    File step_12_subcohort
    String PGS_ID
    String cancer_type
    String analysis_4_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/"
  }

  Int negative_shards = 0

  # Takes in a directory and outputs a Array[File] holding all of the vcf shards for each pathway
  call Tasks.gather_vcfs as gather_vcfs{
    input:
      dir = step_8_output_dir
  }

  call Tasks.sort_vcf_list {
    input:
      unsorted_vcf_list = gather_vcfs.vcf_list
  }

  call T1_get_variant_ids {
    input:
      prs_file = prs_file
  }
  scatter (i in range(length(sort_vcf_list.vcf_arr) - negative_shards)){
    call Tasks.Filter_VCF_USING_ID {
      input:
        vcf = sort_vcf_list.vcf_arr[i],
        variant_list  = T1_get_variant_ids.out1 #This is to filter by ID: CHROM_POS_REF_ALT
    }
  }
  call Tasks.ConcatVcfs {
    input:
      vcfs = Filter_VCF_USING_ID.out1,
      vcf_idxs = Filter_VCF_USING_ID.out2,
      callset_name = "place_holder" 
  }

  call T3_calculate_scores {
    input:
      vcf = ConcatVcfs.merged_vcf,
      prs_file = prs_file
  }
  call T4_control_for_ancestry {
    input:
      raw_prs = T3_calculate_scores.out1,
      sample_data = sample_data,
      step_11_pcs = step_11_pcs,
      step_12_subcohort = step_12_subcohort,
      cancer_type = cancer_type,
      PGS_ID = PGS_ID
  }
  call Tasks.copy_file_to_storage {
    input:
      text_file = T4_control_for_ancestry.out1,
      output_dir = analysis_4_output_dir
  } 
}

task T1_get_variant_ids {
  input{
    File prs_file
  }
  command <<<
  python3 <<CODE
  import pandas as pd

  df = pd.read_csv("~{prs_file}", comment = "#", index_col=False, sep='\t')
  df['ID'] = "chr" + df['hm_chr'].astype(str) + "_" + df['hm_pos'].astype(str) + "_" + df['effect_allele'] + "_" + df['other_allele']
  df['other_ID'] = "chr" + df['hm_chr'].astype(str) + "_" + df['hm_pos'].astype(str) + "_" + df['other_allele'] + "_" + df['effect_allele']
  df['ID'].to_csv('variant_ids.txt',header=False,index=False)
  df['other_ID'].to_csv('variant_ids.txt',index=False,header=False, mode='a')
  CODE
  >>>
  output {
    File out1 = "variant_ids.txt"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}

task T3_calculate_scores {
  input{
    File vcf
    File prs_file
  }
  command <<<
  set -euxo pipefail
  bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GT]\n' ~{vcf} > all_genotypes.txt

  python3 <<CODE
  import pandas as pd

  # Load PRS data
  df = pd.read_csv("~{prs_file}", sep='\t', comment='#')
  df['ID'] = "chr" + df['hm_chr'].astype(str) + "_" + df['hm_pos'].astype(str) + "_" + df['effect_allele'] + "_" + df['other_allele']
  df['other_ID'] = "chr" + df['hm_chr'].astype(str) + "_" + df['hm_pos'].astype(str) + "_" + df['other_allele'] + "_" + df['effect_allele']

  # Combine weights and allele info
  id_to_weight = dict(zip(df['ID'], df['effect_weight']))
  id_to_effect = dict(zip(df['ID'], df['effect_allele']))
  id_to_weight.update(zip(df['other_ID'], df['effect_weight']))
  id_to_effect.update(zip(df['other_ID'], df['effect_allele']))

  scores = {}
  print("Finished Preparing to Assign Scores!")

  with open("all_genotypes.txt") as f:
      for line in f:
          chrom, pos, var_id, ref, alt, *genotypes = line.strip().split('\t')

          if var_id not in id_to_weight:
              continue

          weight = id_to_weight[var_id]
          effect_allele = id_to_effect[var_id]

          for gt_field in genotypes:
              sample, gt = gt_field.split('=')
              if '.' in gt:
                  continue

              alleles = gt.replace('|', '/').split('/')
              if effect_allele == ref:
                  dosage = alleles.count('0')  # 0 = REF
              elif effect_allele == alt:
                  dosage = alleles.count('1')  # 1 = ALT
              else:
                  continue  # effect_allele doesn't match REF or ALT

              scores[sample] = round(scores.get(sample, 0.0) + dosage * weight, 3)

  print("Finished Assigning Scores")

  with open("pgs_scores.tsv", "w") as out:
      out.write("sample\tPGS\n")
      for sample, score in scores.items():
          out.write(f"{sample}\t{score}\n")
 
  CODE
  >>>
  output {
    File out1 = "pgs_scores.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}

task T4_control_for_ancestry {
  input {
    File raw_prs
    File sample_data
    File step_11_pcs
    File step_12_subcohort
    String cancer_type
    String PGS_ID 
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  import statsmodels.api as sm
  from scipy.stats import zscore

  # Load input files
  prs_df = pd.read_csv("~{raw_prs}", sep='\t', index_col=False)
  pcs_df = pd.read_csv("~{step_11_pcs}", sep='\t', index_col=False)
  sample_data = pd.read_csv("~{sample_data}", sep='\t', index_col=False)

  # Convert merge keys to string
  prs_df['sample'] = prs_df['sample'].astype(str)
  pcs_df['#IID'] = pcs_df['#IID'].astype(str)
  sample_data['original_id'] = sample_data['original_id'].astype(str)
  print(sample_data.columns)

  # Merge on sample ID
  merged_df = prs_df.merge(pcs_df, left_on = "sample", right_on="#IID")
  merged_df = merged_df.merge(sample_data, left_on = "sample", right_on="original_id")
  print(merged_df.columns)

  # Filter to subjects of interest
  subcohort_list = []
  with open("~{step_12_subcohort}", 'r') as f:
    subcohort_list = set(line.strip() for line in f)

  merged_df = merged_df[merged_df['sample'].isin(subcohort_list)]
  merged_df[['PC1','PC2','PC3','PC4','age']] = merged_df[['PC1','PC2','PC3','PC4','age']].apply(zscore)
  merged_df['male_sex'] = merged_df['inferred_sex'].apply(lambda x: 1 if x == 'male' else 0)

  # Set up regression: regress PGS ~ PC1 + PC2 + PC3 + PC4 + age + male_sex
  X = merged_df[['PC1', 'PC2', 'PC3', 'PC4','male_sex']]
  X = sm.add_constant(X)
  y = merged_df['PGS']

  model = sm.OLS(y, X).fit()
  merged_df['PGS'] = model.resid

  # Calculate mean and std using only controls
  control_pgs = merged_df.loc[merged_df['cancer'] == 'control', 'PGS']
  control_mean = control_pgs.mean()
  control_std = control_pgs.std()

  # Standardize PGS for all samples using control-derived mean and std
  merged_df['PGS'] = (merged_df['PGS'] - control_mean) / control_std

  # Write final output
  merged_df[['sample', 'PGS']].to_csv("~{cancer_type}.~{PGS_ID}.pgs", sep='\t', index=False)

  CODE
  >>>
  output {
    File out1 = "~{cancer_type}.~{PGS_ID}.pgs"
  }
  runtime {
    docker: "rblancomi/statsmodels"
    preemptible: 3
  }
}
