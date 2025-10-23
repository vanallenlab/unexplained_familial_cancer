# Unexplained Familial Cancer (UFC)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0g

version 1.0
import "Ufc_utilities/Ufc_utilities.wdl" as Tasks

workflow ANALYSIS_4B_PRS {
  input {
    Array[String] matched_cancer_prs = ["breast","breast","breast","breast","bladder","bladder","bladder","cervix","cervix","colorectal","colorectal","colorectal","colorectal","uterus","uterus","uterus","kidney","kidney","kidney","leukemia","lung","lung","lung","lung","melanoma","melanoma","melanoma","non-hodgkins","non-hodgkins","ovary","ovary","ovary","ovary","pancreas","pancreas","pancreas","prostate","prostate","prostate","prostate","thyroid","brain","esophagus"]
    Array[String] PGS_IDS = ["PGS000783","PGS003380","PGS004242","PGS004688","PGS004241","PGS000782","PGS004687","PGS000784","PGS003389","PGS000785","PGS003386","PGS004243","PGS004689","PGS000786","PGS003381","PGS004244","PGS000787","PGS004690","PGS004245","PGS000788","PGS000789","PGS003391","PGS004246","PGS004691","PGS000790","PGS003382","PGS004247","PGS000791","PGS004248","PGS000793","PGS003385","PGS004249","PGS004692","PGS000794","PGS004250","PGS004693","PGS000795","PGS003383","PGS004251","PGS004694","PGS000797","PGS003384","PGS003388"]
    Array[String] cancer_types = ["basal_cell","bladder","breast","colorectal","cervix","esophagus","gastrointestinal","hematologic","kidney","lung","melanoma","neuroendocrine","nervous","non-hodgkins","ovary","pancreas","prostate","sarcoma","squamous_cell","thyroid","uterus","leukemia"]
    String analysis_4_output_dir = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ANALYSIS_4_PRS/"
    String workspace_bucket = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"
  }
  File sample_data = "gs://~{workspace_bucket}/UFC_REFERENCE_FILES/analysis/"
  Int negative_shards = 0


  scatter (i in range(length(PGS_IDS) - negative_shards)) {
    String PGS_ID = PGS_IDS[i]
    scatter( cancer_type in cancer_types){
      File specific_cohort_sample_data = sample_data + cancer_type + "/" + cancer_type + ".metadata"
      call T4_control_for_ancestry {
        input:
          raw_prs = analysis_4_output_dir + PGS_ID + ".raw.pgs",
          sample_data = specific_cohort_sample_data,
          cancer_type = cancer_type,
          PGS_ID = PGS_ID,
          appropriate_cancer_type = matched_cancer_prs[i]
      }

      call Tasks.copy_file_to_storage as copy1{
        input:
          text_file = T4_control_for_ancestry.out1,
          output_dir = analysis_4_output_dir + "ADJUSTED_PRS/"
      }

      call Tasks.copy_file_to_storage as copy1_anon{
        input:
          text_file = T4_control_for_ancestry.out2,
          output_dir = analysis_4_output_dir + "ADJUSTED_PRS_ANON/"
      }

      call T4_control_for_ancestry as T4_control_for_ancestry_euro {
        input:
          raw_prs = analysis_4_output_dir + PGS_ID + ".raw.pgs",
          sample_data = specific_cohort_sample_data,
          cancer_type = cancer_type,
          PGS_ID = PGS_ID,
          appropriate_cancer_type = matched_cancer_prs[i],
          euro = "True"
      }

      call Tasks.copy_file_to_storage as copy2{
        input:
          text_file = T4_control_for_ancestry_euro.out1,
          output_dir = analysis_4_output_dir + "ADJUSTED_PRS/"
      }

      call T5_get_summary_statistics {
        input:
          adjusted_prs = T4_control_for_ancestry.out1,
          sample_data = specific_cohort_sample_data,
          cancer_type = cancer_type,
          PGS_ID = PGS_ID
      }
      call T5_get_summary_statistics as T5_get_summary_statistics_euro{
        input:
          adjusted_prs = T4_control_for_ancestry_euro.out1,
          sample_data = specific_cohort_sample_data,
          cancer_type = cancer_type,
          PGS_ID = PGS_ID
      }
    }
    call T6_Combine_Values as T6_Combine_Pvalues{
      input:
        files = T5_get_summary_statistics.out2_p,
        appropriate_cancer_type = matched_cancer_prs[i],
        PGS_ID = PGS_ID
    }
    call T6_Combine_Values as T6_Combine_OR{
      input:
        files = T5_get_summary_statistics.out3_or,
        appropriate_cancer_type = matched_cancer_prs[i],
        PGS_ID = PGS_ID
    }
    call T6_Combine_Values as T6_Combine_Pvalues_euro{
      input:
        files = T5_get_summary_statistics_euro.out2_p,
        appropriate_cancer_type = matched_cancer_prs[i],
        PGS_ID = PGS_ID + "-euro"
    }
    call T6_Combine_Values as T6_Combine_OR_euro{
      input:
        files = T5_get_summary_statistics_euro.out3_or,
        appropriate_cancer_type = matched_cancer_prs[i],
        PGS_ID = PGS_ID + "-euro"
    }
    call Tasks.concatenateFiles_noheader as concatenateFiles_noheader_Pvalues{
      input:
        files = [T6_Combine_Pvalues.out1,T6_Combine_Pvalues_euro.out1],
        callset_name = PGS_ID + ".pvalue"
    }
    call Tasks.concatenateFiles_noheader as concatenateFiles_noheader_OR {
      input:
        files = [T6_Combine_OR.out1,T6_Combine_OR_euro.out1],
        callset_name = PGS_ID + ".OR"
    }
  }
  call Tasks.concatenateFiles_noheader as concatenateFiles_noheader_Pvalues_meta{
    input:
      files = concatenateFiles_noheader_Pvalues.out2,
      callset_name =  "FINAL_PRS.pvalue"
  }
  call Tasks.concatenateFiles_noheader as concatenateFiles_noheader_OR_meta {
    input:
      files = concatenateFiles_noheader_OR.out2,
      callset_name = "FINAL_PRS.OR"
  }
  #call Tasks.concatenateFiles {
  #  input:
  #    files = T5_get_summary_statistics.out1,
  #    output_name = PGS_ID
  #}

  #call Tasks.concatenateFiles as concatenateFiles_euro{
  #  input:
  #    files = T5_get_summary_statistics_euro.out1,
  #    output_name = PGS_ID + "_euro"
  #}

  #call Tasks.copy_file_to_storage as copy3{
  #  input:
  #    text_file = concatenateFiles.out1,
  #    output_dir = analysis_4_output_dir + "ALL_BY_ALL/"
  #}

  #call Tasks.copy_file_to_storage as copy4{
  #  input:
  #    text_file = concatenateFiles_euro.out1,
  #    output_dir = analysis_4_output_dir + "ALL_BY_ALL/"
  #}
}

task T1_prepare_scores_file {
  input{
    File prs_file
  }
  command <<<
  python3 <<CODE
  import pandas as pd

  usecols = ['hm_source','hm_chr', 'hm_pos', 'effect_allele', 'other_allele', 'effect_weight']

  df = pd.read_csv("~{prs_file}", comment = "#", index_col=False, sep='\t',dtype={'hm_chr': str, 'hm_pos': str},usecols=usecols)
  df = df[df['hm_source'].notna() & (df['hm_source'].str.strip() != "")]

  # Keep only rows where hm_chr is in 1-22, X, or Y
  valid_chroms = set([str(i) for i in range(1, 23)] + ['X', 'Y'])
  df = df[df['hm_chr'].isin(valid_chroms)]

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
  score_df.to_csv("score_file.txt", sep="\t", index=False, header=False)
  CODE
  >>>
  output {
    File out1 = "score_file.txt"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
    memory: "8GB"
  }
}

task T4_control_for_ancestry {
  input {
    File raw_prs
    File sample_data
    String cancer_type
    String PGS_ID
    String appropriate_cancer_type
    String euro = "False"
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  import statsmodels.api as sm
  from scipy.stats import zscore

  # Load input files
  prs_df = pd.read_csv("~{raw_prs}", sep='\t', index_col=False)
  sample_data = pd.read_csv("~{sample_data}", sep='\t', index_col=False)

  subsets = {
    ("brain","nervous"),
    ("esophagus","gastrointestinal"),
    ("leukemia","hematologic"),
    ("pancreas","gastrointestinal")
  }

  if "~{appropriate_cancer_type}" == "~{cancer_type}" or ("~{appropriate_cancer_type}","~{cancer_type}") in subsets:
    pass
  else:
    sample_data = sample_data[~sample_data['original_dx'].str.contains("~{appropriate_cancer_type}", case = False)]

  # Limit to just Europeans for Sensitivity Analysis
  euro = "~{euro}" == "True"
  if euro:
    sample_data = sample_data[sample_data['intake_qc_pop'] == "EUR"]

  # Convert merge keys to string
  prs_df['sample'] = prs_df['sample'].astype(str)
  sample_data['original_id'] = sample_data['original_id'].astype(str)

  # Merge on sample ID
  merged_df = prs_df.merge(sample_data, left_on = "sample", right_on="original_id")

  merged_df[['PC1','PC2','PC3','PC4','age']] = merged_df[['PC1','PC2','PC3','PC4','age']].apply(zscore)
  merged_df['male_sex'] = merged_df['inferred_sex'].apply(lambda x: 1 if x == 'male' else 0)

  covariates = ['PC1', 'PC2', 'PC3', 'PC4']
  if merged_df['male_sex'].nunique() > 1:
    covariates.append('male_sex')
  # Set up regression: regress PGS ~ PC1 + PC2 + PC3 + PC4 + male_sex
  X = merged_df[covariates]
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

  # Choose filename
  suffix = ".euro" if euro else ""
  outname = f"~{cancer_type}.~{PGS_ID}{suffix}.pgs"
  outname2 = f"~{cancer_type}.~{PGS_ID}{suffix}.anon_pgs"
  # Write final output
  merged_df[['sample', 'PGS']].to_csv(outname, sep='\t', index=False)

  # Output Anonymized Files
  merged_df['case'] = (merged_df['cancer'] != 'control').astype(int)
  merged_df[['case', 'PGS']].to_csv(outname2, sep='\t', index=False)
  CODE
  >>>
  output {
    File out1 = glob("~{cancer_type}.~{PGS_ID}*.pgs")[0]
    File out2 = glob("~{cancer_type}.~{PGS_ID}*.anon_pgs")[0]
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
    memory: "8GB"
  }
}


task T5_get_summary_statistics {
  input {
    File adjusted_prs
    File sample_data
    String cancer_type
    String PGS_ID
  }
  command <<<
  python3 <<CODE
  from scipy.stats import ttest_ind
  from sklearn.metrics import roc_auc_score
  import statsmodels.api as sm
  import pandas as pd
  import numpy as np
  from scipy.stats import zscore

  f = open("~{cancer_type}.~{PGS_ID}.stats","w")

  df = pd.read_csv("~{adjusted_prs}",sep='\t',index_col=False)
  df['pgs_score'] = df['PGS']
  sample_data = pd.read_csv("~{sample_data}",sep='\t',index_col=False)
  df['sample'] = df['sample'].astype(str)
  sample_data['original_id'] = sample_data['original_id'].astype(str)
  merged_df = df.merge(sample_data, left_on = "sample", right_on = "original_id", how = "left")

  # 0. Normalize in PCs
  merged_df[['PC1','PC2','PC3','PC4','age']] = merged_df[['PC1','PC2','PC3','PC4','age']].apply(zscore)

  # 1. T-test between case and control PGS scores
  cases = merged_df[merged_df['cancer'] != 'control']['pgs_score']
  controls = merged_df[merged_df['cancer'] == 'control']['pgs_score']
  t_stat, p_val = ttest_ind(cases, controls, equal_var=False)

  f.write("\n\n####~{cancer_type}####\n\n".capitalize())
  f.write(f"T-test p-value (case vs control): {p_val:.4e}\n")

  # 2. Mean and std of the two distributions
  f.write(f"Case mean: {cases.mean():.4f}, std: {cases.std():.4f}\n")
  f.write(f"Control mean: {controls.mean():.4f}, std: {controls.std():.4f}\n")

  # 3. AUC of PRS
  # Binary labels: 1 = case, 0 = control
  y_true = merged_df['cancer'].apply(lambda x: 0 if x == 'control' else 1)
  auc = roc_auc_score(y_true, merged_df['pgs_score'])
  f.write(f"AUC of PGS score: {auc:.4f}\n")

  # 4. Logistic regression of cancer status on PGS with covariates
  # Prepare design matrix
  merged_df['sex_binary'] = merged_df['inferred_sex'].map({'male': 1, 'm': 1, 'female': 0, 'f': 0})

  # Base covariates
  covariates = ['pgs_score', 'PC1', 'PC2', 'PC3', 'PC4']

  # Add 'sex_binary' if it varies
  if merged_df['sex_binary'].nunique() > 1:
      covariates.append('sex_binary')

  Xy = merged_df[covariates + ['cancer']].dropna()

  # Now safely split cancer out to y
  y = Xy['cancer'].apply(lambda x: 0 if x == 'control' else 1)

  X = Xy[covariates]

  # Add constant for intercept
  X = sm.add_constant(X)

  # Fit logistic regression model
  logit_model = sm.Logit(y, X)
  result = logit_model.fit()

  # Get summary table as a DataFrame
  summary_df = result.summary2().tables[1]

  # Format p-values in scientific notation with higher precision
  summary_df['P>|z|'] = summary_df['P>|z|'].apply(lambda p: f"{p:.5e}" if not np.isnan(p) else "nan")

  # Write formatted table to file
  f.write(summary_df.to_string())

  f.close()

  coef_row = summary_df.loc['pgs_score']

  # Compute odds ratio and confidence interval
  or_value = np.exp(coef_row["Coef."])
  ci_low = np.exp(coef_row["[0.025"])
  ci_high = np.exp(coef_row["0.975]"])
  pval = coef_row["P>|z|"]

  # File 1: p-value label file
  pval_filename = f"~{cancer_type}.pvalue"
  with open(pval_filename, "w") as pf:
    pf.write(f"~{cancer_type}\n")
    pf.write(f"{pval}\n")

  # File 2: odds ratio label file
  or_filename = f"~{cancer_type}.or"
  with open(or_filename, "w") as of:
    of.write(f"~{cancer_type}\n")
    of.write(f"{or_value:.3f} ({ci_low:.3f}, {ci_high:.3f})\n")
  CODE
  >>>
  output {
    File out1 = "~{cancer_type}.~{PGS_ID}.stats"
    File out2_p = "~{cancer_type}.pvalue"
    File out3_or = "~{cancer_type}.or"
  }
  runtime {
    docker: "vanallenlab/pydata_stack"
    preemptible: 3
    memory: "8GB"
  }
}





task T6_Combine_Values {
  input {
    Array[File] files
    String PGS_ID
    String appropriate_cancer_type
  }

  command <<<
  python3 <<CODE
  # WDL inputs
  pvalue_files = "~{sep=' ' files}".split(" ")   # WDL will substitute the file paths
  output_file = "~{PGS_ID}.row"

  pairs = []
  for f in pvalue_files:
      with open(f) as fh:
          lines = [line.strip() for line in fh if line.strip()]
          if len(lines) != 2:
              raise ValueError(f"File {f} does not have exactly 2 lines")
          label, value = lines
          pairs.append((label, value))

  pairs.sort(key=lambda x: x[0])

  with open(output_file, "w") as out:
      labels = "\t".join(["Intended_Cancer_Type","PGS_ID"] + [label for label, _ in pairs])
      values = "\t".join(["~{appropriate_cancer_type}","~{PGS_ID}"] + [value for _, value in pairs])
      out.write(f"{labels}\n{values}\n")
  CODE
  >>>

  output {
    File out1 = "~{PGS_ID}.row"
  }

  runtime {
    docker: "vanallenlab/pydata_stack"
    memory: "1G"
    preemptible:3
  }
}
