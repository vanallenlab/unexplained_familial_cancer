# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0
workflow ANALYSIS_1D_LOCAL_GWAS {
  input {
    String cancer_type
    File vcf
    File vcf_idx
    String region
    String? output_dir
  }
  File metadata = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".metadata"
  File subjects_list = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".list"

  call T1_Convert_VCF {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      subjects_list = subjects_list,
      metadata = metadata,
      region = region,
      cancer_type = cancer_type
  }
  call T2_Perform_GWAS {
    input:
      metadata = T1_Convert_VCF.out1,
      cancer_type = cancer_type
  }
}

task T1_Convert_VCF {
  input {
    File vcf
    File vcf_idx
    File subjects_list
    File metadata
    String region
    String cancer_type
  }
  command <<<
  set -euxo pipefail
  
  # Filter VCF range
  bcftools view -r ~{region} ~{vcf} -Oz -o tmp.vcf.gz

  # Filter VCF to common variants
  bcftools view -i 'AN >= 100' tmp.vcf.gz -Oz -o tmp1.vcf.gz

  # Filter to samples of interest
  bcftools view -S ~{subjects_list} tmp1.vcf.gz -Oz -o out.vcf.gz

  # Convert to TSV w/ plink2
  plink2 --vcf out.vcf.gz --recode A --out gt_matrix.~{cancer_type}

  # Trim the file down
  cut -f2,7- gt_matrix.~{cancer_type}.raw > tmp.tsv

  python3 <<CODE
  import pandas as pd
  gt_df = pd.read_csv("tmp.tsv",sep='\t',index_col=False)
  metadata = pd.read_csv("~{metadata}",sep='\t',index_col=False)
  gt_df['original_id'] = gt_df['IID'].astype(str)
  metadata['original_id'] = metadata['original_id'].astype(str)

  metadata = metadata.merge(gt_df,how="left",on="original_id")
  metadata.to_csv("~{cancer_type}.local_gwas.tsv",sep='\t',index=False)
  CODE
  >>>
  output {
    File out1 = "~{cancer_type}.local_gwas.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible:3
  }
}

task T2_Perform_GWAS {
  input {
    File metadata
    String cancer_type
  }
  command <<<
  python3 <<CODE
  import pandas as pd
  import numpy as np
  import statsmodels.api as sm

  # Load
  df = pd.read_csv("~{metadata}", sep="\t")

  # Encode outcome + sex
  df["sex_binary"] = df["inferred_sex"].map({"male": 1, "female": 0})

  # Covariates
  covs = ["PC1", "PC2", "PC3", "PC4", "sex_binary"]

  # SNP columns
  snps = [c for c in df.columns if c.startswith("chr") and "ploidy" not in c.lower()]
  # Ensure SNPs are numeric where possible
  df[snps] = df[snps].apply(pd.to_numeric, errors="coerce").astype("Int8")

  models = {
      "additive": lambda x: x,
      "recessive": lambda x: (x == 2).astype(int),
      "dominant": lambda x: (x >= 1).astype(int),
      "homozygous": lambda x: (x != 1).astype(int)
  }

  results = []

  for snp in snps:
      d0 = df[["case_control"] + covs + [snp]].copy()

      # keep only valid genotypes
      d0 = d0[d0[snp].isin([0,1,2])].dropna()

      if d0.empty:
          continue

      for model_name, transform in models.items():
          d = d0.copy()
          d[snp] = transform(d[snp])

          if d[snp].nunique() < 2:
              continue

          X = sm.add_constant(d[covs + [snp]].astype(float))
          y = d["case_control"].astype(int)

          # remove any bad rows
          mask = X.notna().all(axis=1) & np.isfinite(X).all(axis=1) & y.notna()
          X = X[mask]
          y = y[mask]

          if len(X) < 20 or y.nunique() < 2:
              continue

          try:
              fit = sm.Logit(y, X).fit(disp=0)

              OR = np.exp(fit.params[snp])
              p = fit.pvalues[snp]
              ci_low, ci_high = np.exp(fit.conf_int().loc[snp])

              results.append({
                  "SNP": snp,
                  "model": model_name,
                  "OR": OR,
                  "CI_lower": ci_low,
                  "CI_upper": ci_high,
                  "p_value": p,
                  "N": len(X)
              })

          except Exception:
              continue

  # Save
  pd.DataFrame(results).to_csv("gwas_results.~{cancer_type}.tsv", sep="\t", index=False) #.sort_values("p_value").to_csv("gwas_results.~{cancer_type}.tsv", sep="\t", index=False)  
  CODE
  >>>
  output {
    File out1 = "gwas_results.~{cancer_type}.tsv"
  }
  runtime {
    docker:"vanallenlab/pydata_stack"
    preemptible:3
  }
}
  
