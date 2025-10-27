version 1.0

workflow ANALYSIS_6A_DO_IT_ALL {
  input {
    String cancer_type = "thyroid"
    String workspace_bucket = "fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228"

    # SAIGE_RESULTS
    File analysis_3b_saige_results

    # PRS metrics
    #Float analysis_4_p_cutoff = 0.05
    #File? prs_tsv

    # Damaging Missense
    #Float analysis_5_p_cutoff = 0.05
    #File step_10_cpg_file = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/STEP_10_VISUALIZE_VEP/v2/ufc.cpg.variant_counts.tsv.gz""
    #File cosmic_tsv = "gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/UFC_REFERENCE_FILES/cosmic_ufc.tsv"
    #File analysis5_results 
    #Array[String] tier_of_interest
    #Array[String] analysis_5_afs
    #File step10_results

    # ANALYSIS 1
    #File? homozygosity_tsv
    #Array[File]? analysis_1_genes

    # SV ANALYSIS?
    # File sv_tsv  # optional - commented out
    
  }

  # Metadata for cancer type
  File metadata_tsv = "gs://" + workspace_bucket + "/UFC_REFERENCE_FILES/analysis/" + cancer_type + "/" + cancer_type + ".metadata"

  call T1_normalize_saige {
    input: 
      saige_results_tsv=analysis_3b_saige_results
  }

  #call normalize_analysis5 {
  #  input: analysis5_tsv=analysis5_tsv
  #}

  #call normalize_step10 {
  #  input: step10_tsv=step10_tsv
  #}

  #call normalize_prs {
  #  input: prs_tsv=prs_tsv
  #}

  #call normalize_homozygosity {
  #  input: homoz_tsv=homozygosity_tsv,
  #         gene_of_interest=gene_of_interest
  #}

  # call normalize_svs {
  #   input: sv_tsv=sv_tsv
  # }

  call aggregate {
    input:
      saige_normalized=normalize_saige.saige_normalized,
      analysis5_normalized=normalize_analysis5.analysis5_normalized,
      step10_normalized=normalize_step10.step10_normalized,
      prs_normalized=normalize_prs.prs_normalized,
      metadata_normalized=normalize_metadata.metadata_normalized,
      homozygosity_normalized=normalize_homozygosity.homozygosity_normalized,
      gene_of_interest=gene_of_interest,
      cohort_name=cohort_name
  }

  output {
    File combined_tsv = aggregate.aggregated_tsv
  }
}

task T1_normalize_saige {
  input {
    File saige_results_tsv
  }

  command <<<
  set -euxo pipefail
  python3 <<'CODE'
  import pandas as pd

  # --- Load data ---
  df = pd.read_csv("~{saige_results_tsv}", sep="\t", index_col=False)

  # --- Basic cleanup ---
  df = df[['original_id', 'gene']].dropna().drop_duplicates()

  # --- Pivot to wide format ---
  wide = (
      df.assign(value=1)
        .pivot_table(index='original_id', columns='gene', values='value', fill_value=0)
        .reset_index()
  )

  # --- Save ---
  wide.to_csv("saige_normalized.tsv", sep="\t", index=False)
  CODE
  >>>

  output {
    File out1 = "saige_normalized.tsv"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
}


task normalize_analysis5 {
  input {
    File analysis5_tsv
  }

  command <<~PY
    set -euo pipefail
    python - <<'PYCODE'
import pandas as pd
df = pd.read_csv("~{analysis5_tsv}", sep="\t", dtype=str)
rows = []
# expected: damaging missense per-sample per-variant
for _, r in df.iterrows():
    sample = r.get("sample", r.get("sample_id", "")) or ""
    rows.append({
        "sample": sample,
        "gene": r.get("SYMBOL", r.get("gene","")),
        "variant": r.get("ID", r.get("variant","")),
        "score": r.get("REVEL", r.get("score","")),
        "meta": r.get("annotation",""),
        "source": "damaging_missense"
    })
out = pd.DataFrame(rows)
out.to_csv("analysis5_normalized.tsv", sep="\t", index=False)
PYCODE
  PY
  output {
    File analysis5_normalized = "analysis5_normalized.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    memory: "6G"
    cpu: 1
  }
}

task normalize_step10 {
  input {
    File step10_tsv
  }

  command <<~PY
    set -euo pipefail
    python - <<'PYCODE'
import pandas as pd
df = pd.read_csv("~{step10_tsv}", sep="\t", dtype=str)
rows = []
for _, r in df.iterrows():
    sample = r.get("sample", "") or r.get("sample_id", "")
    rows.append({
        "sample": sample or "",
        "gene": r.get("SYMBOL", r.get("gene","")),
        "variant": r.get("ID", r.get("variant","")),
        "score": r.get("score", r.get("impact_score","")),
        "meta": r.get("extra",""),
        "source": "step10"
    })
pd.DataFrame(rows).to_csv("step10_normalized.tsv", sep="\t", index=False)
PYCODE
  PY
  output {
    File step10_normalized = "step10_normalized.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    memory: "6G"
    cpu: 1
  }
}

task normalize_prs {
  input {
    File prs_tsv
  }

  command <<~PY
    set -euo pipefail
    python - <<'PYCODE'
import pandas as pd
prs = pd.read_csv("~{prs_tsv}", sep="\t", dtype=str)
# Expect columns: sample, prs_trait1, prs_trait2, ...
# Ensure sample column name detection
if "sample" not in prs.columns and "sample_id" in prs.columns:
    prs = prs.rename(columns={"sample_id":"sample"})
prs.to_csv("prs_normalized.tsv", sep="\t", index=False)
PYCODE
  PY
  output {
    File prs_normalized = "prs_normalized.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    memory: "4G"
    cpu: 1
  }
}

task normalize_metadata {
  input {
    File metadata_tsv
  }

  command <<~PY
    set -euo pipefail
    python - <<'PYCODE'
import pandas as pd
meta = pd.read_csv("~{metadata_tsv}", sep="\t", dtype=str)
# ensure sample column standardization
if "sample" not in meta.columns and "original_id" in meta.columns:
    meta = meta.rename(columns={"original_id":"sample"})
meta.to_csv("metadata_normalized.tsv", sep="\t", index=False)
PYCODE
  PY
  output {
    File metadata_normalized = "metadata_normalized.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    memory: "2G"
    cpu: 1
  }
}

task normalize_homozygosity {
  input {
    File homoz_tsv
    String gene_of_interest
  }

  command <<~PY
    set -euo pipefail
    python - <<'PYCODE'
import pandas as pd
df = pd.read_csv("~{homoz_tsv}", sep="\t", dtype=str)
# Expect columns: sample, gene, zygosity or similar
# We'll create a binary column indicating homozygosity for gene_of_interest
df_cols = {c.lower():c for c in df.columns}
sample_col = df_cols.get("sample", df_cols.get("sample_id", None))
gene_col = df_cols.get("gene", df_cols.get("symbol", None))
zyg_col = df_cols.get("zygosity", df_cols.get("genotype", df_cols.get("gt", None)))

if sample_col is None or gene_col is None:
    raise SystemExit("homoz file must contain sample and gene columns")

df['sample'] = df[sample_col].astype(str)
df['gene'] = df[gene_col].astype(str)
if zyg_col:
    df['zygosity'] = df[zyg_col].astype(str)
else:
    df['zygosity'] = ""

# Mark homozygous carriers for the gene_of_interest
gene = "~{gene_of_interest}"
homo = df[(df['gene'].str.lower() == gene.lower()) & df['zygosity'].str.contains('hom', case=False, na=False)]
homo_df = homo[['sample']].drop_duplicates()
homo_df['homo_'+gene] = 1
homo_df.to_csv("homozygosity_normalized.tsv", sep="\t", index=False)
PYCODE
  PY
  output {
    File homozygosity_normalized = "homozygosity_normalized.tsv"
  }
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    memory: "2G"
    cpu: 1
  }
}

# Structural variants normalization task - commented out per request
# task normalize_svs {
#   input {
#     File sv_tsv
#   }
#   command <<~PY
#     set -euo pipefail
#     python - <<'PYCODE'
# import pandas as pd
# df = pd.read_csv("~{sv_tsv}", sep="\t", dtype=str)
# # normalize to sample,gene,sv_type,info
# rows=[]
# for _,r in df.iterrows():
#     sample = r.get("sample","")
#     gene = r.get("gene","")
#     svtype = r.get("sv_type", r.get("type",""))
#     rows.append({"sample":sample,"gene":gene,"svtype":svtype,"info":r.to_json()})
# pd.DataFrame(rows).to_csv("svs_normalized.tsv", sep="\t", index=False)
# PYCODE
#   PY
#   output {
#     File svs_normalized = "svs_normalized.tsv"
#   }
#   runtime {
#     docker: "vanallenlab/g2c_pipeline"
#     memory: "8G"
#     cpu: 1
#   }
# }

task aggregate {
  input {
    File saige_normalized
    File analysis5_normalized
    File step10_normalized
    File prs_normalized
    File metadata_normalized
    File homozygosity_normalized
    # File svs_normalized  # optional - commented out for now
    String gene_of_interest
    String cohort_name = "cohort"
  }

  command <<~PY
    set -euo pipefail
    python - <<'PYCODE'
import pandas as pd
import sys
# Load normalized pieces
saige = pd.read_csv("~{saige_normalized}", sep="\t", dtype=str)
a5 = pd.read_csv("~{analysis5_normalized}", sep="\t", dtype=str)
s10 = pd.read_csv("~{step10_normalized}", sep="\t", dtype=str)
prs = pd.read_csv("~{prs_normalized}", sep="\t", dtype=str)
meta = pd.read_csv("~{metadata_normalized}", sep="\t", dtype=str)
homo = pd.read_csv("~{homozygosity_normalized}", sep="\t", dtype=str) if "__file__" else pd.DataFrame()

# make sure sample column exists
for df in (saige,a5,s10,prs,meta):
    if 'sample' not in df.columns:
        # try common alternatives
        for alt in ['Sample','sample_id','original_id']:
            if alt in df.columns:
                df.rename(columns={alt:'sample'}, inplace=True)
                break

# Start with set of all samples from metadata or prs
samples = set(meta['sample'].astype(str)) if 'sample' in meta.columns else set(prs['sample'].astype(str))
# also include any samples seen in variant tables
for tab in (saige,a5,s10):
    if 'sample' in tab.columns:
        samples.update(tab['sample'].astype(str).unique())

samples = sorted(s for s in samples if s and str(s).strip() != "")

# Create base dataframe
base = pd.DataFrame({"sample": samples})

# Add metadata columns (join)
if 'sample' in meta.columns:
    meta = meta.set_index('sample')
    base = base.join(meta, on='sample')

# Add PRS columns (join)
if 'sample' in prs.columns:
    prs = prs.set_index('sample')
    base = base.join(prs, on='sample', rsuffix='_prs')

# Aggregate variant burdens per gene for analysis files
def collapse_burden(df, prefix):
    if 'sample' not in df.columns:
        return pd.DataFrame(columns=['sample'])
    df2 = df.copy()
    # ensure numeric score if present
    if 'score' in df2.columns:
        df2['score'] = pd.to_numeric(df2['score'], errors='coerce')
    # set indicator of presence
    df2['present'] = 1
    agg = df2.groupby(['sample']).agg(
        **{f"{prefix}_variant_count": ('present','sum')}
    ).reset_index()
    return agg

agg_saige = collapse_burden(saige, 'saige')
agg_a5 = collapse_burden(a5, 'damaging_missense')
agg_s10 = collapse_burden(s10, 'step10')

# Join aggregates
for agg in (agg_saige, agg_a5, agg_s10):
    if not agg.empty:
        base = base.merge(agg, on='sample', how='left')

# Fill zeros for counts
for col in ['saige_variant_count','damaging_missense_variant_count','step10_variant_count']:
    if col in base.columns:
        base[col] = base[col].fillna(0).astype(int)
    else:
        base[col] = 0

# Add homozygosity indicator
if not homo.empty and 'sample' in homo.columns:
    homo = homo.set_index('sample')
    # column name like homo_{gene}
    homo_col = f"homo_{~{gene_of_interest}}"
    # but we can't have dynamic ~{} inside python - replace:
    homo_col = "homo_"+("~{gene_of_interest}")
    base = base.join(homo, on='sample')
    if homo_col not in base.columns:
        # if homo table used different name, try to coerce
        existing = [c for c in homo.columns if c.startswith('homo_')]
        if existing:
            base[existing[0]] = base[existing[0]].fillna(0).astype(int)
    else:
        base[homo_col] = base[homo_col].fillna(0).astype(int)

# optional: annotate presence of gene_of_interest damaging variants
g = "~{gene_of_interest}".lower()
def gene_presence(df, gene_name):
    if 'gene' not in df.columns:
        return pd.DataFrame(columns=['sample'])
    df2 = df[df['gene'].str.lower() == gene_name.lower()]
    if df2.empty:
        return pd.DataFrame(columns=['sample'])
    out = df2[['sample']].drop_duplicates()
    out[f"{gene_name}_carrier"] = 1
    return out

gp_saige = gene_presence(saige, g)
gp_a5 = gene_presence(a5, g)
gp_s10 = gene_presence(s10, g)

for gp, name in zip([gp_saige, gp_a5, gp_s10], ['saige','damaging_missense','step10']):
    if not gp.empty:
        base = base.merge(gp, on='sample', how='left')
        col = f"{g}_{name}_carrier"
        # if columns are named differently, standardize:
        # ensure column exists
        if col not in base.columns:
            existing = [c for c in base.columns if c.endswith('_carrier') and c.startswith(g)]
            if existing:
                base[col] = base[existing[0]].fillna(0).astype(int)
        base[col] = base.get(col, 0).fillna(0).astype(int)

# final cleaning: replace NaN with empty for string columns
base = base.fillna("")

# Write final TSV
base.to_csv("risk_factors_~{cohort_name}.tsv", sep="\t", index=False)
print("WROTE risk_factors_~{cohort_name}.tsv")
PYCODE
  PY
  output {
    File aggregated_tsv = "risk_factors_~{cohort_name}.tsv"
  }

  runtime {
    docker: "vanallenlab/g2c_pipeline"
    memory: "16G"
    cpu: 2
  }
}


